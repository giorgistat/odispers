#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

// Binomial log likelihood with log exposure as an offset, excluding the constant.
static double survey_loglik(double exposure, double predictor, double y, double n, bool logit) {
  // Inactive and held-out villages have n = y = 0. They contribute no
  // binomial term; the caller still evaluates their activity likelihood.
  if (n == 0) return 0.0;
  if (exposure == 0) return (y > 0) ? R_NegInf : 0.0;
  double log_force = std::log(exposure) + predictor;
  if (logit) {
    double tail = std::log1p(std::exp(-std::abs(log_force)));
    double ll = -n * tail;
    if (y > 0 && log_force < 0) ll += y * log_force;
    if (n > y && log_force > 0) ll -= (n-y) * log_force;
    return ll;
  }
  double force = std::exp(log_force);
  double ll = (y < n) ? -(n - y) * force : 0.0;
  if (y > 0) {
    double log_p = log_force < std::log(std::numeric_limits<double>::epsilon()) ? log_force :
      (force < std::log(2.0) ? std::log(-std::expm1(-force)) : std::log1p(-std::exp(-force)));
    ll += y * log_p;
  }
  return ll;
}

// One balanced sum tree per village. Each nonzero of K is one permanent
// leaf. Replacing a site's mass requires only positive pairwise additions
// along its leaf-to-root path, even when the old exposure nearly cancels.
// [[Rcpp::export]]
List lgcp_tree_layout(const S4& K) {
  IntegerVector village = K.slot("i"), dims = K.slot("Dim");
  const int m = dims[0];
  IntegerVector degree(m), size(m), base(m), used(m), leaf(village.size());
  for (int i : village) ++degree[i];
  int total = 0;
  for (int i = 0; i < m; ++i) {
    int capacity = 1;
    while (capacity < degree[i]) {
      if (capacity > INT_MAX/4) stop("LGCP village kernel tree is too large");
      capacity *= 2;
    }
    if (total > INT_MAX - 2*capacity) stop("LGCP kernel layout exceeds integer range");
    base[i] = total;
    size[i] = capacity;
    total += 2*capacity;
  }
  for (int k = 0; k < village.size(); ++k) {
    int i = village[k];
    leaf[k] = size[i] + used[i]++;
  }
  return List::create(_["base"]=base, _["size"]=size, _["leaf"]=leaf, _["total"]=total);
}

// CSC columns identify villages touched by a proposed cell update.
// [[Rcpp::export]]
List lgcp_local_sweep(const S4& K, const List& layout, const IntegerVector& old_counts,
                      const IntegerVector& old_high, const NumericVector& mu,
                      const NumericVector& predictor, const NumericVector& y,
                      const NumericVector& n, const IntegerVector& order,
                      int block_size, double spread, const IntegerVector& relocation_group,
                      bool logit, const NumericVector& gate_base,
                      const IntegerVector& activity_sign, double gate_slope) {
  IntegerVector counts = clone(old_counts), high = clone(old_high);
  const bool marked = high.size() != 0;
  IntegerVector col = K.slot("p"), village = K.slot("i");
  IntegerVector base = layout["base"], size = layout["size"], leaf = layout["leaf"];
  NumericVector weight = K.slot("x");
  std::vector<double> tree(as<int>(layout["total"]), 0.0);
  const int m = y.size(), cells = counts.size();
  if (mu.size() != cells || order.size() != cells || (marked && high.size() != cells) ||
      predictor.size() != m || n.size() != m || base.size() != m || size.size() != m || leaf.size() != weight.size() ||
      col.size() != cells + 1 || block_size < 1 || !(spread >= 0 && spread <= 1))
    stop("incompatible dimensions or controls in LGCP count update");
  if (activity_sign.size() != gate_base.size() || (gate_base.size() && gate_base.size() != m))
    stop("incompatible occurrence predictors in LGCP count update");
  auto village_loglik = [&](int i, double E) {
    double ll = survey_loglik(E, predictor[i], y[i], n[i], logit);
    if (activity_sign.size() && activity_sign[i]) {
      double t = activity_sign[i] * (gate_base[i] + gate_slope*E);
      ll += std::min(0.0,t) - std::log1p(std::exp(-std::abs(t)));
    }
    return ll;
  };

  auto cell_mass = [&](int j) {
    return marked ? (counts[j]-high[j])*(1-spread) + high[j]*(1+spread) : double(counts[j]);
  };
  // Build the current exposures from nonnegative contributions. Subsequent
  // proposals replace leaves instead of subtracting nearly equal totals.
  for (int j = 0; j < cells; ++j) {
    const double mass = cell_mass(j);
    if (mass == 0) continue;
    for (int k = col[j]; k < col[j+1]; ++k)
      tree[base[village[k]] + leaf[k]] = weight[k]*mass;
  }
  NumericVector exposure(m), loglik(m), trial_exposure(m), trial_loglik(m);
  for (int i = 0; i < m; ++i) {
    double* sums = tree.data() + base[i];
    for (int k = size[i]-1; k >= 1; --k) sums[k] = sums[2*k] + sums[2*k+1];
    exposure[i] = sums[1];
    loglik[i] = village_loglik(i, exposure[i]);
    if (!R_FINITE(loglik[i])) stop("non-finite initial survey likelihood in LGCP count update");
  }
  std::vector<int> seen(m, -1), touched, changed, previous_counts, previous_high;
  auto replace_cell = [&](int j, int tag) {
    const double mass = cell_mass(j);
    for (int t = col[j]; t < col[j+1]; ++t) {
      int i = village[t];
      double* sums = tree.data() + base[i];
      int k = leaf[t];
      sums[k] = weight[t]*mass;
      while (k > 1) {
        k /= 2;
        sums[k] = sums[2*k] + sums[2*k+1];
      }
      if (tag >= 0 && seen[i] != tag) {
        seen[i] = tag;
        touched.push_back(i);
      }
    }
  };
  auto accept_changes = [&]() {
    double diff = 0;
    for (int i : touched) {
      trial_exposure[i] = tree[base[i] + 1];
      trial_loglik[i] = village_loglik(i, trial_exposure[i]);
      diff += trial_loglik[i] - loglik[i];
    }
    if (!(std::log(R::runif(0, 1)) < diff)) return false;
    for (int i : touched) {
      exposure[i] = trial_exposure[i];
      loglik[i] = trial_loglik[i];
    }
    return true;
  };
  int proposed = 0, accepted = 0;
  // Refresh counts from their conditional Poisson prior (and marks from
  // their binomial prior), so only the village likelihood enters acceptance.
  for (int start = 0; start < cells; start += block_size) {
    touched.clear(); changed.clear(); previous_counts.clear(); previous_high.clear();
    for (int k = start; k < start + std::min(block_size, cells-start); ++k) {
      int j = order[k] - 1;
      if (j < 0 || j >= cells) stop("invalid cell index in LGCP count update");
      double draw = R::rpois(mu[j]);
      if (!R_FINITE(draw) || draw > INT_MAX) stop("Poisson site count exceeds integer range");
      int new_count = static_cast<int>(draw);
      int new_high = marked ? static_cast<int>(R::rbinom(new_count, 0.5)) : 0;
      if (counts[j] == new_count && (!marked || high[j] == new_high)) continue;
      changed.push_back(j); previous_counts.push_back(counts[j]);
      counts[j] = new_count;
      if (marked) { previous_high.push_back(high[j]); high[j] = new_high; }
      replace_cell(j, start);
    }
    if (changed.empty()) continue;
    ++proposed;
    if (accept_changes()) {
      ++accepted;
    } else {
      for (size_t k = 0; k < changed.size(); ++k) {
        counts[changed[k]] = previous_counts[k];
        if (marked) high[changed[k]] = previous_high[k];
        replace_cell(changed[k], -1);
      }
    }
  }
  int moved = 0, move_proposals = 0;
  if (relocation_group.size()) {
    // Conditional on the type counts within a tile, a labelled site's
    // position has prior probabilities proportional to the Poisson means.
    // Sampling its new cell from these probabilities cancels the prior
    // exactly. Repeated random label selection preserves exchangeability.
    if (relocation_group.size() != cells) stop("invalid relocation partition");
    int groups = max(relocation_group);
    if (groups < 1 || groups > cells) stop("invalid relocation group index");
    std::vector<std::vector<int>> members(groups);
    std::vector<std::vector<double>> cumulative(groups);
    std::vector<int> locations, is_high;
    for (int j = 0; j < cells; ++j) {
      int g = relocation_group[j]-1;
      if (g < 0 || g >= groups) stop("invalid relocation group index");
      members[g].push_back(j);
      double previous = cumulative[g].empty() ? 0.0 : cumulative[g].back();
      cumulative[g].push_back(previous + mu[j]);
      for (int k = 0; k < counts[j]; ++k) {
        locations.push_back(j);
        is_high.push_back(marked && k < high[j]);
      }
    }
    std::fill(seen.begin(), seen.end(), -1);
    for (int step = 0; step < static_cast<int>(locations.size()); ++step) {
      int label = static_cast<int>(R::runif(0,1) * locations.size());
      int old = locations[label], g = relocation_group[old]-1;
      const auto& cdf = cumulative[g];
      double draw = R::runif(0,1) * cdf.back();
      size_t pos = std::upper_bound(cdf.begin(), cdf.end(), draw) - cdf.begin();
      int next = members[g][std::min(pos, cdf.size()-1)];
      if (old == next) continue;
      ++move_proposals;
      --counts[old]; ++counts[next];
      if (is_high[label]) { --high[old]; ++high[next]; }
      touched.clear();
      replace_cell(old, step);
      replace_cell(next, step);
      if (accept_changes()) {
        ++moved;
        locations[label] = next;
      } else {
        ++counts[old]; --counts[next];
        if (is_high[label]) { ++high[old]; --high[next]; }
        replace_cell(old, -1);
        replace_cell(next, -1);
      }
    }
  }
  return List::create(_["counts"]=counts, _["high"]=high, _["exposure"]=exposure,
                      _["proposals"]=proposed, _["accepts"]=accepted,
                      _["move_proposals"]=move_proposals, _["moves"]=moved);
}

// Split the undirected metric graph into alternating breadth-first bands.
// A connected block of one colour has no precision edge to another block
// of that colour. Moving the band boundaries between scans avoids a fixed cut.
// [[Rcpp::export]]
List lgcp_graph_blocks(const IntegerVector& from, const IntegerVector& to, int vertices,
                       int width, int shift) {
  if (vertices <= 0 || width <= 0 || shift < 0 || from.size() != to.size())
    stop("invalid graph dimensions or band width");
  std::vector<std::vector<int>> neighbours(vertices);
  for (int k = 0; k < from.size(); ++k) {
    int i = from[k]-1, j = to[k]-1;
    if (i < 0 || i >= vertices || j < 0 || j >= vertices) stop("invalid graph vertex");
    neighbours[i].push_back(j); neighbours[j].push_back(i);
  }
  std::vector<int> depth(vertices,-1), queue;
  IntegerVector colour(vertices), group(vertices, -1);
  for (int root = 0; root < vertices; ++root) {
    if (depth[root] >= 0) continue;
    queue.clear(); queue.push_back(root); depth[root] = 0;
    for (size_t k = 0; k < queue.size(); ++k) {
      int i = queue[k];
      colour[i] = ((depth[i] + shift) / width) % 2;
      for (int j : neighbours[i]) if (depth[j] < 0) {
        depth[j] = depth[i]+1; queue.push_back(j);
      }
    }
  }
  int label = 0;
  for (int root = 0; root < vertices; ++root) {
    if (group[root] >= 0) continue;
    queue.clear(); queue.push_back(root); group[root] = ++label;
    for (size_t k = 0; k < queue.size(); ++k) {
      int i = queue[k];
      for (int j : neighbours[i]) if (group[j] < 0 && colour[i] == colour[j]) {
        group[j] = label; queue.push_back(j);
      }
    }
  }
  return List::create(_["colour"]=colour, _["group"]=group);
}

// Independent elliptical slices conditional on the other colour. Input
// predictors have already been projected onto the affected river cells.
// [[Rcpp::export]]
NumericMatrix lgcp_block_angles(const IntegerVector& group, int groups,
                                const NumericVector& offset, const NumericVector& current,
                                const NumericVector& direction, const IntegerVector& counts) {
  if (groups < 1 || current.size() != group.size() || offset.size() != group.size() ||
      direction.size() != group.size() || counts.size() != group.size())
    stop("incompatible block-slice dimensions");
  std::vector<std::vector<int>> members(groups);
  for (int j = 0; j < group.size(); ++j) {
    if (group[j] < 1 || group[j] > groups) stop("invalid conditional block index");
    members[group[j]-1].push_back(j);
  }
  NumericMatrix rotate(groups, 2);
  const double circle = 2 * std::acos(-1.0);
  for (int g = 0; g < groups; ++g) {
    if (members[g].empty()) { rotate(g,1) = 1; continue; }  // fresh conditional prior draw
    auto likelihood = [&](double c, double s) {
      double value = 0;
      for (int j : members[g]) {
        double effect = c * current[j] + s * direction[j];
        value += counts[j] * effect - std::exp(offset[j] + effect);
      }
      return value;
    };
    double value = likelihood(1,0);
    if (!R_FINITE(value)) stop("non-finite Poisson likelihood at current graph block");
    double threshold = value + std::log(R::runif(0,1));
    double angle = R::runif(0,circle), lo = angle-circle, hi = angle;
    for (;;) {
      double c = std::cos(angle), s = std::sin(angle);
      if (likelihood(c,s) >= threshold) { rotate(g,0) = c; rotate(g,1) = s; break; }
      if (angle < 0) lo = angle; else hi = angle;
      angle = R::runif(lo,hi);
    }
  }
  return rotate;
}
