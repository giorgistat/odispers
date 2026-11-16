# Draw the field diagnostics from display-thinned traces and precomputed ACFs.
# ACFs and numerical diagnostics must have been calculated before thinning.
draw_lgcp_field_mixing <- function(data, output_file) {
  colours <- c("#2F6596", "#207C82", "#C78A3B", "#963D32")
  png(output_file, width = 2600, height = 2250, res = 200, pointsize = 18)
  on.exit(dev.off())
  par(mfrow = c(3, 2), mar = c(4.1, 4.5, 3, 1), oma = c(1, 0, 2.8, 0), mgp = c(2.6, .8, 0))
  for (row in seq_along(data$variables)) {
    variable <- data$variables[row]
    diagnostic <- data$diagnostics[data$diagnostics$variable == variable, ]
    matplot(data$sweeps, data$traces[, , variable], type = "l", lty = 1, lwd = .6,
            col = adjustcolor(colours, .8), xlab = "Sweep (after burn-in)", ylab = "Field value",
            main = data$titles[row], cex.main = .92)
    mtext(sprintf("R-hat %.3f; bulk ESS %.1f / %s draws", diagnostic$rhat,
                  diagnostic$ess_bulk, format(data$retained_draws, big.mark = ",")),
          side = 3, line = .2, cex = .76)
    correlations <- data$autocorrelations[, , variable]
    matplot(data$lags, correlations, type = "l", lty = 1, lwd = 1.4, col = colours,
            xlab = "Lag (sweeps)", ylab = "Within-chain autocorrelation",
            ylim = range(-.1, 1, correlations), main = "Autocorrelation", cex.main = .92)
    abline(h = 0, col = "grey70", lty = 3)
  }
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = rep(0, 4))
  plot.new()
  legend("top", legend = paste0("Chain ", 1:4, " (initial width ", c(10, 20, 40, 80), " km)"),
         col = colours, lty = 1, lwd = 2, bty = "n", horiz = TRUE, cex = .8)
}
