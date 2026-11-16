# Summarise observed positivity at reported DRC survey coordinates.
# Run from the repository root. Each diagnostic and source dataset remains
# separate: datasets 1.1 and 2.1 explicitly overlap in the source documentation.

root <- "data-raw/external_blackflies"
out <- file.path(root, "maps")
source_file <- file.path(root, "drc_kakoi_koda_2026/tables/S1_epidemiological_data_utf8.csv")
records <- read.csv(source_file, check.names = FALSE, colClasses = "character",
                    na.strings = c("NA", ""), encoding = "UTF-8")

diagnostics <- c("Skin snip" = "microfilarial prevalence",
                 "Ov16 antibodies" = "anti-Ov16 RDT field",
                 "Nodules" = "REMO prevalence")
survey_names <- c(
  "1.1" = "2015: case-control, analytic dataset 1.1",
  "1.2 & 4.1" = "2015: separate case-control survey 1.2 / 4.1",
  "2.1" = "2015: mapping dataset 2.1",
  "5.1" = "2017: epilepsy cohort 5.1",
  "7.1" = "2016: community survey 7.1",
  "7.1.1" = "2018: child survey 7.1.1",
  "8.1" = "2021: child survey 8.1",
  "9.1" = "2017: adult male REA survey 9.1")

identifier <- paste(records[["Dataset order"]], records[["Row Identifier"]], sep = ":")
stopifnot(!anyDuplicated(identifier))
latitude <- suppressWarnings(as.numeric(records$Latitude))
longitude <- suppressWarnings(as.numeric(records$Longitude))
located <- is.finite(latitude) & is.finite(longitude) &
  abs(latitude) <= 90 & abs(longitude) <= 180

site_results <- list()
denominators <- list()
for (diagnostic in names(diagnostics)) {
  field <- diagnostics[[diagnostic]]
  tested <- !is.na(records[[field]]) & records[[field]] %in% c("0", "1")
  for (dataset_id in names(survey_names)) {
    belongs <- records[["Dataset order"]] == dataset_id
    denominator <- belongs & tested
    denominators[[length(denominators) + 1L]] <- data.frame(
      diagnostic = diagnostic, dataset_id = dataset_id,
      survey = survey_names[[dataset_id]],
      tested_all_coordinates = sum(denominator),
      positive_all_coordinates = sum(records[[field]][denominator] == "1"),
      tested_with_coordinates = sum(denominator & located),
      positive_with_coordinates = sum(records[[field]][denominator & located] == "1"),
      tested_without_coordinates = sum(denominator & !located))

    selected <- which(denominator & located)
    if (!length(selected)) next
    # No unweighted averaging of village percentages: each value is sum(y)/n.
    coordinate_key <- paste(latitude[selected], longitude[selected], sep = ":")
    for (indices in split(selected, coordinate_key)) {
      n_tested <- length(indices)
      n_positive <- sum(records[[field]][indices] == "1")
      sites <- unique(na.omit(records$Locality[indices]))
      design <- unique(na.omit(records[["Population-based or Case-control study"]][indices]))
      site_results[[length(site_results) + 1L]] <- data.frame(
        diagnostic = diagnostic, dataset_id = dataset_id,
        survey = survey_names[[dataset_id]],
        site = paste(sites, collapse = "; "),
        latitude = latitude[indices[1]], longitude = longitude[indices[1]],
        n_positive = n_positive, n_tested = n_tested,
        proportion_positive = n_positive / n_tested,
        survey_design = paste(design, collapse = "; "),
        year = paste(unique(na.omit(records[["Year of study implementation"]][indices])),
                     collapse = "; "),
        source_field = field, source_doi = "10.1371/journal.pntd.0014164",
        stringsAsFactors = FALSE)
    }
  }
}

site_results <- do.call(rbind, site_results)
denominators <- do.call(rbind, denominators)
stopifnot(all(site_results$n_positive <= site_results$n_tested),
          all(site_results$n_tested > 0),
          all(site_results$proportion_positive >= 0 & site_results$proportion_positive <= 1))
write.csv(site_results, file.path(out, "drc_infection_proportions.csv"),
           row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(denominators, file.path(out, "drc_infection_denominators.csv"),
           row.names = FALSE, na = "")
cat("Saved", nrow(site_results), "site / survey / diagnostic summaries.\n")
print(denominators[denominators$tested_with_coordinates > 0, ], row.names = FALSE)
