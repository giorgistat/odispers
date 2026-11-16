# Blackfly data for Sudan, South Sudan and neighboring countries

Downloaded on 2026-10-31. The collection contains 32 original files and 75 CSV
tables: 73 tables extracted from study files and two regional GBIF exports.
Original downloads are retained in each study's `originals/` folder. CSV
conversions are in `tables/`.

[Rivers and locations maps](maps/README.md) include an interactive HTML map,
regional PNG and DRC close-ups, with site coordinates and an omissions log.

## Start with these files

| Country / area | Data and period | Local file | Source and interpretation |
|---|---|---|---|
| Sudan: Atbara / Galabat | 40 breeding sites, coordinates, eggs/larvae/pupae in 2009, 2013 and 2015 | [Breeding-site table](sudan_atbara_2019/tables/breeding_sites_table_01.csv) | [Zarroug et al. 2019](https://doi.org/10.1186/s12879-019-4113-1). Ordinal abundance: `+` rare, `++` few, `+++` abundant. Original DOCX included. |
| Sudan: Abu Hamed and Galabat | Four villages, adult catches and aquatic-stage observations, 2007–2010 | [Study tables](sudan_seasonality_2016/tables/) | [Zarroug et al. 2016](https://doi.org/10.1371/journal.pone.0150309). Monthly series are presented as figures; the complete PDF is downloaded. No separate raw monthly-catch file was found. |
| Sudan: Merowe / Abu Hamed | Breeding-site presence and monthly adult counts before/after the dam, 2007/08 and 2010/11 | [Monthly counts and biting rates](sudan_merowe_2014/tables/article_table_02.csv) | [Zarroug et al. 2014](https://doi.org/10.1186/1756-3305-7-168). Breeding-site observations are in article table 1. |
| South Sudan: Naam and Yei | 16 surveyed sites with GPS coordinates and immature/adult abundance categories, 2023–2024 | [Site table](south_sudan_naam_yei_2025/tables/article_table_01.csv) | [Lakwo et al. 2025](https://doi.org/10.1371/journal.pntd.0013309). This is the 16-site table from the Results; the abstract reports 13 breeding sites. Monitoring series are figures in the downloaded PDF; DOCX supplements give statistical summaries. |
| South Sudan: Maridi | Childhood Ov16 surveys, five villages, 2019 and 2023 | [Original Excel data](south_sudan_maridi_2023/originals/Seroprevalence_dataset.xlsx) | [B2SHARE dataset](https://b2share.eudat.eu/records/bjjcy-s4c86), CC0. Original R script and CSV exports of all four sheets included. Human seroprevalence, not fly counts or adult nodule prevalence. |
| Ethiopia: Metema and Metekel | 11 sampling entries with coordinates, river, elevation, date and specimens, 2013 and 2017–2018 | [Site table](ethiopia_northwest_2024/tables/S2_sample_sites_table_01.csv) | [Hedtke et al. 2024](https://doi.org/10.1371/journal.pntd.0011868). Sample numbers are specimens used for cytotaxonomic/genetic study, not standardized abundance measurements. Three original DOCX supplements included. |
| Uganda: northern districts, with Tanzania comparison | Trap-level catches by date, trap type and species, 2015–2016; site coordinates and distance to breeding sites | [Original trap workbook](uganda_tanzania_traps_2017/originals/S1_trap_data.xlsx) and [CSV tables](uganda_tanzania_traps_2017/tables/) | [Hendy et al. 2017](https://doi.org/10.1371/journal.pntd.0005688). Includes Moyo sites 13–16 km from known Nile breeding sites in South Sudan. The workbook's Uganda, Tanzania and within-trap distribution sheets have different sampling units. |
| Uganda: western foci | 16 sentinel sites, human microfilarial prevalence, focus-level biting rates and intervention histories | [Study and supplement tables](uganda_western_2020/tables/) | [Michael et al. 2020](https://doi.org/10.1038/s41598-020-61194-w). Full DOCX supplement downloaded. The biting rates are not independent measurements at every village. |
| DRC and adjacent Uganda: Kakoi–Koda / Ituri | Coordinates, species, life-stage counts, human landing collections and crab trapping, 2009–2018 | [Original supplementary tables](drc_uganda_kakoi_koda_2022/originals/) and [CSV conversions](drc_uganda_kakoi_koda_2022/tables/) | [Post et al. 2022](https://doi.org/10.1371/journal.pntd.0010684). S2 gives collections by date/locality; S4 gives cross-border collections. Includes species other than the S. damnosum complex, notably S. dentulosum and S. vorax. |
| DRC: Kakoi–Koda / Ituri | Epidemiological records from several studies, with geographical and survey metadata | [UTF-8 CSV, 2176 records and 49 columns](drc_kakoi_koda_2026/tables/S1_epidemiological_data_utf8.csv) | [Amaral et al. 2026](https://doi.org/10.1371/journal.pntd.0014164). Original semicolon-delimited CSV and appendix downloaded. The current S1 Data contains epidemiological data; published entomological data are referenced separately. Mixed study designs and diagnostics must be separated before analysis. |
| Chad: Oulibangala and Kaitia | Monthly catches at two sites on the Lim/Logone rivers, June 2021–February 2022 | [Monthly catch table](chad_oulibangala_2026/tables/article_table_01.csv) | [Oumar et al. 2026](https://doi.org/10.1093/inthealth/ihaf141). Counts plus published infection-screening results in the source XML. |
| Regional GBIF: Kenya, Uganda, Sudan and CAR | 160 georeferenced Simuliidae occurrence records | [Occurrence records](gbif_regional_occurrences_2026/tables/georeferenced_simuliidae.csv) and [country counts](gbif_regional_occurrences_2026/tables/country_counts.csv) | [GBIF occurrence search](https://www.gbif.org/occurrence/search?taxon_key=3522). Kenya 151, Uganda 4, Sudan 3, CAR 2. Museum, observation and molecular records; many identified only to family or molecular cluster. These are presence records, not standardized catch counts or confirmed vector observations. Record licenses, source dataset keys and GBIF issue flags are retained. |

## Sources found but not downloaded

- **Kenya: Gura and Mara rivers.** [Dryad 10.5061/dryad.k3j9kd5nc](https://doi.org/10.5061/dryad.k3j9kd5nc)
  describes abundance and habitat data from seven Gura and nine Mara sites.
  Dataset metadata and the 13-file listing were downloaded. The public file
  server returned HTTP 403, including a retry outside the sandbox; API file
  downloads returned HTTP 401. The 13 data/code files are marked `blocked` in
  the download manifest. The Kenyan GBIF records above are a separate source.
- **Central African Republic: Boali and Zinga.** [Enyong et al. 2006](https://doi.org/10.1051/parasite/2006131035)
  reports baseline catches, biting and transmission indices. The publisher's
  PDF endpoint returned HTTP 403. The paper is recorded as `blocked`.
- The regional GBIF query also included South Sudan, Ethiopia, DRC, Chad,
  Eritrea, Egypt and Libya, and returned no georeferenced records for those
  countries under the selected taxon key. This does not establish absence of
  blackflies or absence of other datasets.

## Provenance and use

- [downloads.csv](downloads.csv) records every attempted file, its URL, status,
  byte size and SHA-256 checksum. A metadata record is not counted as a
  downloaded biological dataset.
- [tables.json](tables.json) links each extracted table to its source and
  retains captions, footnotes or document notes. Keep these alongside the CSVs.
- CSV conversions preserve source categories and missing markers. Multiple
  header rows, section labels and source spelling remain. Vertically merged
  cells are repeated; horizontally merged labels occupy the first cell.
- The DRC epidemiological CSV was converted from Windows-1252 and semicolons
  to UTF-8 and commas, preserving embedded newlines and all 49 fields.
- PDF table conversions retain life-stage strings such as `5 P + 13 L` and
  geographical strings. Page-specific CSVs require their original species
  headers. PDF extraction was spot-checked against a rendered source page;
  these are source transcriptions, not harmonized analysis datasets.
- Check coordinates against maps before matching rivers. Apparent source
  coordinate errors have not been silently corrected. Align survey dates with
  dam construction, vegetation clearing, larviciding and ivermectin treatment.
- Source licenses remain in original article XML or repository metadata. The
  GBIF export preserves individual record licenses. Cite the underlying study
  and dataset when reusing data.

To repeat the document/spreadsheet conversions, run `extract_tables.py` with
Python, `openpyxl` and `pdfplumber`. This regenerates the 73 extracted CSV
tables. The two GBIF CSVs were exported separately from the retained JSON.
