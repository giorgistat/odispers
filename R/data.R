#' REMO onchocerciasis nodule survey data
#'
#' Village-level results of the Rapid Epidemiological Mapping of
#' Onchocerciasis (REMO): the number of adults examined for palpable
#' onchocercal nodules and the number found positive, in 16 African
#' countries. The Sudan records (\code{country == "SUDAN"}) cover present-day
#' Sudan and South Sudan and are those analysed in the package scripts.
#'
#' @format A data frame with 14,126 rows and 6 columns:
#' \describe{
#'   \item{country}{Country name as recorded in REMO.}
#'   \item{village}{Village name.}
#'   \item{lat, long}{Latitude and longitude (decimal degrees, WGS 84).}
#'   \item{n}{Number of individuals examined.}
#'   \item{y}{Number with palpable nodules. A few values are not whole numbers
#'     in the source file; the analysis scripts round them.}
#' }
#' @source REMO data compiled by the African Programme for Onchocerciasis
#'   Control; see Zoure et al. (2014), Parasites & Vectors 7, 326.
#' @examples
#' data(remo)
#' sudan <- subset(remo, country == "SUDAN")
#' summary(sudan$y / sudan$n)
"remo"

#' River network of Sudan and South Sudan (HydroRIVERS)
#'
#' Reaches of HydroRIVERS v1.0 in the study area (longitude 21--40, latitude
#' 2.5--23.5) with a long-term average discharge of at least 0.1 m\eqn{^3}/s.
#'
#' @format An \code{sf} object with 93,094 LINESTRING features (WGS 84),
#'   each digitised from upstream to downstream, and columns:
#' \describe{
#'   \item{HYRIV_ID}{Reach identifier.}
#'   \item{NEXT_DOWN}{Identifier of the reach downstream (0 at an outlet);
#'     used by \code{\link{river_graph}} to join reaches at confluences.}
#'   \item{DIS_AV_CMS}{Modelled long-term average natural discharge (m\eqn{^3}/s).}
#' }
#' @source Lehner, B. and Grill, G. (2013). Global river hydrography and
#'   network routing: baseline data and new approaches to study the world's
#'   large river systems. Hydrological Processes 27, 2171--2186.
#'   \url{https://www.hydrosheds.org/products/hydrorivers}. Distributed under
#'   the HydroSHEDS licence, which requires this citation.
"sudan_rivers"

#' Aridity index for Sudan and South Sudan (Global Aridity Index v3.1)
#'
#' Annual aridity index, the ratio of mean annual precipitation to reference
#' evapotranspiration, at 30 arc-second resolution over the study area
#' (longitude 21--40, latitude 2.5--23.5). Values are stored as integers equal
#' to the index multiplied by 10,000, as in the original GeoTIFF.
#'
#' @format A packed \code{terra} raster (\code{PackedSpatRaster}); use
#'   \code{terra::unwrap(sudan_aridity)} to obtain a \code{SpatRaster}.
#'   \code{\link{extract_aridity}} accepts it directly and rescales the values.
#' @source Zomer, R. J., Xu, J. and Trabucco, A. (2022). Version 3 of the
#'   Global Aridity Index and Potential Evapotranspiration Database.
#'   Scientific Data 9, 409. Licence CC BY 4.0.
#' @examples
#' \dontrun{
#' data(sudan_aridity)
#' terra::plot(terra::unwrap(sudan_aridity) / 10000)
#' }
"sudan_aridity"

#' Borders of Sudan and South Sudan
#'
#' @format An \code{sf} object with two MULTIPOLYGON features (WGS 84) and a
#'   column \code{admin} with the country name.
#' @source Natural Earth, 1:50m admin-0 countries (public domain),
#'   \url{https://www.naturalearthdata.com}.
"sudan_borders"

#' RiverATLAS attributes of the Sudan and South Sudan river network
#'
#' Hydro-environmental attributes of the reaches in \code{\link{sudan_rivers}}
#' (same rows, same order), extracted from RiverATLAS, the river-reach layer of
#' HydroATLAS v1.0. They describe the larval habitat of the vector and are used
#' as breeding covariates: channel slope (a proxy for water speed), dry-season
#' discharge (flow permanence), inundation, lakes, regulation, forest cover and
#' lithology. Missing values are \code{NA}.
#'
#' @format A data frame with 93,094 rows and 15 columns:
#' \describe{
#'   \item{HYRIV_ID}{Reach identifier, as in \code{sudan_rivers}.}
#'   \item{sgr_dk_rav}{Channel slope of the reach (decimetres per km).}
#'   \item{dis_m3_pyr, dis_m3_pmn, dis_m3_pmx}{Natural discharge (m\eqn{^3}/s):
#'     annual mean, minimum month and maximum month.}
#'   \item{run_mm_cyr}{Land surface runoff in the reach catchment (mm per year).}
#'   \item{inu_pc_umn, inu_pc_umx, inu_pc_ult}{Inundated area in the upstream
#'     catchment (percent): minimum month, maximum month, long-term maximum.}
#'   \item{lka_pc_use}{Lake area in the upstream catchment (percent).}
#'   \item{dor_pc_pva}{Degree of regulation by dams (percent).}
#'   \item{for_pc_use}{Forest cover in the upstream catchment (percent).}
#'   \item{lit_cl_cmj}{Lithology class, majority in the reach catchment (HydroATLAS codes).}
#'   \item{ele_mt_cav}{Mean elevation of the reach catchment (m).}
#'   \item{riv_tc_usu}{River volume upstream (thousand m\eqn{^3}).}
#' }
#' @source Linke, S. et al. (2019). Global hydro-environmental sub-basin and
#'   river reach characteristics at high spatial resolution. Scientific Data 6,
#'   283. \url{https://www.hydrosheds.org/hydroatlas}. Distributed under the
#'   HydroSHEDS licence (CC BY 4.0), which requires this citation.
#' @examples
#' data(sudan_rivers); data(sudan_river_covariates)
#' stopifnot(identical(sudan_rivers$HYRIV_ID, sudan_river_covariates$HYRIV_ID))
#' summary(sudan_river_covariates$sgr_dk_rav)
"sudan_river_covariates"
