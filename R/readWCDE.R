#' @title readWCDE
#' @description It reads the population projections by age, sex and education
#' of the Wittgenstein Centre Human Capital Data Explorer as downloaded by
#' [downloadWCDE()]. The output format follows [readLutz2014()]:
#' a magpie object with the dimensions country, year, scenario, sex, age and
#' education in million people. In addition to the education categories of
#' Lutz2014, the datasets contain a finer split of the highest levels
#' ("Short Post Secondary", "Bachelor" and "Master and higher"). These belong
#' to a separate classification level and do not sum to "Post Secondary"; in
#' "epop_v2" they are reported as NA. World and continental aggregates
#' (UN M49 codes >= 900) are removed.
#'
#' @param subtype "epop_v3" (projections 2020-2100, K.C. et al. 2024) or
#' "epop_v2" (SSP2 only, including the historical reconstruction 1950-2015,
#' Lutz et al. 2018)
#' @return magpie object with missing countries, which can be filled with
#' the function convertWCDE.
#'
#' @seealso
#' [downloadWCDE()], [convertWCDE()], [readLutz2014()]
#'
#' @importFrom reshape2 acast
readWCDE <- function(subtype = "epop_v3") {
  merge <- list()
  for (file in sort(list.files(pattern = "\\.rds$"))) {
    scenario <- sub("\\.rds$", "", file)
    d <- as.data.frame(readRDS(file), stringsAsFactors = FALSE)

    # remove world and continental/regional aggregates (UN M49 codes >= 900)
    d <- d[as.integer(as.character(d$country_code)) < 900, ]

    # change country names to iso codes; the Channel Islands (M49 code 830)
    # have no ISO 3166-1 code and are removed. Unlike Lutz2014, the WCDE
    # datasets contain real data for Taiwan, which is mapped to TWN.
    d$Area <- toolCountry2isocode(d$name, mapping = c("Taiwan Province of China" = "TWN"),
                                  ignoreCountries = "Channel Islands")
    d <- d[!is.na(d$Area), ]

    # add "y" in front of the years
    d$Year <- paste0("y", d$year)

    # transform into magpie object
    out <- acast(d, Area ~ Year ~ sex ~ age ~ education, value.var = "epop")
    out <- as.magpie(out)
    out <- add_dimension(out, dim = 3.1, add = "scenario", nm = scenario)
    merge[[length(merge) + 1]] <- out
  }
  merge <- do.call(mbind, merge)
  getSets(merge) <- c("country", "year", "scenario", "sex", "age", "education")
  merge <- merge / 1000
  return(merge)
}
