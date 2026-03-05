#' @title calcNitrogenFixationFreeliving
#' @description calculates fixation from freeliving bacteria
#' @param cellular cellular estimates optional
#' @param irrigation if TRUE, distinguishes irrigated and non-irrigated crops
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Benjamin Leon Bodirsky
#' @examples
#' \dontrun{
#' calcOutput("NitrogenFixationFreeliving")
#' }
#' @importFrom magpiesets findset
#'

calcNitrogenFixationFreeliving <- function(cellular = FALSE, irrigation = FALSE) {

  area <- collapseNames(calcOutput("Croparea", cellular = cellular, aggregate = FALSE,
                                   physical = TRUE, irrigation = irrigation))
  fallow <- calcOutput("FallowLand", aggregate = FALSE)
  commonYears <- intersect(getYears(area), getYears(fallow))
  area <- mbind(area[, commonYears, ], fallow[, commonYears, ])
  freeliving <- setYears(readSource("Herridge", subtype = "freeliving", convert = FALSE), NULL)
  freeliving <- mbind(
    freeliving,
    setNames(freeliving[, , "tece"], "fallow") # use value of temperate cereals for fallow
  )
  freeliving <- area * freeliving
  freeliving <- dimSums(freeliving, dim = 3)
  fixFreeliving <- setNames(freeliving, "fixation_freeliving")

  return(list(x = fixFreeliving,
              weight = NULL,
              unit = "Mt Nr",
              description = "Nitrogen fixation by freeliving microorganisms",
              isocountries = !cellular))
}
