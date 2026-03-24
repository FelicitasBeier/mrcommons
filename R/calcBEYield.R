#' @title calcBEYield
#' @description Calculates bioenergy crop yields from Li et al. (2020) for the two
#' MAgPIE bioenergy types begr and betr. begr takes the maximum yield across
#' Miscanthus and Switchgrass per cell. betr takes the maximum yield across
#' Eucalypt, Poplar, and Willow per cell.
#'
#' @param cellular if TRUE returns 0.5deg grid level, if FALSE aggregates to country ISO level
#' @return List of magpie objects with bioenergy yields in t DM ha-1 yr-1
#' @author Kristine Karstens
#' @seealso [readLi2020()]
#' @examples
#' \dontrun{
#' calcOutput("BEYield")
#' calcOutput("BEYield", cellular = TRUE, aggregate = FALSE)
#' }
calcBEYield <- function(cellular = FALSE) {

  x <- readSource("Li2020", convert = FALSE)

  # begr: maximum yield across herbaceous bioenergy crops
  begr <- setNames(pmax(x[, , "Miscanthus"], x[, , "Switchgrass"], na.rm = TRUE), "begr")

  # betr: maximum yield across all woody bioenergy crops
  betr <- setNames(pmax(x[, , "Eucalypt"], x[, , "Poplar"], x[, , "Willow"], na.rm = TRUE), "betr")

  out <- mbind(begr, betr)

  if (!cellular) {
    weight  <- calcOutput("Croparea", sectoral = "kcr", physical = TRUE,
                          cellular = TRUE, cells = "lpjcell", aggregate = FALSE)
    weight  <- setYears(dimSums(weight[, "y1995", ], dim = 3), "y2010")
    mapping <- toolGetMappingCoord2Country()
    mapping$coordiso <- paste(mapping$coords, mapping$iso, sep = ".")
    out <- toolAggregate(out, rel = mapping, weight = weight + 10^-10,
                         from = "coordiso", to = "iso", dim = 1)
    out    <- toolCountryFill(out, fill = NA)
    weight <- toolAggregate(weight, rel = mapping, from = "coordiso", to = "iso", dim = 1)
    weight <- toolCountryFill(weight, fill = 0)
  } else {
    weight <- calcOutput("Croparea", sectoral = "kcr", physical = TRUE,
                         cellular = TRUE, cells = "lpjcell", aggregate = FALSE)
    weight <- setYears(dimSums(weight[, "y1995", ], dim = 3), "y2010")
  }

  return(list(x            = out,
              weight       = weight + 10^-10,
              unit         = "t DM ha-1 yr-1",
              description  = "Bioenergy crop yields from Li et al. (2020)",
              isocountries = !cellular))
}
