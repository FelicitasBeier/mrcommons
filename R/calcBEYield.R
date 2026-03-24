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
                          cellular = TRUE, aggregate = FALSE)
    weight  <- setYears(dimSums(weight[, "y1995", ], dim = 3), "y2010")
    mapping <- toolGetMappingCoord2Country()
    mapping$coordiso <- paste(mapping$coords, mapping$iso, sep = ".")
    weightBegr <- weight
    weightBegr[is.na(out[, , "begr"])] <- 0
    weightBetr <- weight
    weightBetr[is.na(out[, , "betr"])] <- 0
    weightExp <- mbind(setNames(weightBegr, "begr"), setNames(weightBetr, "betr"))
    out[is.na(out)] <- 0
    out <- toolAggregate(out, rel = mapping, weight = weight,
                         from = "coordiso", to = "iso", dim = 1, zeroWeight = "setNA")
    out    <- toolCountryFill(out, fill = NA)
    weight <- toolAggregate(weight, rel = mapping, from = "coordiso", to = "iso", dim = 1)
    weight <- toolCountryFill(weight, fill = 0)
    weightOut <- mbind(setNames(weight, "betr"), setNames(weight, "begr"))
    weightOut[is.na(out)] <- 0
    out[is.na(out)]       <- 0 
  } else {
    weight <- calcOutput("Croparea", sectoral = "kcr", physical = TRUE,
                         cellular = TRUE, aggregate = FALSE)
    weight <- setYears(dimSums(weight[, "y1995", ], dim = 3), "y2010")
  }

  return(list(x            = out,
              weight       = weight,
              unit         = "t DM ha-1 yr-1",
              description  = "Bioenergy crop yields from Li et al. (2020)",
              isocountries = !cellular))
}
