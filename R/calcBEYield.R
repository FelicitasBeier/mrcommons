#' @title calcBEYield
#' @description Calculates bioenergy crop yields from Li et al. (2020) for the two
#' MAgPIE bioenergy types begr and betr. begr takes the maximum yield across
#' Miscanthus and Switchgrass per cell. betr takes the maximum yield across
#' Eucalypt, Poplar, and Willow per cell.
#'
#' @param returnWeights if TRUE returns crop-specific cellular cropland area weights
#'   (begr/betr) instead of the yields
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Kristine Karstens
#' @seealso [readLi2020()]
#' @examples
#' \dontrun{
#' calcOutput("BEYield", aggregate = FALSE)
#' calcOutput("BEYield", returnWeights = TRUE, aggregate = FALSE)
#' }
calcBEYield <- function(returnWeights = FALSE) {

  x <- setYears(readSource("Li2020", convert = FALSE), NULL)

  # begr: cell-wise maximum across herbaceous crops; na.rm = TRUE retains cells
  # where only one of the two species has a valid yield
  begr <- setNames(pmax(x[, , "Miscanthus"], x[, , "Switchgrass"], na.rm = TRUE), "begr")

  # betr: cell-wise maximum across woody crops; same na.rm rationale
  betr <- setNames(pmax(x[, , "Eucalypt"], x[, , "Poplar"], x[, , "Willow"], na.rm = TRUE), "betr")

  out <- mbind(begr, betr)

  # Cropland area (summed over all crops, 1995 base year) used as aggregation weight.
  weight <- calcOutput("Croparea", sectoral = "kcr", physical = TRUE, irrigation = TRUE,
                       cellular = TRUE, aggregate = FALSE)
  weightDetail <- setYears(dimSums(weight[, "y1995", ], dim = "kcr"), NULL)
  weight       <- setYears(dimSums(weight[, "y1995", ], dim = 3), NULL)

  # Crop-specific weights: cells with NA yield are set to 0 so they are
  # excluded from weighted aggregation and flagged via zeroWeight = "setNA".
  weightBegr <- weight
  weightBegr[is.na(out[, , "begr"])] <- 0
  weightDetailBegr <- weightDetail
  weightDetailBegr[is.na(out[, , "begr"])] <- 0
  weightBetr <- weight
  weightBetr[is.na(out[, , "betr"])] <- 0
  weightDetailBetr <- weightDetail
  weightDetailBetr[is.na(out[, , "betr"])] <- 0
  weightExp <- mbind(setNames(weightBegr, "begr"), setNames(weightBetr, "betr"))
  weightDetailExp <- mbind(add_dimension(weightDetailBegr, dim = 3.1, add = "kbe", nm = "begr"),
			   add_dimension(weightDetailBetr, dim = 3.1, add = "kbe", nm = "betr"))
  out[is.na(out)] <- 0

  if (returnWeights) {
    return(list(x            = weightDetailExp,
                weight       = NULL,
                unit         = "Mha",
                description  = paste("Cellular cropland area weights (y1995 cropland area,",
                                     "timeless) used for aggregation of bioenergy crop",
                                     "yields from Li et al. (2020)"),
                isocountries = FALSE))
  }

  return(list(x            = out,
              weight       = weightExp,
              unit         = "t DM ha-1 yr-1",
              description  = "Bioenergy crop yields from Li et al. (2020)",
              isocountries = FALSE))
}
