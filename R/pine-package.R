
#' The Pine Model Package
#' 
#' Full suite of all user-defined functions, subroutines, and data for the
#' high-elevation white pine/WPBR project.
#' 
#'
#' \tabular{ll}{
#' Package:  \tab pine \cr
#' Type:     \tab Package \cr
#' Version:  \tab 1.0 \cr
#' Date:     \tab 2017-01-07 \cr
#' License:  \tab GNU GPL 3 \cr
#' LazyLoad: \tab yes \cr
#' }
#' 
#' @name pine-package
#' @aliases pine-package pine
#' @docType package
#' @author Stu Field <wiki.kimi@gmail.com> \cr Maintainer: Stu Field
#' @seealso \code{\link{pine12}} \cr \code{\link{pine36}} \cr
#' Data Objects: \cr
#' *   \code{\link{pine12_solutions}} -- Generation solutions for pine12 \cr
#' *   \code{\link{pine36_solutions}} -- Generation solutions for pine36 \cr
#' @references This research is the result of a collaboration between the 
#' Colorado State University Biology Department, the Mathematics Department, 
#' and the United States Department of Agriculture, Forest Service: \cr
#'
#' Stu Field and Mike Antolin \cr
#' Department of Biology \cr
#' Colorado State University \cr
#' Fort Collins, CO  80523-1878 \cr
#'
#' Simon Tavener \cr
#' Department of Mathematics \cr
#' Colorado State University \cr
#' Fort Collins, CO  80523-1874 \cr
#'
#' Jen Klutsch and Anna Schoettle \cr
#' USDA Forest Service \cr
#' 240 West Prospect Road \cr
#' Fort Collins, CO  80526 \cr
#' 
#' @keywords whitebark pine package
#' @examples
#' 
#' pine12(Gen=10)
#' pine36(Gen=10)
#' 
NULL




#' Output Projections from pine12
#' 
#' All output including solutions for numerous time steps throughout the
#' projection. A list including the following time steps:
#' Gen2, Gen10, Gen100, Gen101, Gen1000, Gen1001
#' 
#' @name pine12_solutions
#' @aliases pine12_solutions
#' @docType data
#' @format List containing numerous objects from \code{\link{pine12}}
#' @references Stu Field
#' @source Stu Field
#' @keywords Pine12 Infection Model
NULL

#' Output Projections from pine36
#' 
#' All output including solutions for numerous time steps throughout the
#' projection. A list including the following time steps:
#' Gen2, Gen10, Gen100, Gen101, Gen1000, Gen1001
#' 
#' @name pine36_solutions
#' @aliases pine36_solutions
#' @docType data
#' @format List containing numerous objects from \code{\link{pine36}}
#' @references Stu Field
#' @source Stu Field
#' @keywords Pine36 Infection Model with Genetic Resistance
NULL

