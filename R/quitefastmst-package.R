# This file is part of the quitefastmst package for R.

# ############################################################################ #
#                                                                              #
#   Copyleft (C) 2025-2025, Marek Gagolewski <https://www.gagolewski.com>      #
#                                                                              #
#                                                                              #
#   This program is free software: you can redistribute it and/or modify       #
#   it under the terms of the GNU Affero General Public License                #
#   Version 3, 19 November 2007, published by the Free Software Foundation.    #
#   This program is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               #
#   GNU Affero General Public License Version 3 for more details.              #
#   You should have received a copy of the License along with this program.    #
#   If this is not the case, refer to <https://www.gnu.org/licenses/>.         #
#                                                                              #
# ############################################################################ #



#' @title Euclidean and Mutual Reachability Minimum Spanning Trees
#'
#' @description
#' See \code{\link{mst_euclid}()} for more details.
#'
#' @details
#' For best speed, consider building the package from sources
#' using, e.g., \code{-O3 -march=native} compiler flags and with OpenMP
#' support on.
#'
#' @useDynLib quitefastmst, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @encoding UTF-8
#' @keywords internal
"_PACKAGE"
