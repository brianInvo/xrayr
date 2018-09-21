
#' Error measurements of images (from redr)
#'
#' This function calculates error between two images
#'
#' @param x,y \code{antsImage} objects
#'
#' @examples
#' library( ANTsRCore )
#' MSE( ri(1), ri(2) )
#' #alternatively it can be done like:
#' MSE( ri(1) - ri(2) )
#' @name error
NULL

#' @export
#' @describeIn error Mean Squared Error
MSE <- function(x, y = NULL){

  if(is.null(y))
    return(mean(x^2))
  else
    return(mean((x - y)^2))
}

#' @export
#' @describeIn error Mean Absolute Error
MAE <- function(x, y = NULL){

  if(is.null(y))
    return(mean(abs(x)))
  else
    return(mean(abs(x - y)))
}

#' @export
#' @describeIn error Peak Signal-to-Noise Ratio
PSNR <- function(x, y){
  20*log10(max(x)) - 10*log10(MSE(x, y))
}



#' Suggested parameters for deformable intrasubject xray registration.
#'
#' This function, by default, regularizes the total deformation field to produce a
#' a constrained diffeomorphic mapping.  Full custom parameters are possible
#' by passing the same parameters that one would pass to \code{antsRegistration}.
#'
#' @param fixed fixed target image
#' @param moving image to map to the fixed space
#' @param totalSigma regularization for total deformation field
#' @param flowSigma regularization for update deformation field
#' @param regIterations deformable registration iterations / multi-resolution parameters
#' @param ... parameters to pass to antsRegistration
#' @return deformed image
#' @author Avants BB
#' @examples
#'
#' library( magrittr )
#' i1 = ANTsRCore::ri( 1 ) %>% ANTsRCore::resampleImage( 4 )
#' i2 = ANTsRCore::ri( 2 ) %>% ANTsRCore::resampleImage( 4 )
#' reg = xrayMotionCorrect( i1, i2 )
#'
#' @export xrayMotionCorrect
#' @importFrom stats rnorm
#' @importFrom magrittr %>%
#' @importFrom ANTsRCore antsRegistration antsImageRead
#' @importFrom ANTsRCore resampleImage antsAverageImages
#' @importFrom ANTsRCore ri iMath
#' @importFrom ANTsR mergeListToNDImage splitNDImageToList
#' @importFrom pracma detrend
xrayMotionCorrect <- function( 
  fixed, 
  moving, 
  totalSigma, 
  flowSigma, 
  regIterations = c( 50, 50, 10 ), 
  ... ) {
  if ( missing( totalSigma ) ) totalSigma = 6
  if ( missing( flowSigma ) ) flowSigma = 6
  reg = ANTsRCore::antsRegistration( fixed, moving, 'SyN', 
    regIterations = regIterations, 
    flowSigma = flowSigma, totalSigma = totalSigma, ... )
  return( reg$warpedmovout )
}




#' Suggested parameters for deformable motion correction for a time series.
#'
#' This function, by default, regularizes the total deformation field to produce a
#' a constrained diffeomorphic motion correction of a time series.  Full custom
#' by passing the same parameters that one would pass to \code{antsRegistration}.
#'
#' @param fixedStaticImage fixed static image
#' @param movingTimeSeries time series image to map to the fixed space
#' @param totalSigma regularization for total deformation field
#' @param flowSigma regularization for update deformation field
#' @param regIterations deformable registration iterations / multi-resolution parameters
#' @param ... parameters to pass to antsRegistration
#' @return deformed images
#' @author Avants BB
#' @examples
#'
#' library( magrittr )
#' library( ANTsRCore )
#' library( ANTsR )
#' temp1 = resampleImage( ri(1), 4 )
#' temp2 = resampleImage( ri(2), 4 )
#' temp3 = resampleImage( ri(1), 2 )
#' tarspace = makeImage( c( dim( temp1 ), 2 ), spacing = c( 4, 4, 1 ) )
#' pdcr = mergeListToNDImage( tarspace, list( temp1, temp2 ) )
#' reg = timeSeriesMotionCorrect( temp3, pdcr )
#'
#' @export timeSeriesMotionCorrect
timeSeriesMotionCorrect <- function( 
  fixedStaticImage, 
  movingTimeSeries, 
  totalSigma, 
  flowSigma,   
  regIterations = c( 50, 50, 10 ), 
  ... ) {
  if ( missing( totalSigma ) ) totalSigma = 6
  if ( missing( flowSigma ) ) flowSigma = 6
  motlist = splitNDImageToList( movingTimeSeries )
  for ( k in 1:length( motlist ) ) {

    reg = ANTsRCore::antsRegistration( fixedStaticImage, motlist[[k]], 'SyN', 
      regIterations = regIterations, 
      flowSigma = flowSigma, totalSigma = totalSigma, ... )

    motlist[[ k ]] = reg$warpedmovout

    }
  return( motlist )
}



#' Time series super resolution algorithms.
#'
#' Implements several different strategies for super-resolution from time series data.
#'
#' @param fixedStaticImage fixed target image
#' @param mappedTimeSeries pre-aligned to the fixed space
#' @param method one of average, sharpen, pysr, neurenhance, majon
#' @param weight generic weight term
#' @return super resolution image
#' @author Avants BB
#' @examples
#'
#' library( magrittr )
#' library( ANTsRCore )
#' library( ANTsR )
#' temp1 = resampleImage( ri(1), 8 )
#' temp2 = resampleImage( ri(1), 8 )
#' temp3 = resampleImage( ri(1), 4 )
#' tarspace = makeImage( c( dim( temp1 ), 2 ) , spacing = c( 8, 8, 1 ) )
#' pdcr = mergeListToNDImage( tarspace, list( temp1, temp2 ) )
#' reg = timeSeriesMotionCorrect( temp3, pdcr )
#' sup = timeSeriesSuperResolution( temp3, reg )
#'
#' @export timeSeriesSuperResolution
timeSeriesSuperResolution <- function( 
  fixedStaticImage, 
  mappedTimeSeries, 
  method = 'average', 
  weight = 0.5 ) {
  allowableMethods = c( 'average', 'sharpen', 'pysr', 'neurenhance', 'majon' )
  if ( ! ( method %in% allowableMethods ) )
    stop( paste( 'method', method, 'not allowed -- see help.' ) )
  if ( method == 'average' )
    return( antsAverageImages( mappedTimeSeries ) )
  if ( method == 'sharpen' ) {
    avg = antsAverageImages( mappedTimeSeries )
    return( avg * weight + iMath( avg, 'Sharpen' ) * ( 1.0 - weight )  )
    }
  if ( method == 'majon' ) {
    majicfn = 'x'
    majic = antsImageRead( majicfn )
    return( majic )
    }
  if ( method == 'pysr' | method == 'neurenhance' ) {
    avg = antsAverageImages( mappedTimeSeries )
    return( avg )
    }

}
