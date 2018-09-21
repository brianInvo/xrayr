
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
#' @importFrom ANTsRCore ri iMath antsImageWrite antsGetSpacing
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
  if ( class( movingTimeSeries ) != 'list' )
    motlist = splitNDImageToList( movingTimeSeries ) else motlist = movingTimeSeries
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
#' @param fixedStaticImage target image space
#' @param mappedTimeSeries pre-aligned to the fixed space, not necessarily the same resolution
#' @param method one of average, sharpen, majon
#' @param param1 generic weight term for sharpen or mesh size for majon
#' @param param2 number of levels for majon
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
#' sup = timeSeriesSuperResolution( temp3, pdcr )
#'
#' @export timeSeriesSuperResolution
timeSeriesSuperResolution <- function( 
  fixedStaticImage,
  mappedTimeSeries, 
  method = 'average', 
  param1, param2 ) {
  allowableMethods = c( 'average', 'sharpen', 'majon' )
  if ( ! ( method %in% allowableMethods ) )
    stop( paste( 'method', method, 'not allowed -- see help.' ) )
  if ( class( mappedTimeSeries ) != 'list' )
    mappedTimeSeries = splitNDImageToList( mappedTimeSeries )
  if ( method == 'average' | method == 'sharpen' ) {
    if ( all( antsGetSpacing( mappedTimeSeries[[1]] ) == antsGetSpacing( fixedStaticImage )  )  ) 
      {
      avg = ( antsAverageImages( mappedTimeSeries ) ) 
      } else {
      avg = timeSeriesMotionCorrect( fixedStaticImage, mappedTimeSeries )
      avg = antsAverageImages( avg )
      }
      if ( method == 'sharpen' ) {
        return( avg * param1 + iMath( avg, 'Sharpen' ) * ( 1.0 - param1 )  )
      } else return( avg )
    }
  if ( method == 'majon' ) {
    majicfn = '/tmp/majic.nii.gz'
    reffn='/tmp/ref.nii.gz'
    antsImageWrite( fixedStaticImage, reffn )
    fns=paste0( '/tmp/reg',1:length(mappedTimeSeries),'.nii.gz')
    for ( fn in 1:length( fns ) ) antsImageWrite( mappedTimeSeries[[fn]], fns[fn] )
    fns=paste( paste0( '/tmp/reg',1:length(mappedTimeSeries),'.nii.gz'), collapse=' ' )
    cmd=paste('SuperResolution ', fixedStaticImage@dimension,' ',majicfn,' ',reffn,' -1 ',param1, param2, ' ',fns )
    print( cmd )
    system( cmd )
    majic = antsImageRead( majicfn )
    return( majic )
    }

}



#' Static super resolution algorithms.
#'
#' Implements several different strategies for super-resolution from static data.
#'
#' @param inputImage image to be expanded
#' @param expansionFactor float value
#' @param method one of sharpen, pysr, neurenhance
#' @return super resolution image
#' @author Avants BB
#' @examples
#'
#' library( magrittr )
#' library( ANTsRCore )
#' library( ANTsR )
#' temp1 = resampleImage( ri(1), 8 )
#' sup = staticSuperResolution( temp1, 2 )
#'
#' @export staticSuperResolution
staticSuperResolution <- function( 
  inputImage, 
  expansionFactor, 
  method = 'sharpen' ) {
  allowableMethods = c( 'sharpen', 'pysr', 'neurenhance' )
  if ( ! ( method %in% allowableMethods ) )
    stop( paste( 'method', method, 'not allowed -- see help.' ) )
  if ( method == 'sharpen' ) return( 
    iMath( resampleImage( inputImage, 
      antsGetSpacing( inputImage )/expansionFactor, useVoxels=FALSE ), "Sharpen" ) )
  if ( method == 'neurenhance' | method == 'pysr' ) {
    # must first call 
#    system( 'source ~/code/tensorflow/venv/bin/activate' )
    ANTsRCore::antsImageWrite( ANTsRCore::antsImageClone( iMath( inputImage, "Normalize") * 255, 'unsigned char' ), '/tmp/temp.jpg')
    cmd='python3 /home/bavants/code/Image-Super-Resolution/main.py /tmp/temp.jpg --scale=2'
    system( cmd )
    return( ANTsRCore::antsImageRead( '/tmp/temp_scaled(2x).jpg' ) )
#    cmd=paste0('enhance --zoom=',expansionFactor,' --model=repair /tmp/temp.jpg ')
    }
  return( NA )
}
