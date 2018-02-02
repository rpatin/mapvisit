#' Calculate number of visit and residence time for movement data
#'
#' Creates 4 maps with : number of visit in a given radius, total time spent, mean and median time per visit.
#' @param df movement data
#' @param radius for finding residence event
#' @param time_out maximum time (in s) individuals are allowed to leave the circle without considering a new visit.
#' @param basename rasters will be saved as basename_nvisit.tif,basename_totRT.tif,basename_meanRT.tif,basename_medianRT.tif
#' @param driver for rasters
#' @param coord.names coordinate names in df data.frame. default x,y
#' @param timecol name of time column in df data.frame
#' @param grid if you want to provide a raster for calculations.
#' @param cellsize if you dont provide any grid argument, size for the cell used for calculations. Else, will generate a grid with 50 cells for the highest dimension.
#' @param proj4string projection string. Default UTM35S.
#' @param progress (boolean) should the function display a progress bar
#' @useDynLib mapvisit
#' @importFrom Rcpp evalCpp
#' @return  a named list of raster object
#' @examples
#' map_residence(df, radius = 500, maxtime = 2, time_out = 1, basename = "file",
#' coord.names = c("x","y"), timecol = "dateTime",units ="hour", proj4string =
#' "+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84
#' +towgs84=0,0,0",cellsize = 250, progress = T)
#' @export
#'
map_residence <- function(df, radius = NULL,time_out = NULL, basename = NULL, driver = "GeoTIFF", coord.names = c("x","y"), timecol = "dateTime", grid = NULL, proj4string = "+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",cellsize = NULL, progress = T){
  x = df[,coord.names[1]]
  y = df[,coord.names[2]]
  spdf <- df
  sp::coordinates(spdf) <- coord.names
  sp::proj4string(spdf) <- proj4string
  sfdf <- sf::st_as_sf(spdf)
  sfbuffer <- sf::st_buffer(sfdf,dist=radius)
  time = df[,timecol]
  if(any(is.na(time))) stop("time should not contain NA")
  if(any(is.na(x))) stop("x should not contain NA")
  if(any(is.na(y))) stop("y should not contain NA")
  if(is.null(grid)){
    current_extent <- raster::extent(spdf)*1.25
    xmin = current_extent[1,]
    ymin = current_extent[3,]
    xmax = current_extent[2,]
    ymax = current_extent[4,]
    if(is.null(cellsize)) cellsize = max((xmax-xmin)/50,(ymax-ymin)/50)
    xcell = ceiling((xmax-xmin)/cellsize)
    ycell = ceiling((ymax-ymin)/cellsize)
    spGrid <- sp::SpatialGrid(sp::GridTopology(c(xmin,ymin),c(cellsize,cellsize),c(xcell,ycell)))
    # rast = raster(spGrid)
  } else {
    tmptopo <- sp::GridTopology(sp::bbox(grid)[,1] + res(grid)/2, res(grid), dim(grid)[2:1])
    spGrid <- sp::SpatialGrid(tmptopo)
  }

  # spdf_buffer <- rgeos::gBuffer(spdf,byid=)
  spCenter <- as.data.frame(sp::coordinates(spGrid))
  # sp::coordinates(spCenter) <- c("s1","s2")
  # sp::proj4string(spCenter) <- proj4string

  # sfCenter = sf::st_as_sf(spCenter)

  lgrid = length(spGrid)
  spGridDF = sp::SpatialGridDataFrame(spGrid,data = data.frame("nvisit" = rep(NA,lgrid), "totRT" = rep(NA,lgrid), "meanRT" = rep(NA,lgrid), "medianRT" = rep(NA,lgrid)),proj4string = CRS(proj4string))

  # if(progress) {
  #   pb <- progress::progress_bar$new(
  #   format = "  Calculations [:bar] :percent elapsed time: :elapsed ; estimated time left: :eta",
  #   total = length(spCenter), clear = FALSE, width= 120)
  # }

  # for(i in 1:length(spCenter)){
  #   if(progress) pb$tick()
  # message(i)
  # current_center = sfCenter[i,]
  # SEXP calc_visit(SEXP xyt, SEXP xygrid, SEXP distr, SEXP maxt)
  #
  resRT <-   calc_visit_cpp(x,y,df[,timecol], spCenter[,1], spCenter[,2], radius, time_out)

  # resRT <-   .Call("calc_visit_c",data.frame(x,y,time=df[,timecol]), spCenter, radius, time_out)


  # resRT <- calc_visit(current_center,sfdf,radius = radius, maxtime = maxtime, timecol = timecol, units = units,proj4string = proj4string,time_out = time_out)
# test
  spGridDF@data$nvisit <- resRT[["nvisit"]]
  spGridDF@data$totRT <- resRT[["residtime"]]
  spGridDF@data$meanRT <- resRT[["nvisit"]] / resRT[["residtime"]]
  spGridDF@data$medianRT <- NA
  # spGridDF@data$nvisit[i] <- i
  # spGridDF@data$totRT[i] <- i
  # spGridDF@data$meanRT[i] <- i
  # spGridDF@data$medianRT[i] <- i
  # }
  # spGridData <- lapply(1:nrow(sfCenter),function(c){
  #   current_center <- sfCenter[c,]
  #   resRT <- calc_visit(current_center,sfdf,radius = radius, maxtime = maxtime, timecol = timecol, units = units,proj4string = proj4string,time_out = time_out)
  #   return(data.frame("nvisit" = resRT[["nvisit"]],
  #                     "totRT" = resRT[["totRT"]],
  #                     "meanRT" = resRT[["meanRT"]],
  #                     "medianRT" = resRT[["medianRT"]]))
  # })
  # spGridDF@data <- do.call('rbind',spGridData)

  variables <- c("nvisit","totRT","meanRT")
  listRaster <- list()
  for(var in 1:length(variables)){
    rast = raster::raster(spGridDF,layer=var)
    sp::proj4string(rast) <- proj4string
    raster::writeRaster(rast,paste(basename,"_",variables[var],".tif",sep=""),driver=driver,overwrite=T)
    listRaster[[variables[var]]] <- rast
  }
  return(listRaster)
}

# mapvisit
