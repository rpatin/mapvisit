



#' Calculation of visits around a point.
#'
#' Takes a points and calculates the number, total, mean and median time of visits inside a given distance of this point.
#' @param current_center spatialPoints on which visit calculations should be done
#' @param spdf spatialPointsDataFrame with movement data.
#' @param timecol name of spdf@data time column
#' @param radius radius for visit
#' @param maxtime maximum time to consider relocation as successive and add time before first point and after last point proportional to the part of the segment inside the circle.
#' @param time_out maximum time individuals are allowed to leave the circle without considering a new visit
#' @param units unit used for maxtime and time_out
#' @param proj4string projection string. Default UTM35S.
#' @return  a named list with data.
#'
#' @examples
#' calc_visit(current_center, spdf ,radius = 500, maxtime = 5, time_out = 2, timecol = "dateTime", units = "hour", proj4string = "+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs")
#' @importFrom magrittr "%>%"
#' @useDynLib mapvisit
#' @export
# current_center = sfCenter[53,]
# current_center = sfCenter[1,]

calc_visit <- function(current_center, sfdf ,radius = NULL, maxtime = NULL, time_out = NULL, timecol = NULL, units = NULL, proj4string = NULL){
  # current_center = spCenter[11,]
  # sp_circle <- as(circle,"Spatial")
  circle = sf::st_buffer(current_center,dist=radius)
  # sp::proj4string(circle) <- proj4string
  pointsIntersection <- sf::st_intersects(circle,sfdf)[[1]]
  # sppointsIntersection <- rgeos::gIntersects(sp_circle,spdf,byid=T)
  if(length(pointsIntersection) == 0){
    return(list("nvisit"=0, "totRT" = 0, "meanRT" = 0, "medianRT" = 0))
  } else {
    tmp <- as.data.frame(sfdf)
    ldata= nrow(sfdf)
    pointsTF <- 1:ldata %in% pointsIntersection
    rupt = which(pointsTF[1:(ldata-1)] != pointsTF[2:(ldata)])
    sfdf$inside <- pointsTF

    seqList <- rle(pointsTF)
    start <- c(0,cumsum(seqList$lengths)[1:length(rupt)])+1
    names(start) <- NULL
    end <- cumsum(seqList$lengths)
    names(end) <- NULL
    value <- seqList$values
    inside <- which(value)
    names(inside)  <- NULL
    # merge visits less than time_out apart
    start_inside <- start[dplyr::first(inside)]
    end_inside <- c()
    if(length(inside) > 1){
      for(i in 2:length(inside)){
        index_leave <- end[inside[i-1]]
        time_leave <- tmp[index_leave,timecol]
        index_reenter <- start[inside[i]]
        time_reenter <- tmp[index_reenter,timecol]
        difft <- as.numeric(difftime(time_reenter,time_leave,units=units))
        # message(i)
        if(difft > time_out){
            end_inside <- c(end_inside,index_leave)
            start_inside <- c(start_inside,index_reenter)
        }
      }
    }
    end_inside <- c(end_inside,end[dplyr::last(inside)])
    nvisit <-  length(start_inside)
    times <- as.numeric(difftime(tmp[end_inside,timecol],tmp[start_inside,timecol],units = units))

    for(i in 1:nvisit){
      index_enter <- start_inside[i]
      index_leave <- end_inside[i]
      if(index_enter > 1){
        difft_enter = as.numeric(difftime(tmp[index_enter,timecol],tmp[index_enter-1,timecol],units = units))
        if( difft_enter <= maxtime ){
          line <- sfdf[c(index_enter,index_enter-1),] %>% sf::st_coordinates() %>% sf::st_linestring() %>% sf::st_sfc() %>% sf::st_sf(data.frame("id"=1,"geometry"=.)) %>% sf::st_set_crs(proj4string)

          # plot(circle,add=T)
          line_int <- sf::st_intersection(line,circle)
          percent_time <- sf::st_length(line_int)/sf::st_length(line)
          times[i]<- times[i] + difft_enter * percent_time[[1]]
        }
      }
      if(index_leave < ldata){
        difft_leave = as.numeric(difftime(tmp[index_leave+1,timecol],tmp[index_leave,timecol],units = units))
        if( difft_leave <= maxtime ){
          line <- sfdf[c(index_leave+1,index_leave),] %>% sf::st_coordinates() %>% sf::st_linestring() %>% sf::st_sfc() %>% sf::st_sf(data.frame("id"=1,"geometry"=.)) %>% sf::st_set_crs(proj4string)

          # plot(circle,add=T)
          line_int <- sf::st_intersection(line,circle)
          percent_time <- sf::st_length(line_int)/sf::st_length(line)
          times[i]<- times[i] + difft_leave * percent_time[[1]]
        }
      }
    }
    totRT <- sum(times)
    meanRT = totRT/nvisit
    medianRT = median(times)
    # trouver des régions consécutives de "inside". Les numéroter. Regarder celles qui sont distantes de moins de time_out. Fusionner les régions ainsi distante. Compter le temps pour chaque visite en rediscrétisant les segments d'entrée et de fin. Renvoyer le nombre de visite et le temps de residence total, median et moyen.
    return(list("nvisit"=nvisit, "totRT" = totRT, "meanRT" = meanRT, "medianRT" = medianRT))
  }
}

# mapvisit
