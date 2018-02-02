#'  Generate an empty grid based on total extent of data
#'
#' Returns a grid containing all data provided.
#' @param spdf movement data
#' @param cellsize resoluto
#' @param proj4string projection string. Default UTM35S.
#' @return  a grid
#' @examples
#' generate_grid(df, cellsize = 1000)
#' @export
#'
generate_grid <- function(spdf, cellsize = 1000, expand = 1.05){
  sp::coordinates(spdf) <- c("x","y")
  sp::proj4string(spdf) <- UTMstring
  tot_extent <- raster::extent(spdf)*1.05

  xmin = tot_extent[1,]
  ymin = tot_extent[3,]
  xmax = tot_extent[2,]
  ymax = tot_extent[4,]
  xcell = ceiling((xmax-xmin)/cellsize)
  ycell = ceiling((ymax-ymin)/cellsize)

  spGrid <- sp::SpatialGrid(sp::GridTopology(c(xmin,ymin),c(cellsize,cellsize),c(xcell,ycell)))
  totgrid <- raster::raster(spGrid)
  totgrid
}

#' Runs map_residence for a list individuals
#'
#' Returns a grid containing all data provided.
#' @param spdf movement data
#' @param cellsize resoluto
#' @param proj4string projection string. Default UTM35S.
#' @return  a grid
#' @examples
#' generate_grid(df, cellsize = 1000)
#' @importFrom magrittr "%>%"
#' @export
#'
run_mapvisit <- function(df, listIndiv, folder, grid, IDcol = "id", timecol = "dateTime", coord.names = c("x","y"), ...){
  df <- df[!duplicated(dplyr::select_at(df,c(IDcol,timecol))),]

  for(individu in listIndiv){
    subdata <- dplyr::filter_at(df, vars(IDcol), any_vars(. == individu))
    subdata <- dplyr::arrange_at(subdata,timecol)
    pb <- difftime(subdata[2:(nrow(subdata))  , timecol],
                   subdata[1:(nrow(subdata)-1), timecol],units='min') < 10
    if(any(pb))  subdata <- subdata[-which(pb),]
    # subdata <- subdata[!]
    # pb <- which(diff(subdata$dateTime) < 30)
    # head(subdata[rep(pb,each=3)+c(-1,0,1),])
    message(individu)
    tryCatch(listRaster <- map_residence(subdata, basename = paste0(folder,"/",individu), grid = grid, timecol = timecol,  ...))
  }
  return(invisible(NULL))
}

#' Plot number of visit vs total time spent for a dataset
#'
#' Returns and save a plot
#' @param df movement data
#' @param listIndiv list of individuals
#' @param folder folder where mapvisit data are stored
#' @param filename name of pdf file
#' @param width width of pdf
#' @param height height of pdf
#' @param IDcol names of id columne
#' @param coord.names names of coordinates columns
#' @return  a ggplot2 object
#' @examples
#' generate_grid(df, cellsize = 1000)
#' @importFrom magrittr "%>%"
#' @export
#'
plot_mapvisit <- function(df,listIndiv, folder, filename,width = 8, height = 8, IDcol = "id", coord.names = c("x","y"), timecol = "dateTime"){
  df <- df[!duplicated(dplyr::select_at(df,c(IDcol,timecol))),]

  pdf(file = filename,width = width, height = height)
  ID = first(df$id)
  for(individu in unique(df[,IDcol])){
    message(individu)
    subdata <- dplyr::filter_at(df, vars(IDcol), any_vars(. == individu))
    subdata <- dplyr::arrange_at(subdata,timecol)
    nvisit_std = raster::raster(paste0(folder,individu,"_nvisit_std.tif"))
    totRT_std  = raster::raster(paste0(folder,individu,"_totRT_std.tif"))
    meanRT     = raster::raster(paste0(folder,individu,"_meanRT.tif"))
    dvisit = nvisit_std[]
    visited = which(dvisit != 0)
    dfplot = data.frame(nvisit = dvisit[visited],totRT = totRT_std[visited])
    g <- ggplot2::ggplot(dfplot)+
      ggplot2::geom_point(ggplot2::aes(y=nvisit, x = totRT))+
      ggplot2::ggtitle(individu)+
      ggplot2::ylab("Number of visit per year")+
      ggplot2::xlab("Total Resident time (day/year)")+
      ggplot2::ylim(0,max(dfplot$nvisit))
    gridExtra::grid.arrange(g)
  }
  dev.off()
  return(g)
}
#' Merge nvisit and totRT for several individuals
#'
#' Returns and save a plot
#' @param df movement data
#' @param listIndiv list of individuals
#' @param folder folder where mapvisit data are stored
#' @param filename name of pdf file
#' @param width width of pdf
#' @param height height of pdf
#' @param IDcol names of id columne
#' @param coord.names names of coordinates columns
#' @param basename name for merged object
#' @return  a ggplot2 object
#' @examples
#' generate_grid(df, cellsize = 1000)
#' @importFrom magrittr "%>%"
#' @export
#'
calc_meanmapvisit <- function(df,listIndiv, basename, folder, filename,width = 8, height = 8, IDcol = "id", coord.names = c("x","y"), timecol = "dateTime", burstcol = IDcol){


  for(individu in listIndiv){
    message(individu)
    subdata <- dplyr::filter_at(df, vars(IDcol), any_vars(. == individu))

    subdata %>%
      dplyr::group_by_at(burstcol) %>%
      dplyr::summarise_at(vars(timecol),function(t){last(t)-first(t)}) %>%
      dplyr::rename_at(timecol,function(x){"duration"}) -> df.duration

    duration_tot = sum(as.numeric(df.duration$duration))
    time_factor <- 365/duration_tot

    nvisit = raster::raster(paste0(folder,"/",individu,"_nvisit.tif"))
    nvisit_std = calc(nvisit,fun = function(x) x*time_factor,filename = paste0(folder,"/",individu,"_nvisit_std.tif"),overwrite = T)

    totRT = raster::raster(paste0(folder,"/",individu,"_totRT.tif"))
    totRT_std = calc(totRT,fun = function(x) x*time_factor/60/60/24,filename= paste0(folder,"/",individu,"_totRT_std.tif"),overwrite = T)

    meanRT = raster::raster(paste0(folder,"/",individu,"_meanRT.tif"))

    ntot = length(listIndiv)
    if(individu == first(listIndiv)){
      mean_nvisit = calc(nvisit_std, fun = function(x) x/ntot, filename = paste0(basename,"_nvisit.tif"),overwrite = T)
      mean_totRT  = calc(totRT_std,  fun = function(x) x/ntot, filename = paste0(basename,"_totRT.tif"),overwrite = T)
      mean_meanRT = calc(meanRT,     fun = function(x) x/ntot, filename = paste0(basename,"_meanRT.tif"),overwrite = T)
    } else {
      mean_nvisit = raster::overlay(mean_nvisit,nvisit_std, fun = function(x,y) x+y/ntot, filename = paste0(basename,"_nvisit.tif"),overwrite = T)
      mean_totRT  = raster::overlay(mean_totRT,totRT_std,   fun = function(x,y) x+y/ntot, filename = paste0(basename,"_totRT.tif"),overwrite = T)
      mean_meanRT = raster::overlay(mean_meanRT,meanRT,     fun = function(x,y) x+y/ntot, filename = paste0(basename,"_meanRT.tif"),overwrite = T)
    }
  }


  dvisit = mean_nvisit[]
  visited = which(dvisit != 0)
  dfplot = data.frame(nvisit = dvisit[visited],totRT = mean_totRT[visited])
  g <- ggplot2::ggplot(dfplot)+
    ggplot2::geom_point(aes(y=nvisit, x = totRT))+
    ggplot2::ggtitle(individu)+
    ggplot2::ylab("Number of visit per year")+
    ggplot2::xlab("Total Resident time (day/year)")+
    ggplot2::ylim(0,max(dfplot$nvisit))
  cowplot::save_plot(filename,g,base_width = 8,base_height = 8)
  return(g)
}
