% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{plot_mapvisit}
\alias{plot_mapvisit}
\title{Plot number of visit vs total time spent for a dataset}
\usage{
plot_mapvisit(df, listIndiv, folder, filename, width = 8, height = 8,
  IDcol = "id", coord.names = c("x", "y"))
}
\arguments{
\item{df}{movement data}

\item{listIndiv}{list of individuals}

\item{folder}{folder where mapvisit data are stored}

\item{filename}{name of pdf file}

\item{width}{width of pdf}

\item{height}{height of pdf}

\item{IDcol}{names of id columne}

\item{coord.names}{names of coordinates columns}
}
\value{
a ggplot2 object
}
\description{
Returns and save a plot
}
\examples{
generate_grid(df, cellsize = 1000)
}
