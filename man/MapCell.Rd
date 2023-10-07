\name{MapCell}
\alias{MapCell}
\encoding{latin1}
\title{Cell maps}
\description{It allows to depict on a map any of the variables  (records richness, observed richess, predicted richness, completeness and the slope) exported by the function \link[KnowBR]{KnowB} using a CSV, RData or raster file, and with the spatial resolution (cell size) specified in the file.}
\usage{
MapCell(data, Area="World", minLon, maxLon, minLat, maxLat, colbg="#FFFFFF",
colcon="#C8C8C8", colf="black", pro=TRUE, inc=0.005, exclude=NULL,
colexc=NULL, colfexc="black", colscale=c("#C8FFFFFF","#64FFFFFF","#00FFFFFF","#64FF64FF",
"#C8FF00FF","#FFFF00FF","#FFC800FF","#FF6400FF","#FF0000FF"),
legend.pos="y", breaks=9, xl=0, xr=0, yb=0, yt=0, asp, lab=NULL, xlab="Longitude",
ylab="Latitude", main=NULL,cex.main=1.2, cex.lab=1, cex.axis=0.9, cex.legend=0.9,
family="sans", font.main=2, font.lab=1, font.axis=1, lwdP=0.6, lwdC=0.1, trans=c(1,1),
log=c(0,0), ndigits=0, ini=NULL, end=NULL, jpg=FALSE, filejpg="Map.jpg")
}
\arguments{
  \item{data}{A CSV or RData file exported by the function \link[KnowBR]{KnowB} (see details) or an ESRI ASCII raster file with the variable (richness, records, etc.).}
  \item{Area}{Only if using RWizard. A character with the name of the administrative area or a vector with several administrative areas (countries, regions, etc.) or river basins. If it is "World" (default) the entire world is plotted. For using administrative areas or river basins, in addition to use RWizard, it is also necessary to replace data(adworld) by @_Build_AdWorld_ (see example 2).}
  \item{minLon, maxLon}{Optionally it is possible to define the minimum and maximum longitude.}
  \item{minLat, maxLat}{Optionally it is possible to define the minimum and maximum latitude.}
  \item{colbg}{Background color of the map (in some cases this is the sea).}
  \item{colcon}{Background color of the administrative areas.}
  \item{colf}{Color of administrative areas border.}
  \item{pro}{If it is TRUE an automatic calculation is made in order to correct the aspect ratio y/x along latitude.}
  \item{inc}{Adds some room along the map margins with the limits x and y thus not exactly the limits of the selected areas.}
  \item{exclude}{A character with the name of the administrative area or a vector with several administrative areas that may be plotted with a different color on the map (only if using RWizard).}
  \item{colexc}{Background color of areas selected in the argument exclude.}
  \item{colfexc}{Color of borders of the areas selected in the argument exclude.}
  \item{colscale}{Palette color.}
  \item{legend.pos}{Whether to have a horizontal "x" or vertical "y" color scale.}
  \item{breaks}{Number of breakpoints of the color legend.}
  \item{xl,xr,yb,yt}{The lower left and upper right coordinates of the color legend in user coordinates.}
  \item{asp}{The y/x aspect ratio.}
  \item{lab}{A numerical vector of the form c(x, y) which modifies the default way that axes are annotated. The values of x and y give the (approximate) number of tickmarks on the x and y axes.}
  \item{xlab}{A title for the X axis.}
  \item{ylab}{A title for the Y axis.}
  \item{main}{An overall title for the plot.}
  \item{cex.main}{The magnification to be used for main titles relative to the current setting of cex.}
  \item{cex.lab}{The magnification to be used for X and Y labels relative to the current setting of cex.}
  \item{cex.axis}{The magnification to be used for axis annotation relative to the current setting of cex.}
  \item{cex.legend}{The magnification to be used for the color scale relative to the current setting of cex.}
  \item{family}{The name of a font family for drawing text.}
  \item{font.main}{The font to be used for plot main titles.}
  \item{font.lab}{The font to be used for x and y labels.}
  \item{font.axis}{The font to be used for axis annotation.}
  \item{lwdP}{Line width of the plot.}
  \item{lwdC}{Line width of the borders.}
  \item{trans}{It is possible to multiply or divide the dataset by a value. For a vector with two values, the first may be 0 (divide) or 1 (multiply), and the second number is the value of the division or multiplication.}
  \item{log}{It is possible to apply a logarithmic transformation to the dataset. For a vector with two values, the first may be 0 (do not log transform) or 1 (log transformation), and the second number is the value to be added in case of log transformation.}
  \item{ndigits}{Number of decimals in legend of the color scale.}
  \item{ini}{Minimum to be considered in the color scale.}
  \item{end}{Maximum to be considered in the color scale.}
  \item{jpg}{If TRUE the plots are exported to jpg files instead of using the windows device.}
  \item{filejpg}{Name of the jpg file.}
}
\details{

It allows to depict on a map any of the files (CSV or RData) exported by the function \link[KnowBR]{KnowB}: records, observed richness, completeness and slope.

\strong{FUNCTIONS}

The function \link[plotrix]{color.legend} of the package plotrix (Lemon et al., 2014) is used for building the map.
}
\value{A map is obtained.}
\references{
Lemon, J. (2006) Plotrix: a package in the red light district of R. \emph{R-News}, 6(4):8-12.

Lemon, J., Bolker, B., Oom, S., Klein, E., Rowlingson, B., Wickham, H., Tyagi, A., Eterradossi, O., Grothendieck, G., Toews, M., Kane, J., Turner, R., Witthoft, C., Stander, J., Petzoldt, T., Duursma, R., Biancotto, E., Levy, O., Dutang, C., Solymos, P., Engelmann, R., Hecker, M., Steinbeck, F., Borchers, H., Singmann, H., Toal, T. & Ogle, D. (2015). Various plotting functions. R package version 3.6-1. Available at: https://CRAN.R-project.org/package=plotrix.
}
\examples{

\dontrun{

#Example 1. Observed pecies richness of freshwater fishes around the world.

data(RFishes)
data(adworld)
MapCell(data=RFishes, main= "Species richness of freshwater fishes")



#Example 2. Only to be used with RWizard.

data(RFishes)
@_Build_AdWorld_
MapCell(data = RFishes , Area = c("Argentina", "Bolivia", "Brazil", "Chile",
"Colombia", "Ecuador", "French Guiana", "Guyana", "Paraguay", "Peru", "Suriname",
"Uruguay", "Venezuela","Panama","Nicaragua","Costa Rica"),
main = "Species richness of freshwater fishes in South America")
}
}
\keyword{MapCell}
