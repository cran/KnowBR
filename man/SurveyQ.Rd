\name{SurveyQ}
\alias{SurveyQ}
\encoding{latin1}
\title{Survey quality}
\description{Discriminations among good, fair and poor quality of surveys in cells and polygons.}
\usage{
SurveyQ(data, Longitude=NULL, Latitude=NULL, cell=60, Areas=NULL,
variables=c("Slope","Completeness","Ratio"), completeness=c(50,90),
slope=c(0.02,0.3), ratio=c(3,15), shape=NULL, shapenames=NULL, admAreas=TRUE,
Area="World", minLon, maxLon, minLat, maxLat, main=NULL, PLOTP=NULL,
PLOTB=NULL, POINTS=NULL, XLAB=NULL, YLAB=NULL, XLIM=NULL, YLIM=NULL,
palette=c("blue","green","red"), COLOR=c("red","green","blue"), colm="black",
labels=TRUE, sizelabels=1, LEGENDP=NULL, LEGENDM=NULL, file="Polar coordinates.csv",
na="NA", dec=",", row.names=FALSE, jpg=FALSE, filejpg="Map.jpg")
}
\arguments{
  \item{data}{Data file exported by the functions \link[KnowBR]{KnowBPolygon} or \link[KnowBR]{KnowB} named "Estimators" with the values of records, observed richness, predicted richness, completeness and slope for each area polygon.}
  \item{Longitude}{Variable with the longitude of the cells if the file "Estimators" was obtained with the function \link[KnowBR]{KnowB}, so it is the case when working with cells.}
  \item{Latitude}{Variable with the latitude of the cells if the file "Estimators" was obtained with the function \link[KnowBR]{KnowB}, so it is the case when working with cells.}
  \item{cell}{Resolution of the cells (spatial units) in minutes on which calculations were carried out. In the present version the user can select any resolution between 1 and 60 minutes.}
  \item{Areas}{Variable with the names of the polygons, if the file "Estimators" was obtained with the function \link[KnowBR]{KnowBPolygon}, so it is the case when working with polygons.}
  \item{variables}{The slope, completeness and ratio obtained in the file "Estimators", in that order.}
  \item{completeness}{Values of the completeness to define the thresholds for poor, fair and good quality surveys of the cells or polygons.}
  \item{slope}{Values of the slope to define the thresholds for poor, fair and good quality surveys of the cells or polygons.}
  \item{ratio}{Values of the ratio to define the thresholds for poor, fair and good quality surveys of the cells or polygons.}
  \item{shape}{If the estimators in the function link[KnowBR]{KnowBPolygon} were calculated using an external shape file, it is necessary to indicate the file in this argument. It is not necessary to select any polygon within the file, just to load the whole shape file.}
  \item{shapenames}{Variable in the shapefile with the names of the polygons.}
  \item{admAreas}{If it is TRUE the border lines of the countries are depicted in the map.}
  \item{Area}{Only if using RWizard (\url{http://www.ipez.es/RWizard}). A character with the name of the administrative area or a vector with several administrative areas (countries, regions, etc.) or river basins. If it is "World" (default) the entire world is plotted. For using administrative areas or river basins, in addition to use RWizard, it is also necessary to replace data(world) by @_Build_AdWorld_ (see examples of function \link[KnowBR]{KnowBPolygon}).}
\item{minLon, maxLon}{Optionally it is possible to define the minimum and maximum longitude.}
  \item{minLat, maxLat}{Optionally it is possible to define the minimum and maximum latitude.}
  \item{main}{Main title of the map.}
  \item{PLOTP}{It allows to specify the characteristics of the function \link[graphics]{plot.default} of the polar coordinates plot.}
  \item{PLOTB}{It allows to specify the characteristics of the function \link[graphics]{plot.default} of the bubble chart.}
  \item{POINTS}{It allows to modify the points of the bubble chart with the function \link[graphics]{points}.}
  \item{XLAB}{Legend of the X axis.}
  \item{YLAB}{Legend of the Y axis.}
  \item{XLIM}{Vector with the limits of the X axis.}
  \item{YLIM}{Vector with the limits of the Y axis.}
  \item{palette}{The color gradient of the bubble chart may be one of these palettes: "heat.colors", "terrain.colors", "gray.colors", "topo.colors" or "cm.colors".}
  \item{COLOR}{It allows to modify the colors of the map and polar coordinates plot. It must be three colors.}
  \item{colm}{Color of the polygons without information when using when using an external shape file.}
  \item{labels}{If it is FALSE, points are depicted instead of the labels of the polygons in the polar coordinates plot.}
  \item{sizelabels}{Text size of the labels of the polygons in the polar coordinates plot.}
  \item{LEGENDP}{It allows to modify the legend of the polar coordinates plot.}
  \item{LEGENDM}{It allows to modify the legend of the map.}
  \item{file}{CSV FILES. Filename with the polar coordinates.}
  \item{na}{CSV FILE. Text that is used in the cells without data.}
  \item{dec}{CSV FILE. It defines if the comma "," is used as decimal separator or the dot ".".}
  \item{row.names}{CSV FILE. Logical value that defines if identifiers are put in rows or a vector with a text for each of the rows.}
  \item{jpg}{If TRUE the map is exported to jpg files instead of using the windows device.}
  \item{filejpg}{Name of the jpg file.}
}
\value{It is depicted a plot with the polar coordinates of each polygon or cell, a map with the quality survey of the cells or polygons and a file with the polar coordinates of the cells or polygons.}
\details{
This function has been designed to identify and plot the cells or polygons with good, fair and poor quality surveys. This function uses the file called "Estimators" obtained from the functions \link[KnowBR]{KnowBPolygon} or \link[KnowBR]{KnowB} to estimate the polar coordinates of all cells or polygons and to discriminate among cells or polygons according to the quality of the survey.

The variables used by this function are slope, completeness and ratio (number of records/species observed). The default values to identify the cells or polygons with good, fair and poor quality surveys are: slope lower than 0.02, completeness higher than 90\% and ratio higher than 15 for good quality surveys, and slope higher than 0.3, completeness lower than 50\% and ratio lower than 3 for poor quality surveys.

The order of the variables is important for the estimation of the polar coordinates because a different angle is assigned to each variable. Therefore, the variables must be introduced in this order: slope, completeness and ratio.

All variables are transformed to a scale ranged between -1 and 1. For each value the X and Y polar coordinates are estimated using the following equations:
\tabular{ll}{
\eqn{X=\displaystyle\sum_{i=1}^{3}|z_j|cos(\alpha)} \tab \eqn{Y=\displaystyle\sum_{i=1}^{3}|z_j|sin(\alpha)}
}
where \emph{z} is the value of the variable \emph{j}.

Each variable is assigned an angle (\eqn{\alpha}). The increment value of the angle is always 60. Therefore, the first variable (slope) if the transformed value is \eqn{\geq0} the \eqn{\alpha} value is 60 and if the transformed value is \eqn{<0} the value is 240.

For the second variable (completeness) if the transformed value is \eqn{\geq0} the \eqn{\alpha} value is 120 and if the value is \eqn{<0} the value is 300.

For the third variable (ratio) if the value is \eqn{\geq0} the \eqn{\alpha} value is 180 and if the transformed value is \eqn{<0} the value is 360.

Degrees to radians angle conversion is carried out assuming that 1 degree = pi/180 radians.

\strong{EXAMPLES}
Polar coordinates of the records of freshwater fish species in all countries of the world.
\tabular{c}{
\if{html}{\figure{F6.jpg}{options: width=5.0in}}
\if{latex}{\figure{F6.jpg}{options: width=5.0in}}
}
Quality survey of the records of freshwater fish species in all countries of the world.
\tabular{c}{
\if{html}{\figure{F7.jpg}{options: width=5.5in}}
\if{latex}{\figure{F7.jpg}{options: width=5.5in}}
}
Bubble chart of the relationship between log(Ratio) and completeness, being the color gradient the slope value.
\tabular{c}{
\if{html}{\figure{F8.jpg}{options: width=4in}}
\if{latex}{\figure{F8.jpg}{options: width=4in}}
}
}
\references{
Guisande, C., Manjarrés-Hernández, A., Pelayo-Villamil, P., Granado-Lorencio, C., Riveiro, I., Acuña, A., Prieto-Piraquive, E., Janeiro, E., Matías, J.M., Patti, C., Patti, B., Mazzola, S., Jiménez, S., Duque, V. & Salmerón, F. (2010) IPez: An expert system for the taxonomic identification of fishes based on machine learning techniques. \emph{Fisheries Research}, 102, 240-247.
}
\examples{
\dontrun{

data(adworld)
data(Estimators)
SurveyQ(data=Estimators, Areas="Area")
}

}
\keyword{SurveyQ}
