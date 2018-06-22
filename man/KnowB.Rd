\name{KnowB}
\alias{KnowB}
\title{Discriminating well surveyed cell units from exhaustive biodiversity databases}
\description{
Advances during the last decades in information technology allow us to store, retrieve, transmit and manipulate an unprecedented magnitude of massive information about species distributions (Guralnick \emph{et al}., 2007). Unfortunately, this compilation process suffers from three main shortcomings:

i) \emph{Unknown survey effort}. A lack of knowledge of the effort devoted to survey each territorial unit that is due to most occurrence records lacking any associated measure of the effort carried out to obtain them.

ii) \emph{Unknown absences}. As almost all the available information involves only species occurrences (i.e., the localities in which a species has been collected), without any indication of the likelihood that a species is actually absent from the localities where it was not collected (whether these have been surveyed or not).

iii) \emph{Unknown recurrence}. Which results from the incomplete compilation of species occurrences in many biodiversity databases, as multiple records of the same species in the same site or territorial unit are considered redundant and not reported (Hortal \emph{et al}. 2007); this prevents teasing apart occasional records from the continued presence of the species in an area.

These three limitations are mutually interrelated, so when all known occurrences are compiled exhaustively it is possible to estimate survey effort with some reliability. Therefore, a biodiversity database that compiles exhaustively all available information on the identity and distribution of a group of species would enable both identifying well-surveyed areas (e.g. Hortal and Lobo 2005) and obtaining estimates of the repeated occurrence and/or the probability of absence of particular species (e.g. Guillera-Arroita \emph{et al}. 2010).

Employing statistical shortcuts on data with unknown levels of error and bias can generate unreliable results. Consequently, good practice in biodiversity informatics requires knowledge about the number, location and degree of completeness of surveys for those territorial units that have been, at least relatively, well inventoried. Such knowledge would facilitate identifying localities where the lack of records for a target species can be reliably assumed to correspond to its actual absence. Nonetheless, it can be used to guide the location of future surveys and/or determine uncertain or ignorance areas in which biodiversity data are insufficiently consistent (Hortal and Lobo 2005, Ladle and Hortal 2013, Hortal \emph{et al}. 2015, Ruete 2015, Meyer \emph{et al}. 2015, 2016).

Despite the widely recognized importance of evaluating data quality as a preliminary step in any biodiversity study, this process is often neglected. Arguably, this is in part because such evaluation process is highly time-consuming, for it requires using analyses spread over several software applications and/or R packages, and repeating the same process for each one of the territorial units or sites considered (or, in general, for any type of spatial unit). Here we present KnowBR, a freely available R package to estimate the survey coverage of species inventories across an unlimited number of territorial units or sites simultaneously. Starting with any biodiversity database, KnowBR calculates the survey coverage per spatial unit as the final slope of the relationship between the number of collected species and the number of database records, which is used as a surrogate of the survey effort. To do this, KnowBR estimates the accumulation curve (the accumulated increase in the number of species with the addition of database records) for each one of the spatial units according to the \emph{exact} estimator of Ugland \emph{et al}. (2003), as well as performing 200 permutations of the observed data (\emph{random} estimator) to obtain a smoothed accumulation curve. This curve is subsequently adjusted to four different functions with three or less parameters, and the obtained extrapolated asymptotic value used to obtain a completeness percentage (the percentage representing the observed number of species against the predicted one) that also may be used to estimate the territorial units with probable reliable inventories.

These territorial units can be regular cells of any resolution (\emph{cell} option) but also irregular polygons (\emph{polygon} option) according to user preferences. RWizard includes in the "Area" argument the possibility of select the administrative spatial units (countries, regions, departments and/or provinces) or the rivers basins of different levels in which to perform the calculations. Instead of using the polygons available in RWizard, the user may also include any shapefile containing the desired irregular polygons (e.g. protected areas, countries, etc) by means of the "shape" argument.
}
\usage{
KnowB(data, format="A", cell=60, curve= "Rational", estimator=1, cutoff=1,
cutoffCompleteness= 0, cutoffSlope= 1, largematrix=FALSE, Area="World",
extent=TRUE, minLon, maxLon, minLat, maxLat, colbg="transparent",
colcon="transparent", colf="black", pro=TRUE, inc=0.005, exclude=NULL,
colexc=NULL, colfexc="black", colscale=c("#C8FFFFFF","#64FFFFFF","#00FFFFFF",
"#64FF64FF","#C8FF00FF","#FFFF00FF","#FFC800FF","#FF6400FF","#FF0000FF"),
legend.pos="y", breaks=9, xl=0, xr=0, yb=0, yt=0, asp, lab=NULL, xlab="Longitude",
ylab="Latitude", main1="Observed richness", main2="Records", main3="Completeness",
main4="Slope", cex.main=1.6, cex.lab=1.4, cex.axis=1.2, cex.legend=1.2,
family="sans", font.main=2, font.lab=1, font.axis=1, lwdP=0.6, lwdC=0.1,
trans=c(1,1), log=c(0,0), ndigits=0, save="CSV", file1="Observed richness",
file2="List of species", file3="Species per site", file4="Estimators",
file5="Species per record", file6="Records", file7="Completeness", file8="Slope",
file9="Standard error of the estimators", na="NA", dec=",", row.names=FALSE,
jpg=TRUE, jpg1="Observed richness.jpg", jpg2="Records.jpg", jpg3="Completeness.jpg",
jpg4="Slope.jpg",cex=1.5, pch=15, cex.labels=1.5, pchcol="red", ask=FALSE)
}
\arguments{
  \item{data}{The data is introduced as a CSV, TXT or RData file following two simple formats: one in which only four columns are included (see format A; species name, longitude, latitude and a number reflecting the incidence of the species) and another one including the longitude and latitude of each spatial unit and as many columns as species (see format B in the following table). The CSV file with the format A may be obtained using ModestR (see details).
\tabular{c}{
\figure{F1}{options: width=5.3in}
}
The primary matrix used in \emph{KnowBR} has a special characteristic - it must be derived from an exhaustive database including all the available georeferenced information including even those apparently redundant records of a species from the same locality provided that is a difference in some of the collection conditions for a species at a locality (i.e. date of capture, food source, collector, type of microhabitat, etc.). Thus, any difference in any database field value yields a new database record regardless of the number of individuals (see for example Lobo & Martín-Piera, 2002). As biodiversity data can derive from heterogeneous sources with different collector methodologies, no universal sampling effort measure capable of offering  reliable comparisons exists and the number of database records is used as a surrogate (see Soberón \emph{et al}., 2007; Lobo, 2008). This approach is particularly appropriate for poorly surveyed groups and/or regions lacking sufficient information to correct unequal sampling efforts arising from standardized survey protocols.}
  \item{format}{If it is "A" (default), the format of the data frame is species, longitude, latitude and a count value (format A of the table showed above). If it is "B" the format of the data frame is longitude, latitude and the rest of columns are the presence of the species in each site (format B of the table showed above). If numeric values higher than 1 are included in these data a (Count or Sp columns), a database record is considered for each unit. This in the example of format A, four different records are included for the Sp3 with same geographical coordinates.}
  \item{cell}{Resolution of the cells (spatial units) in minutes on which calculations were carried out. In the present version the user can select any resolution between 1 and 60 minutes.}
  \item{curve}{The smoothed accumulation curve generated by the accumulation curve can be adjusted to a "Clench", "Exponential", "Saturation" or "Rational" function (see equations in details section), calculating the asymptotic extrapolated values to further derive a completeness percentage (the percentage representing the observed number of species against the predicted one).}
  \item{estimator}{Vector that defines the used estimator:

0 The data for the estimation of the accumulation curve and the final slope are obtained with both the "exact" and "random" procedures. When the predicted richness is estimated with the type of curve selected by the user ("Clench", "Nexponential", "Saturation" or "Rational") using the data generated by the methods "exact" and "random" at the same time, the mean of both richness values is used to calculate completeness (the percentage representing the observed number of species against the predicted one).

1 It is the "exact" estimator of Ugland \emph{et al.} (2003) (default option) to obtain a smoothed accumulation curve.

2 If the chosen option is "random". It adds records at random performing 200 permutations in the order of records entry to generate the accumulation curve.

}
  \item{cutoff}{This number reflects the ratio between the number of database records and the number of species. If this ratio is lower than the selected threshold value in each considered spatial unit, any one of the estimators will be calculated and these spatial units are considered as lacking information.}
  \item{cutoffCompleteness}{If the value of completeness is lower than this threshold, the completeness is not calculated.}
  \item{cutoffSlope}{If the slope is higher than this threshold, the completeness is not calculated.}
  \item{largematrix}{When there many species and/or many records resulting in a species per record matrix with more than \eqn{2^{31}} cells, it is impossible to create the CSV or RData file with the species per record due to memory limits in R. If this argument is TRUE, the function creates a TXT file with the species per record, but the process is computationally intensive, so it may takes several hours and it may create a large TXT file. The default value is FALSE.}
  \item{Area}{A character with the name of the administrative area or a vector with several administrative areas. If a vector with several administrative areas are used, it is necessary to use RWizard (see details).}
  \item{extent}{If TRUE the minimum and maximum longitudes and latitudes are delimited by the minimum and maximum of the data (default). If FALSE the minimum and maximum longitudes and latitudes are delimited by the arguments Area and, minLat, maxLat, minLon and maxLon.}
  \item{minLon, maxLon}{Optionally it is possible to define the minimum and maximum longitude (see details).}
  \item{minLat, maxLat}{Optionally it is possible to define the minimum and maximum latitude (see details).}
  \item{colbg}{Background color of the map (in some cases this is the sea).}
  \item{colcon}{Background color of the administrative areas.}
  \item{colf}{Color of administrative areas border.}
  \item{pro}{If it is TRUE an automatic calculation is made in order to correct the
aspect ratio y/x along latitude.}
  \item{inc}{Adds some room along the map margins with the limits x and y thus not exactly the limits of the selected areas.}
  \item{exclude}{A character with the name of the administrative area or a vector with several administrative areas that can be plotted with a different color on the map.}
  \item{colexc}{Background color of areas selected in the argument exclude.}
  \item{colfexc}{Color of borders of the areas selected in the argument exclude.}
  \item{colscale}{Palette color.}
  \item{legend.pos}{Whether to have a horizontal (x) or vertical (y) color gradient.}
  \item{breaks}{Number of breakpoints of the color legend.}
  \item{xl,xr,yb,yt}{The lower left and upper right coordinates of the color legend in user coordinates.}
  \item{asp}{The y/x aspect ratio.}
  \item{lab}{A numerical vector of the form c(x. Y) which modified the default method by which axes are annotated. The values of x and y give the (approximate) number of tick marks on the x and y axes.}
  \item{xlab}{A title for the x axis.}
  \item{ylab}{A title for the y axis.}
  \item{main1}{An overall title for the plot of the observed species richness.}
  \item{main2}{An overall title for the plot of the records.}
  \item{main3}{An overall title for the plot of the completeness.}
  \item{main4}{An overall title for the plot of the slope between the last species richness value and the previous value for each one of the accumulation methods.}
  \item{cex.main}{The magnification to be used for main titles relative to the current setting of cex.}
  \item{cex.lab}{The magnification to be used for x and y labels relative to the current setting of cex.}
  \item{cex.axis}{The magnification to be used for axis annotation relative to the current setting of cex.}
  \item{cex.legend}{The magnification to be used in the numbers of the color legend relative to the current setting of cex.}
  \item{family}{The name of a font family for drawing text.}
  \item{font.main}{The font to be used for plot main titles.}
  \item{font.lab}{The font to be used for x and y labels.}
  \item{font.axis}{The font to be used for axis annotation.}
  \item{lwdP}{Line width of the plot.}
  \item{lwdC}{Line width of the borders.}
  \item{trans}{It is possible to multiply or divide the dataset by a value. For a vector with two values, the first may be 0 (divide) or 1 (multiply), and the second number is the value of the division or multiplication.}
  \item{log}{It is possible to apply a logarithmic transformation to the dataset. For a vector with two values, the first may be 0 (do not log transform) or 1 (log transformation), and the second number is the value to be added in case of log transformation.}
  \item{ndigits}{Number of decimals in legend of the color scale.}
  \item{save}{If "CSV" the files are save as CSV and if "RData" the files are save as RData.}
  \item{file1}{RData or CSV file. A character string naming the file with the observed richness.}
  \item{file2}{RData or CSV file. A character string naming the file with the list of species.}
  \item{file3}{RData or CSV file. A character string naming the file with the species incidences per site.}
  \item{file4}{RData or CSV file. A character string naming the file with the estimators per site.}
  \item{file5}{RData or CSV file. A character string naming the file with the species per records.}
  \item{file6}{RData or CSV file. A character string naming the file with the records.}
  \item{file7}{RData or CSV file. A character string naming the file with the completeness.}
  \item{file8}{RData or CSV file. A character string naming the file with the slopes of the accumulation analyses.}
  \item{file9}{RData or CSV file. A character string naming the file with the standard error of the estimators.}
  \item{na}{CSV FILE. Text that is used in the cells without data.}
  \item{dec}{CSV FILE. It defines if the comma "," is used as decimal separator or the dot ".".}
  \item{row.names}{CSV FILE. Logical value that defines if identifiers are put in rows or a vector with a text for each of the rows.}
  \item{jpg}{If TRUE the plots are exported to jpg files instead of using the windows device.}
  \item{jpg1}{Name of the jpg file with the values of the observed richness.}
  \item{jpg2}{Name of the jpg file with the records.}
  \item{jpg3}{Name of the jpg file with the completeness.}
  \item{jpg4}{Name of the jpg file with the slopes of the accumulation analyses.}
  \item{cex}{A numerical value giving the amount by which plotting symbols should be magnified relative to the default in the correlation matrix plot.}
  \item{pch}{Either an integer specifying a symbol or a single character to be used as the default in plotting points in the correlation matrix plot.}
  \item{cex.labels}{Size of labels in the correlation matrix plot.}
  \item{pchcol}{Color of the symbols in the correlation matrix plot.}
  \item{ask}{If TRUE (and the R session is interactive) the user is asked for input before a new figure is drawn.}
}
\details{
The CSV file required in the argument \emph{data} with the format A (species, longitude, latitude and count) may be obtained using ModestR (available at the web site www.ipez.es/ModestR) juast selecting Export/Export maps of the select branch/To RWizard Applications/To KnowBR.

In ModestR is possible to export the valid samples or pseudosamples. The pseudosamples are grid cells for instance of 5' x 5', 30' x 30', 1º x 1º, etc. Therefore, the output of ModestR is a list of species within each of the grid cells with the cell size defined by the user. It is therefore possible to obtain the number of records for each species within the grid cell or just the records available for all the species, with the format described above.

Area = "World" to plot the entire world. If the coordinates minLon, maxLon, minLat and maxLat are not specified, they are calculated automatically based on the selected administrative areas. If some administrative areas are selected, e.g. some countries, so the argument is not "World", it only works with RWizard.

It is important to emphasize that the quality of the geographical records of the administrative areas is lower if it is used the entire world(\emph{Area="World"}, because the file adworld is used), than if it is selected some countries, departments, etc., because the geographical records of the administrative areas available in RWizard are used. It means that the records inside the polygons may vary depending on the selection specified in the argument \emph{Area}.

The type of curves are:

1) The curve of Clench (Clench, 1979), which is a modification of the function of Monod (Monod, 1950), and was proposed to butterflies.

2) The exponential (Miller & Wiegert, 1989) that was proposed for rare plant species.

3) The saturation curve that was used to show the relationship between growth of phytoplankton, a toxic algae of the genus \emph{Alexandrium} and the concentration of phosphate (Frangópulos et al., 2004), which is similar to von Bertalanffy growth curve but adapting the coefficients to better explain the pattern of accumulation function.

4) The rational function (Ratkowski, 1990) that can be used when there is no clear criterion which model to use (Falther, 1996).
\tabular{lll}{
\strong{Name} \tab \strong{Function} \tab \strong{Reference}\cr
Clench \tab \eqn{y=\frac{ax}{1+bx}} \tab (Clench, 1979)\cr
Negative exponential \tab \eqn{y=a\left(1-e^{-bx}\right)} \tab (Miller & Wiegert, 1989)\cr
Saturation \tab \eqn{y=a\left(1-e^{-b(x-c)}\right)} \tab (Frangópulos et al., 2004)\cr
Rational \tab \eqn{y=\frac{(a+bx)}{(1+cx)}} \tab (Ratkowski, 1990)\cr
}
\strong{FUNCTIONS}

The estimators exact and random were estimated with the function \link[vegan]{specaccum} of the package vegan (Oksanen \emph{et al}., 2014).

The color legend of the maps is depicted with the function \link[plotrix]{color.legend} of the package plotrix (Lemon et al., 2014).

\strong{EXAMPLE}

The database of the example includes 15,142 records for the 54 Iberian species of the Scarabaeidae (Coleoptera) previously compiled in the so called BANDASCA database (Lobo & Martín-Piera, 2002). The following map show the slopes obtained in cells of 60ºx 60º using the estimator exact and the Rational's curve.
\tabular{c}{
\figure{F2}{options: width=5.0in}\cr
}
The maps may be easily modified using the function \link[KnowBR]{MapCell} using the exported CSV or RData files detailing the observed species richness (with alias ObservedRichness), the records (with alias Records), the completeness (with alias Completeness) and the slope (with alias Slope).
}
\value{
RData or CSV files: 1) Observed richness, 2) List of species, 3) Species per site, 4) Estimators, 5) Species per record, 6) Records, 7) Completeness, Slope and 9) Standard error of the estimators.

JPG files with maps: 1) Observed richness, 2) Records, 3) Completeness and 4) Slope.
}
\source{Spatial database of the location of the world's administrative areas (or administrative boundaries) was obtained from the Web Site http://www.openstreet.org/.}
\references{
Clench, H.K. (1979) How to make regional lists of butterflies: some invoking empirically based criteria in selecting among thoughts. \emph{The Journal of the Lepidopterists' Society}, 33: 216-231.

Flather, C.H. (1996) Fitting species-accumulation functions and assessing regional land use impacts on avian diversity. \emph{Journal of Biogeography}, 23: 155-168.

Frangópulos, M., Guisande, C., deBlas, E. y Maneiro, I. (2004) Toxin production and competitive abilities under phosphorus limitation of \emph{Alexandrium} species. \emph{Harmful Algae}, 3: 131-139.

Guillera-Arroita, G., Ridout, M.S & Morgan, B.J.T. (2010)Design of occupancy studies with imperfect detection. \emph{Methods in Ecology and Evolution}, 1: 131-139. 

Guralnick, R.P., Hill, A.W. & Lane, M. 2007. Towards a collaborative, global infrastructure for biodiversity assessment. \emph{Ecology Letters} 10: 663-672.

Hortal, J., de Bello, F., Diniz-Filho, J.A.F., Lewinsohn, T.M., Lobo, J.M. & Ladle, R.J. (2015) Seven shortfalls that beset large-scale knowledge of biodiversity. \emph{Annual Review Ecology and Systematics}, 46: 523-549.

Hortal, J. & Lobo, J.M. 2005. An ED-based protocol for the optimal sampling of biodiversity. \emph{Biodiversity and Conservation}, 14: 2913-2947.

Hortal, J., Lobo, J.M. & Jiménez-Valverde, A., 2007. Limitations of biodiversity databases: case study on seed-plant diversity in Tenerife (Canary Islands). \emph{Conservation Biology} 21, 853-863.

Ladle, R. & Hortal, J. (2013) Mapping species distributions: living with uncertainty. \emph{Frontiers of Biogeography}, 5: 8-9. 

Lemon, J., Bolker, B., Oom, S., Klein, E., Rowlingson, B., Wickham, H., Tyagi, A., Eterradossi, O., Grothendieck, G., Toews, M., Kane, J., Turner, R., Witthoft, C., Stander, J., Petzoldt, T., Duursma, R., Biancotto, E., Levy, O., Dutang, C., Solymos, P., Engelmann, R., Hecker, M., Steinbeck, F., Borchers, H., Singmann, H., Toal, T. & Ogle, D. (2017). Various plotting functions. R package version 3.6-5. Available at: http://CRAN.R-project.org/package=plotrix.

Lobo, J.M., Baselga, A., Hortal, J., Jiménez-Valverde, A. & Gómez, J.F. 2007. How does the knowledge about the spatial distribution of Iberian dung beetle species accumulate over time? \emph{Diversity and Distributions} 13:772-780.

Lobo, J.M. 2008. Database records as a surrogate for sampling effort provide higher species richness estimations. \emph{Biodiversity and Conservation} 17: 873-881.

Meyer, C., Kreft, H., Guralnick, R. & Jetz, W. (2015) Global priorities for an effective information basis of biodiversity distributions. \emph{Nature Communications}6: 8221.

Meyer, C., Weigelt, P. & Kreft, H. (2016) Multidimensional biases, gaps and uncertainties in global plant occurrence information. \emph{Ecoly Letters}, 19: 992-1006.

Miller, R.I. & Wiegert, R.G. (1989) Documenting completeness species-area relations, and the species-abundance distribution of a regional flora. \emph{Ecology}, 70: 16-22.

Oksanen, J., Blanchet, F.G., Kindt, R., Legendre, P., Minchin, P.R., O'Hara, R.B., Simpson, G.L., Solymos, P., Henry, M., Stevens, H. & Wagner, H. 2014. Community Ecology Package. R package version 2.0-10. Available at: https://CRAN.R-project.org/package=vegan.

Ratkowski, D.A. (1990) \emph{Handbook of nonlinear regression models}. Marcel Dekker, New York, 241 pp.

Ruete, A. (2015) Displaying bias in sampling effort of data accessed from biodiversity databases using ignorance maps.  \emph{Biodiversity Data Journal}3: e5361.

Soberón, J., Jiménez, R., Golubov, J. & Koleff, P., 2007. Assessing completeness of biodiversity databases at different spatial scales. \emph{Ecography} 30, 152-160.

Ugland, K.I., Gray J. S. & Ellingsen, K.E. 2003. The species-accumulation curve and estimation of species richness. \emph{Journal of Animal Ecology} 72: 888-897.

}
\examples{

#Example 1. Default conditions using estimator 1 (method exact)
#but only slopes lower than 0.1 are selected for depicting
#and, therefore, only the completeness is depicted for those
#cells with the slope lower than 0.1.
#If using RWizard, for a better quality of the geographic
#coordinates, replace data(adworld) by @_Build_AdWorld_

data(adworld)
data(Beetles)
KnowB(data=Beetles, save="RData", jpg=FALSE, cutoffSlope=0.1, xl=6.1, xr=6.3)

\dontrun{

#Only to be used with RWizard. 
#Example 2. Using @_Build_AdWorld_

data(Beetles)
@_Build_AdWorld_
KnowB(Beetles, cell=15, save="RData")

}
}
\keyword{KnowB}
