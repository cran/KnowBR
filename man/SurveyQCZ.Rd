\name{SurveyQCZ}
\alias{SurveyQCZ}
\encoding{latin1}
\title{Survey quality of climate zones}
\description{Estimation of the survey quality of different climate zones.}
\usage{
SurveyQCZ(data, Longitude="Longitude", Latitude="Latitude", cell=NULL, hull=TRUE,
Area="World", shape=NULL, shapenames=NULL, aprox=TRUE, VIDTAXA=NULL, por=80, k=NULL,
VIF=FALSE, VARSEDIG=FALSE, BUBBLE=FALSE, variables=c("Slope","Completeness","Ratio"),
completeness=c(50,90), slope=c(0.02,0.3), ratio=c(3,15), minLon=NULL,
maxLon=NULL, minLat=NULL, maxLat=NULL, xlab="Longitude", ylab="Latitude",
colscale=c("#C8FFFFFF","#64FFFFFF","#00FFFFFF","#64FF64FF","#C8FF00FF",
"#FFFF00FF","#FFC800FF","#FF6400FF","#FF0000FF"), colcon="transparent",
breaks=10, ndigits=0, xl=0, xr=0, mfrowBOXPLOT=NULL, mfrowMAP=NULL,
main="Percentage of ignorance/poor\n cells in each cluster", cexCM=0.5,
legpos="bottomleft", jpg=FALSE, filejpg="Survey Quality CZ.jpg", dec=",")
}
\arguments{
  \item{data}{Data file exported by the function \link[KnowBR]{KnowB} named "Estimators" with the values of records, observed richness, predicted richness, completeness and slope for each area polygon. This file may include the values of the environmental variables for each cell.}
  \item{Longitude}{Variable with the longitude.}
  \item{Latitude}{Variable with the latitude.}
  \item{cell}{Resolution of the cells (spatial units) in minutes. In the present version the user can select any resolution between 1 and 60 minutes. If it is NULL (default), it is used the resolution of the ASC files with the environmental variables. It is not NULL, the ASC files are rescaled to the resolution selected by the user.}
  \item{hull}{If it is TRUE, the extent is estimated as the convex null of the data. It is FALSE, it is used as extent the polygons defined in the arguments \emph{shape} of \emph{Area}.}
  \item{Area}{Only if using RWizard. It allows to plot countries, regions, river basins, etc., as extentd. A character with the name of the administrative area or a vector with several administrative areas (countries, regions, etc.) or river basins. If it is "World" (default) the entire world is plotted. For using administrative areas or river basins, in addition to use RWizard, it is also necessary to replace data(world) by @_Build_AdWorld_.}
    \item{shape}{Optionally it may be used a shape file with the information of the polygons.}
  \item{shapenames}{Variable in the shapefile with the names of the polygons.}
  \item{aprox}{If it is TRUE and there is not environmental data available in the cell, it is used the nearest cell available in the environmental data set.}
  \item{VIDTAXA}{It accesses the VIDTAXA function of the VARSEDIG package.}
  \item{por}{Cut-off threshold specifying the cumulative variance percentage, to determine how many axes are selected from the Principal Components or Correspondence analyses. By default it is 80\%, which means that the axes are selected until reaching an accumulated variance percentage of 80\%.}
  \item{k}{Number of clusters in which the Dendrogram is divided. If it is NULL, the algorithm select automatically the maximum number of clusters in which the Dendrogram can be divided, which are those groups that are statistically different in at least one variable according to the U Mann-Whitney test. If the are is large, there may be many different climate zones, so with NULL option running time may be long.}
  \item{VIF}{If it is TRUE, the inflation factor of the variance (VIF) is used to select the highly correlated variables and, therefore, not correlated variables are excluded from the Principal Components analysis.}
  \item{VARSEDIG}{If it is TRUE, the VARSEDIG algorithm is performed.}
  \item{variables}{The slope, completeness and ratio obtained in the file "Estimators", in that order.}
  \item{BUBBLE}{If it is TRUE, the BUBBLE plot the VARSEDIG function is depicted.}
  \item{completeness}{Values of the completeness to define the thresholds for poor, fair and good quality surveys of the cells.}
  \item{slope}{Values of the slope to define the thresholds for poor, fair and good quality surveys of the cells.}
  \item{ratio}{Values of the ratio to define the thresholds for poor, fair and good quality surveys of the cells.}
  \item{minLon, maxLon}{Optionally it is possible to define the minimum and maximum longitude.}
  \item{minLat, maxLat}{Optionally it is possible to define the minimum and maximum latitude.}
  \item{xlab}{Legend of the X axis in the map.}
  \item{ylab}{Legend of the Y axis in the map.}
  \item{colscale}{Color of the bar scale.}
  \item{colcon}{Background color of the administrative areas.}
  \item{breaks}{Number of breakpoints of the color legend.}
  \item{ndigits}{Number of decimals in legend of the color scale.}
  \item{xl,xr}{The lower left and right coordinates of the color legend in user coordinates.}
  \item{mfrowBOXPLOT}{It allows to specify the boxplot panel. It is a vector with two numbers, for example c(2,5) which means that the boxplots are put in 2 rows and 5 columns.}
  \item{mfrowMAP}{It allows to specify the map panel. It is a vector with two numbers, for example c(3,2) which means that the map panel are put in 3 rows and 2 columns.}
  \item{main}{Main title of the map with the percentage of ignorance/poor cells.}
  \item{cexCM}{Size of the points in the maps of the climate clusters.}
  \item{legpos}{Legend position with the number of the cluster in the maps of the climate clusters.}
  \item{jpg}{If TRUE the plots are exported to jpg files instead of using the windows device.}
  \item{filejpg}{Name of the jpg file.}
  \item{dec}{CSV FILE. It defines if the comma "," is used as decimal separator or the dot ".".}
}
\value{
It is obtained:

1. A TXT file with the VIF (if the argument \emph{VIF=TRUE}), the correlations between variables, the Kaiser-Meyer-Olkin (KMO) test, the Bartlett sphericity test and the results of the Principal Components or Correspondence analyses. The file is called by default "Output.TXT".

2. A CSV FILE with the coordinates for each case of the Principal Components or Correspondence analyses. The file is called by default "Cat loadings.CSV".

3. A CSV FILE with the descriptive statistics of each variable for each of the clusters obtained in the Dendrogram. The file is called by default "Descriptive statistics of clusters.CSV".

4. A CSV FILE with the original data of the variables and the cluster to which each case belongs. The file is called by default "Original data and cluster number.CSV".

5. A CSV FILE with the coordinates of the variables in the Linear Discriminant Analysis plot. The file is called by default "Var loadings-Linear.csv"

6. A CSV FILE with the coordinates of the categories in the Linear Discriminant Analysis plot. The file is called by default "Cat loadings-Linear.csv".

7. A CSV FILE with the predictions table using the cross-validation of Linear Discriminant Analysis. The file is called by default "Table cross-validation-Linear.csv".

8. A CSV FILE with the group to which each case belongs and the prediction of the Discriminant Analysis using the cross-validation of the Linear Discriminant Analysis. The file is called by default "Cases cross-validation-Linear.csv".

9. A CSV file with the predictions table using the cross-validation of the Quadratic Discriminant Analysis. The file is called by default "Table cross-validation-Quadratic.csv".

10. A CSV file with the group to which each case belongs and the prediction of the Discriminant Analysis using the cross-validation of the Quadratic Discriminant Analysis. The file is called by default "Cases cross-validation-Quadratic.csv".

11. A CSV file with the obtained probabilities of comparing all the variables among all the clusters with the U Mann-Whitney test. The file is called by default "U Mann-Whitney test.csv".

12. A CSV file called by default "Priorization.csv", with the clusters of the different climate zones arranged from the cluster with the higher percentage of cells catalogued as ignorance/poor (lowest quality survey) to the cluster with the lower percentage of cells catalogued as ignorance/poor (highest quality survey).

13. A scatterplot of the Principal Components or Correspondence analyses.

14. A Dendrogram grouping by clusters according to the scores of the Principal Components or Correspondence analyses.

15. A graphic panel with a boxplot for each variable comparing the values of these variables between each of the clusters obtained in the Dendrogram.

16. A Graph of the Discriminant Analysis showing the influence of the variables on the discriminant axis I, differentiating the different clusters.

17. A graph of the Discriminant Analysis showing the scores of the discriminant axes I and II, differentiating the different clusters.

18. If the argument \emph{BUBBLE=TRUE}, a bubble chart with the number of variables that are statistically different between clusters.

19. The maps with the climate zones.

20. The map with the percentage of ignorance/poor cells.

}
\details{
The aim of this algorithm is to identify different climate zones and to estimate survey quality of these zones. The climate zones are classified from those with higher percentage of ignorance/poor survey quality cells (poor surveyed climate zones), to those with the lower percentage of ignorance/poor survey quality cells (better surveyed climate zones). This function uses the algorithm of the VARSEDIG function (Guisande et al., 2016; 2019; Guisande, 2019, which is briefly summarized as follows.

This function only works if there are ASC files with the environmental variables in the working directory.

In the first step of the algorithm, a Principal Components analysis is performed, being the cases the different cells and the variables the environmental variables. The aim is to determine the environmental variables responsible for the variability observed among the cells.

To detect the potential groups being formed in the Principal Components analysis, a Dendrogram is applied to the scores obtained from the axes that absorb a greater variance. By default, the axes that absorb 80\% of the variability are chosen, but this value can be modified by the user.

Subsequently, a Discriminant Analysis is carried out to determine if the clusters that have been generated are well discriminated, that is, to determine the number of correctly identified cases in each cluster.

Next, a U Mann-Whitney test is performed to determine if there are significant differences in the variables between the clusters.

The idea of this function is to find the largest possible number of clusters with the highest discrimination percentage. To do this, the user should perform tests modifying the cut-off threshold by specifying the cumulative variance percentage to determine how many axes are selected from the Main Components (by default \emph{por=80}) and the variables to be included, eliminating those that are not correlated and are not useful in the Principal Components analyses, as well as those that have little discrimination power in the Discriminant Analysis.

Once the different climate zones have been identified, with the file Estimators.CSV obtained from the function KnowB (Lobo et al., 2018; Guisande & Lobo, 2019), with the values of the slope, completeness and the ratio between the number of records and the observed species (R/S) for different cells, it is estimated the percentage of ignorance/poor cells in each climate zone. Ignorance cell are those for which was not possible to estimate the slope, completeness and/or R/S, because there are not records or the number of records in small. Poor quality surveys cells are those that slope > 0.3, completeness < 50% or R/S ratio < 3, although these thresholds may be modified by the user in the function. With this information, one of the outputs of the function is a map showing all climate zones identified accordantly with the percentage of ignorance/poor cells in each climate zone.

\strong{FUNCTIONS}

This function uses the algorithm of the VIDTAX function of the VARSEDIG package (Guisande et al., 2016; 2019). The Principal Components Analysis was performed with the \link[stats]{prcomp} function of the stats package. The \link[usdm]{vif} function of the usdm package was used for the calculation of VIF (Naimi et al., 2014; Naimi, 2017). To perform the \emph{biplot} graph the \link[car]{scatterplot} function of the car package was used (Fox et al., 2018). The arrows are depicted with the function \link[IDPmisc]{Arrows} of the package IDPmisc (Locher & Ruckstuhl, 2014). The convex hull is estimated with the function \link[grDevices]{chull} of the package grDevices. KMO test was performed with the function \link[psych]{KMO} of the package psych (Revelle, 2018).  The U Mann-Whitney test is performed with the \emph{wilcox.test} function of the base stats package. The comparison between clusters with the VARSEDIG algorithm is done with the VARSEDIG function  (Guisande et al., 2016). The Linear Discriminant Analysis was performed with the functions \link[candisc]{candisc} of the candisc package (Friendly, 2007; Friendly & Fox, 2017) and \link[MASS]{lda} of the MASS package (Venables & Ripley, 2002; Ripley et al., 2018). The Quadratic Discriminant Analysis was performed with the function \link[MASS]{qda} of the MASS package (Venables & Ripley, 2002; Ripley et al., 2018). The graph with one dimension in the Discriminant analysis was performed with the function \link[candisc]{plot.cancor} of the candisc package (Friendly, 2007; Friendly & Fox, 2017).

\strong{EXAMPLE}

We used the file Estimators obtained with the function \link[KnowBR]{KnowB} using a database that includes 121,709 records of freshwater fishes obtained from Pelayo-Villamil et al. (2015), with the Clench estimator and a cell resolution of 5'.

In the working directory, there are the environmental variables BIO1, BO2, BIO4, BIO8, BIO12, BIO14, BIO15, BIO18 and BIO19 of the WorldClim data set (Hijmans et al., 2005), as ASC raster files. 

As the argument \emph{VIF=FALSE}, there is not information about VIF, and the first statistic obtained is the KMO test, which tells us if the variables are adequate for the Principal Components. The value must be greater than 0.5. Therefore, all variables that do not have a value greater than 0.5, could be eliminated from the analysis. In the case that the value is exactly 0.5, it means that it is not possible to estimate the KMO. The minimum value was 0.73 for BIO12, so all variables are adequate for the Principal Components.

The next statistic that appears is Bartlett's test of sphericity, which tests whether the correlation matrix is an identity matrix, which would indicate that the factor model is inappropriate. A value \emph{p} of the contrast smaller than the level of significance allows rejecting the hypothesis and concluding that there is correlation. Therefore, for the Principal Components analysis to be valid, the probability must be less than 0.05, as it is in this case, with a p<0.001.

The first figure is the Principal Components. The first axis accounts for 62.9\%, the second for 16.5\% and the third for 7\% of the variance observed. The first three axes explain 86.5\% of the variance. Since the default value of \emph{por=80} was selected, these three Principal Component axes are selected.

In the Dendrogram there were 22 clusters statistically different, because the argument \emph{k=NULL} in the script (default option), which means that the algorithm finds the maximum number of statistically different climate zones. The results obtained in the files U Mann-Whitney test.CSV and Descriptive statistics of clusters.CSV show that there are statistical differences among all 22 clusters in at least one variable (U Mann-Whitney test, p<0.001).

Other plots and statistics are obtained, which are fully explained in the manuscript that describes the algorithm VIDTAXA (Guisande et al., 2019).

The final plot is the map with the percentage of ignorance/poor cells (see map below). It is clear that the area with a better survey quality is in South-East of United Kingdom.
\tabular{c}{
\if{html}{\figure{F9.jpg}{options: width="50"}}
\if{latex}{\figure{F9.jpg}{options: width=2.7in}}
}
}
\references{
Fox, J., Weisberg, S., Adler, D., Bates, D., Baud-Bovy, G., Ellison, S., Firth, D., Friendly, M., Gorjanc, G., Graves, S., Heiberger, R., Laboissiere, R., Monette, G., Murdoch, D., Nilsson, H., Ogle, D., Ripley, B., Venables, W. & Zeileis, A. (2018) Companion to Applied Regression. R package version 3.0-0. Available at: \url{https://CRAN.R-project.org/package=car}.

Friendly, M. & Fox, J. (2017) Visualizing Generalized Canonical Discriminant and Canonical Correlation Analysis. R package version 0.8-0. Available at: \url{https://CRAN.R-project.org/package=candisc}.

Friendly, M. (2007). HE plots for Multivariate General Linear Models. \emph{Journal of Computational and Graphical Statistics}, 16: 421-444.

Guisande, C., Vari, R.P., Heine, J., García-Roselló, E., González-Dacosta, J., Pérez-Schofield, B.J., González-Vilas, L. & Pelayo-Villamil, P. (2016) VARSEDIG: an algorithm for morphometric characters selection and statistical validation in morphological taxonomy. \emph{Zootaxa}, 4162: 571-580.

Guisande, C., Rueda-Quecho, A.J., Rangel-Silva, F.A., Heine, J., García-Roselló, E., González-Dacosta, J. & Pelayo-Villamil, P. (2019) VIDTAXA: an algorithm for the identification of statistically different groups based on variability obtained in factorial analyses. \emph{Ecological Informatics}, 46: 62-68.

Hijmans, R.J., Cameron, S.E., Parra, J.L., Jones, P.G. & Jarvis, A. (2005) Very high resolution interpolated climate surfaces for global land areas. \emph{International Journal of Climatology}, 25: 1965-1978.

Locher, R. & Ruckstuhl, A. (2014) Utilities of Institute of Data Analyses and Process Design. R package version 1.1.17. Available at: \url{https://CRAN.R-project.org/package=IDPmisc}.

Naimi, B. (2017) Uncertainty analysis for species distribution models. R package version 1.1-18. Available at: \url{https://CRAN.R-project.org/package=usdm}.

Naimi, B., Hamm, N.A.S., Groen, T.A., Skidmore, A.K., & Toxopeus, A.G. (2014) Where is positional uncertainty a problem for species distribution modelling? \emph{Ecography}, 37: 191-203.

Pelayo-Villamil, P., Guisande, C., Vari, R.P., Manjarrés-Hernández, A., García-Roselló, E., González-Dacosta, J., Heine, J., González-Vilas, L., Patti, B., Quinci, E.M., Jiménez, L.F., Granado-Lorencio, C., Tedesco, P.A., Lobo, J.M. (2015) Global diversity patterns of freshwater fishes-Potential victims of their own success. \emph{Diversity and Distributions}, 21: 345-356.

Revelle,W. (2018) Procedures for Psychological, Psychometric, and Personality Research. R package version 1.8.4. Available at: \url{https://CRAN.R-project.org/package=psych}.

Ripley, B., Venables, B., Bates, D.M., Hornik, K., Gebhardt, A. & Firth, D. (2018) Support Functions and Datasets for Venables and Ripley's MASS. R package version 7.3-50. Available at: \url{https://CRAN.R-project.org/package=MASS}.

Venables, W.N. & Ripley, B.D. (2002) \emph{Modern Applied Statistics with S}. Springer, fourth edition, New York. \url{https://www.stats.ox.ac.uk/pub/MASS4/}.

}
\examples{
\dontrun{
####This script only works if there are ASC files, with
####environmental variables, in the working directory

data(FishIrelandUK)

data(adworld)

SurveyQCZ(data=FishIrelandUK, maxLon=3, mfrowBOXPLOT=c(3,3), cexCM=0.2)
}
}
\keyword{SurveyQCZ}


