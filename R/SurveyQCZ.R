SurveyQCZ<-function(data, Longitude="Longitude", Latitude="Latitude",  cell=NULL, hull=TRUE, Area="World", shape=NULL, shapenames=NULL, aprox=TRUE, VIDTAXA=NULL, por=80,
k=NULL, VIF=FALSE, VARSEDIG=FALSE,  BUBBLE=FALSE, variables = c("Slope","Completeness","Ratio"), completeness=c(50,90), slope=c(0.02,0.3), ratio=c(3,15),  minLon=NULL,
maxLon=NULL, minLat=NULL, maxLat=NULL, xlab="Longitude", ylab="Latitude", colscale=c("#C8FFFFFF","#64FFFFFF","#00FFFFFF","#64FF64FF",
"#C8FF00FF","#FFFF00FF","#FFC800FF","#FF6400FF","#FF0000FF"),  colcon="transparent", 
breaks=10, ndigits=0, xl=0, xr=0,  mfrowBOXPLOT=NULL, mfrowMAP=NULL, main="Percentage of ignorance/poor\n cells in each cluster", cexCM=0.5, legpos="bottomleft", jpg=FALSE, filejpg="Survey Quality CZ.jpg", dec=","){

data<-na.exclude(data)

#####Checking data required

if(exists("adworld")==FALSE){
adworld<-1
stop("It is necessary to load data(adworld)")
}

if(Area!="World" & exists("adworld1")==FALSE){
stop("It is necessary to use RWizard and replace data(adworld) by @_Build_AdWorld_, for using administative areas")
}

if(Area!="World" & exists("adworld2")==FALSE){
stop("It is necessary to use RWizard and replace data(adworld) by @_Build_AdWorld_, for using administative areas")
}

if(!is.null(shape) & is.null(shapenames)){
stop("It is necessary to specify in the argument 'shapenames' the variable with the names of the polygons in the shape")
}



ASC<-TRUE



ZZ<-matrix(c("","","",""), nrow=2)

if(ASC==TRUE){

####Adding variable to file Estimators from ASC files

ZZ[1,1]<-"ADDING VARIABLES TO THE FILE ESTIMATORS FROM ASC FILES"

lista<-list.files(pattern=".ASC$") 
if(length(lista)==0){
lista<-list.files(pattern=".asc$") 
}

leng<-length(lista)

if(leng==0){
stop("There are not ASC files in the working directory")
}

vars<-sapply(strsplit(lista, split='.', fixed=TRUE), function(x) (x[1]))

####Estimation of convex hull and grid

if(hull==TRUE){
ch <- chull(data[,Longitude], data[, Latitude])
coords <- data[c(ch, ch[1]), c(Longitude,Latitude)]  
maxpoLon<-max(data[, Longitude], na.rm=TRUE)
minpoLon<-min(data[, Longitude], na.rm=TRUE)
maxpoLat<-max(data[, Latitude], na.rm=TRUE)
minpoLat<-min(data[, Latitude], na.rm=TRUE)
}
else{

if(is.null(shape)){
extentS<-adworld

if(exists("adworld1")==TRUE){
if(!is.null(adworld1)){
extentS<-adworld1
}
}

if(exists("adworld2")==TRUE){
if(!is.null(adworld2)){
extentS<-adworld2
}
}

maxpoLon<-max(extentS[, "Lon"], na.rm=TRUE)
minpoLon<-min(extentS[, "Lon"], na.rm=TRUE)
maxpoLat<-max(extentS[, "Lat"], na.rm=TRUE)
minpoLat<-min(extentS[, "Lat"], na.rm=TRUE)

}####end NULL shape


if(!is.null(shape)){
if(inherits(shape, "list")){
extentS<-shape[[1]]
lsh<-length(shape)
if(lsh>1){
ss<-seq(2,lsh)
hh<-as.character(shape[ss])
extentS<-eval(parse(text=paste("subset(data,",noquote(shapenames), " %in% hh)", sep="")))
}
}
else{
extentS<-shape
if(inherits(extentS,"character")){
extentS<-eval(parse(text=paste(".GlobalEnv$", extentS, sep="")))
}

}

numero<-length(extentS)
Areas<-eval(parse(text=paste("extentS$",noquote(shapenames),sep="")))

for(z in 1:numero){
polygontemp<-as.data.frame(extentS@polygons[[z]]@Polygons[[1]]@coords)

if(z==1) {
prueba<-polygontemp
}
else{
prueba<-rbind(prueba,c(NA,NA), polygontemp)
}
}

names(prueba)<-c("Lon","Lat")

maxpoLon<-max(prueba[, "Lon"], na.rm=TRUE)
minpoLon<-min(prueba[, "Lon"], na.rm=TRUE)
maxpoLat<-max(prueba[, "Lat"], na.rm=TRUE)
minpoLat<-min(prueba[, "Lat"], na.rm=TRUE)

}####end not NULL shape

}

if(exists("adworld1")==FALSE){
adworld1<-1
}

if(exists("adworld2")==FALSE){
adworld2<-1
}

Lat<-data[,Latitude]
Lon<-data[,Longitude]

nvar<-1


for(hh in 1:leng){

var<-raster::raster(lista[hh])
variable<-vars[hh]

if(hh>1) gg<-1 else gg<-2

ZZ[gg,1]<-paste(variable, ", variable",hh, "of", leng,"variables")
write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)

ZZ<-matrix(c("","","",""), nrow=2)

####ASC to matrix

if(inherits(var, "RasterLayer")){




###Format of matrix

reso<-raster::res(var)

if(is.null(cell)){
cell<-reso[1]*60
}


if(round(raster::xmin(var))==-180 & round(raster::ymin(var))==-90 & round(raster::xmax(var))==180 & round(raster::ymax(var))==90){
m1<-raster::as.matrix(var)
dimm<-dim(m1)
long<-seq(from=(-180+360/dimm[2]), to = 180 , by = 360/dimm[2])
m1<-rbind(long,m1)

lat<-seq(from=(90-180/dimm[1]), to = -90 , by = -180/dimm[1])
lat<-c(0,lat)
var<-cbind(lat,m1, deparse.level=0)
}
else{

reso<-raster::res(var)

if(is.null(cell)){
cell<-reso[1]*60
}

r1<-raster::raster(xmn=-180, xmx=180, ymn=-90, ymx=90, resolution=reso)

var<-raster::resample(var,r1)

m1<-raster::as.matrix(var)
dimm<-dim(m1)

long<-seq(from=(raster::xmin(var)+(raster::xmax(var)-raster::xmin(var))/dimm[2]), to = raster::xmax(var) , by = (raster::xmax(var)-raster::xmin(var))/dimm[2])
m1<-rbind(long,m1)

lat<-seq(from=(raster::ymax(var)+(raster::ymin(var)-raster::ymax(var))/dimm[1]), to = raster::ymin(var) , by = (raster::ymin(var)-raster::ymax(var))/dimm[1])
lat<-c(0,lat)
var<-cbind(lat,m1, deparse.level=0)

}


#####Rescale the matrix


if(reso[1]<cell/60){

va<-60/cell*180
var<-var[-1,-1]
dimf<-dim(var)
co<-dimf[1]/va-1
matemp<-matrix(0,60/cell*180, 60/cell*360)


fila<-0
col<-0
for(rr in seq(from=1, to=dimf[1], by=co+1)){
fila<-fila+1
for(cc in seq(from=1, to=dimf[2], by=co+1)){
col<-col+1
valorr<-var[rr:(rr+co),cc:(cc+co)]
matemp[fila,col]<-mean(valorr[valorr!=-9999], na.rm=TRUE) 
}
col<-0
}

long<-seq(-180+cell/60, 180, cell/60)
lati<-c(0,seq(90-cell/60, -90, -cell/60))

matemp<-rbind(long,matemp)
var<-cbind(lati,matemp)

}###End loop rescale



####Selection of grid points inside polygons

if(hull==TRUE){

if(hh==1){

seq1<-(seq(min(Lon, na.rm=TRUE)-abs(0.3*min(Lon, na.rm=TRUE)/100), max(Lon, na.rm=TRUE)+abs(0.3*max(Lon, na.rm=TRUE)/100),cell/60))
seq2<-(seq(min(Lat, na.rm=TRUE)-abs(0.2*min(Lat, na.rm=TRUE)/100), max(Lat, na.rm=TRUE)+abs(0.2*max(Lat, na.rm=TRUE)/100),cell/60))
data.grid<-expand.grid(Lon=seq1, Lat=seq2)


fila<-which((data.grid[,"Lon"] >= minpoLon) &  (data.grid[,"Lon"] <= maxpoLon) & (data.grid[,"Lat"] >= minpoLat) &  (data.grid[,"Lat"] <= maxpoLat))
borrar<-data.grid[fila,]
log<-mgcv::in.out(as.matrix(coords[,c(Longitude, Latitude)]),as.matrix(borrar[, c("Lon", "Lat")]))
mm<-cbind(borrar,log)
dhh<-dim(mm)
mm<-mm[mm[,dhh[2]],]
mm<-mm[,c(-dhh[2])]
}
}###end hull TRUE



##Format Lon, Lat and variable

Lon<-var[1,-1]
lenLon<-length(Lon)


Lat<-var[-1,1]
lenLat<-length(Lat)



Longi<-matrix(t(replicate(lenLat,Lon)),ncol=1)
Lati<-rep(Lat,lenLon)
var<-matrix(var[-1,-1],ncol=1)


matriz<-data.frame(Longi,Lati, var)
names(matriz)<-c("Longitude","Latitude",variable)

}

fila<-which((matriz[,"Longitude"] >= minpoLon) &  (matriz[,"Longitude"] <= maxpoLon) & (matriz[,"Latitude"] >= minpoLat) &  (matriz[,"Latitude"] <= maxpoLat))
matriz<-matriz[fila,]


if(nvar==1){
envar<-matriz
nvar<-2
}
else{
nombres<-names(envar)
envar<-cbind(envar,matriz[,3])
names(envar)<-c(nombres,variable)
}



}#####End for add variables

}####End ASC TRUE




#####Selection of the var coordinates nearest to grid coordinates

ZZ[1,1]<-"SELECTION OF VAR COORDINATES INSIDE THE POLYGONS"
write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)

if(hull==TRUE){
pos<-0

dimr<-dim(mm)
for(z in 1:dimr[1]){
valor<-(mm[z,"Lon"]-envar[,"Longitude"])^2+(mm[z,"Lat"]-envar[,"Latitude"])^2
pos1<-which(valor==min(valor,na.rm=TRUE))
if(is.na(envar[pos1,3]) & aprox==TRUE){
post<-pos1-1
if(is.na(envar[post,3])){
post<-pos1+1
}
pos1<-post
}
pos<-append(pos,pos1)
}
pos<-pos[-1]
envar<-unique(na.exclude(envar[pos,]))
dimv<-dim(envar)
ID<-seq(1,dimv[1],1)
envar<-cbind(ID,envar)
}

else{

if(is.null(shape)){
if(any(is.na(extentS))==TRUE){
dim<-dim(extentS)
L<-extentS[,"Lon"]
nas<-which(is.na(L))
lenA<-length(nas)
}
else{
polygontemp<-extentS
lenA<-0
}

zh<-(-1)

fila<-which((envar[,"Longitude"] >= minpoLon) &  (envar[,"Longitude"] <= maxpoLon) & (envar[,"Latitude"] >= minpoLat) &  (envar[,"Latitude"] <= maxpoLat))
borrar<-envar[fila,]


while(zh<lenA){

if(any(is.na(extentS))==TRUE){

if(zh==-1){
polygontemp<-extentS[1:(nas[zh+2]-1),]
}
if(zh==(lenA-1)){
polygontemp<-extentS[(nas[zh+1]+1):dim[1],]
}
if(zh>-1 & zh<(lenA-1)){
polygontemp<-extentS[(nas[zh+1]+1):(nas[zh+2]-1),]
}
}

log<-mgcv::in.out(as.matrix(polygontemp[,c("Lon", "Lat")]),as.matrix(borrar[, c(Longitude, Latitude)]))

mm<-cbind(borrar,log)
dhh<-dim(mm)
mm<-mm[mm[,dhh[2]],]
mm<-mm[,c(-dhh[2])]
if(zh==(-1)) prueba<-mm else prueba<-rbind(prueba,mm)


zh<-zh+1


}#End while

envar<-prueba
remove(prueba)
envar<-unique(na.exclude(envar))
dimv<-dim(envar)
ID<-seq(1,dimv[1],1)
envar<-cbind(ID,envar)

}####end NULL shape

if(!is.null(shape)){

fila<-which((envar[,"Longitude"] >= minpoLon) &  (envar[,"Longitude"] <= maxpoLon) & (envar[,"Latitude"] >= minpoLat) &  (envar[,"Latitude"] <= maxpoLat))
borrar<-envar[fila,]

for(z in 1:numero){
polygontemp<-as.data.frame(extentS@polygons[[z]]@Polygons[[1]]@coords)
log<-mgcv::in.out(as.matrix(polygontemp[,1:2]),as.matrix(borrar[, c(Longitude, Latitude)]))

mm<-cbind(borrar,log)
dhh<-dim(mm)
mm<-mm[mm[,dhh[2]],]
mm<-mm[,c(-dhh[2])]
if(z==1) prueba<-mm else prueba<-rbind(prueba,mm)
}

envar<-prueba
remove(prueba)
envar<-unique(na.exclude(envar))
dimv<-dim(envar)
ID<-seq(1,dimv[1],1)
envar<-cbind(ID,envar)

}#####end not NULL shape

}


#####Quality

dimda<-dim(data)
Quality<-rep("Fair",dimda[1])
data<-cbind(data,Quality)
data[,"Quality"]<-as.character(data[,"Quality"])
datosL<-subset(data,data[,variables[1]]<slope[1] & data[,variables[2]]>completeness[2] & data[,variables[3]]>ratio[2])
dimv<-dim(datosL)
pos<-0
for(z in 1:dimv[1]){
valor<-(datosL[z,Longitude]-data[,Longitude])^2+(datosL[z,Latitude]-data[,Latitude])^2
pos1<-which(valor==min(valor,na.rm=TRUE))
pos<-append(pos,pos1)
}
pos<-pos[-1]
data[pos,"Quality"]<-"Good"


datosL<-subset(data,data[,variables[1]]>slope[2] & data[,variables[2]]<completeness[1] & data[,variables[3]]<ratio[1])
dimv<-dim(datosL)
pos<-0
for(z in 1:dimv[1]){
valor<-(datosL[z,Longitude]-data[,Longitude])^2+(datosL[z,Latitude]-data[,Latitude])^2
pos1<-which(valor==min(valor,na.rm=TRUE))
pos<-append(pos,pos1)
}
pos<-pos[-1]
data[pos,"Quality"]<-"Poor"



###VIDTAXA

ZZ[1,1]<-"PERFORMING VIDTAXA ALGORITHM"
write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)

if(!is.null(VIDTAXA)){
exe<-paste("VARSEDIG::VIDTAXA(",toString(x=VIDTAXA), ")")
eval(parse(text=exe))
}
else{
exe<-paste("VARSEDIG::VIDTAXA(","data=envar,", "var=vars,", "labels='ID',", "VARSEDIG=VARSEDIG,",
"BUBBLE=BUBBLE,", "por=por,", "k=k,", "VIF=VIF,", "mfrowBOXPLOT=mfrowBOXPLOT,", "dec=dec",  ")")
eval(parse(text=exe))
}



###Adding quality

if(dec=="."){
dataf<-read.csv("Original data and cluster number.csv",header=TRUE,encoding="latin1")
}
else{
dataf<-read.csv2("Original data and cluster number.csv",header=TRUE,encoding="latin1")
}

dataL<-merge(envar[,c("ID","Longitude", "Latitude")], dataf, by="ID")


dimv<-dim(dataL)

pos<-0

dimd<-dim(data)
for(z in 1:dimd[1]){
valor<-(data[z,Longitude]-dataL[,"Longitude"])^2+(data[z,Latitude]-dataL[,"Latitude"])^2
done<-0
while(done==0 | length(valor)==0){
pos1<-which(valor==min(valor,na.rm=TRUE))
if(any(pos==pos1)==TRUE){
valor<-valor[-pos1]
}
else{
done<-1
}
}###end while
pos<-append(pos,pos1)
}

pos<-pos[-1]


Slope<-rep(NA, dimv[1])
Completeness<-rep(NA, dimv[1])
Ratio<-rep(NA, dimv[1])
Quality<-rep(NA, dimv[1])


for(z in 1:dimv[1]){
Slope[pos[z]]<-data[z,variables[1]]
Completeness[pos[z]]<-data[z,variables[2]]
Ratio[pos[z]]<-data[z,variables[3]]
Quality[pos[z]]<-as.character(data[z,"Quality"])
}

dataf<-data.frame(dataL,Slope,Completeness, Ratio, Quality)


###Order climatic areas by quality

dimda<-dim(dataf)


tabla<-table(dataf[,"Quality"], dataf[,"Cluster"], exclude=NULL)


names<-rownames(tabla)

if(!is.na(any(names=="Poor"))==TRUE & !is.na(any(is.na(names)))==TRUE){
p<-which(names=="Poor")
na<-which(is.na(names))
tabla[p,]<-tabla[p,]+tabla[na,]
tabla<-tabla[-na,]
names<-rownames(tabla)
names<-replace(names, names=="Poor", "IgnPoor")
rownames(tabla)<-names
}

if(!is.na(any(is.na(names)))==TRUE){
na<-which(is.na(names))
tabla[na,]<-tabla[na,]
names<-rownames(tabla)
names<-replace(names, names=="Poor", "IgnPoor")
rownames(tabla)<-names
}

if(!is.na(any(names=="Poor"))==TRUE){
p<-which(names=="Poor")
tabla[p,]<-tabla[p,]
names<-rownames(tabla)
names<-replace(names, names=="Poor", "IgnPoor")
rownames(tabla)<-names
}


tadim<-dim(tabla)
suma<-apply(X = tabla , MARGIN = 2 , FUN = sum , na.rm=TRUE)
for(ss in 1:tadim[2]){
tabla[,ss]<-tabla[,ss]*100/suma[ss]
}

IgnPoorPercentage<-rep(100,dimda[1])

dataf<-cbind(dataf,IgnPoorPercentage)

Poor<-tabla["IgnPoor",]

clusters<-as.numeric(names(Poor))

len<-length(clusters)

for(ff in 1:len){
pos<-which(dataf[,"Cluster"]==clusters[ff])
dataf[pos,"IgnPoorPercentage"]<-Poor[ff]
}


dataf<-dataf[order(dataf[,"IgnPoorPercentage"] , decreasing = TRUE),]



###Saving file

if(dec=="."){
write.csv(dataf,"Priorization.csv", row.names=FALSE)
}
else{
write.csv2(dataf,"Priorization.csv", row.names=FALSE)
}


####Plotting map with climate clusters

inc<-0.05

if(is.null(minLon)){
if(min(dataf$Longitude)<0) minLon<-(min(dataf$Longitude, na.rm=TRUE)+min(dataf$Longitude, na.rm=TRUE)*inc) else minLon<-(min(dataf$Longitude, na.rm=TRUE)-min(dataf$Longitude, na.rm=TRUE)*inc)
}
if(is.null(maxLon)){
if(max(dataf$Longitude)<0) maxLon<-(max(dataf$Longitude, na.rm=TRUE)-max(dataf$Longitude, na.rm=TRUE)*inc) else maxLon<-(max(dataf$Longitude, na.rm=TRUE)+max(dataf$Longitude, na.rm=TRUE)*inc)
}

if(is.null(maxLat)){
if(max(dataf$Latitude)<0) maxLat<-(max(dataf$Latitude, na.rm=TRUE)-max(dataf$Latitude, na.rm=TRUE)*inc) else maxLat<-(max(dataf$Latitude, na.rm=TRUE)+max(dataf$Latitude, na.rm=TRUE)*inc)
}

if(is.null(minLat)){
if(min(dataf$Latitude)<0) minLat<-(min(dataf$Latitude, na.rm=TRUE)+min(dataf$Latitude, na.rm=TRUE)*inc) else minLat<-(min(dataf$Latitude, na.rm=TRUE)-min(dataf$Latitude, na.rm=TRUE)*inc)
}


########### function written by Greg Snow
squishplot <- function(xlim,ylim,asp=1){
   if(length(xlim) < 2) stop('xlim must be a vector of length 2')
   if(length(ylim) < 2) stop('ylim must be a vector of length 2')

  tmp <- par(c('plt','pin','xaxs','yaxs'))

  if( tmp$xaxs == 'i' ){ # not extended axis range

        xlim <- range(xlim)
  } else { # extended range

	tmp.r <- diff(range(xlim))
	xlim <- range(xlim) + c(-1,1)*0.04*tmp.r

  }

  if( tmp$yaxs == 'i' ){ # not extended axis range

        ylim <- range(ylim)
  } else { # extended range

	tmp.r <- diff(range(ylim))
	ylim <- range(ylim) + c(-1,1)*0.04*tmp.r

  }

  tmp2 <- (ylim[2]-ylim[1])/(xlim[2]-xlim[1])

  tmp.y <- tmp$pin[1] * tmp2 * asp

  if(tmp.y < tmp$pin[2]){ # squish vertically
	par(pin=c(tmp$pin[1], tmp.y))
	par(plt=c(tmp$plt[1:2], par('plt')[3:4]))
  } else { # squish horizontally
	tmp.x <- tmp$pin[2]/tmp2/asp
	par(pin=c(tmp.x, tmp$pin[2]))
	par(plt=c(par('plt')[1:2], tmp$plt[3:4]))

  }

  return(invisible(tmp['plt']))
} # end of function
###################


if (Area=="World") {
AA<-data.frame(adworld$Lon,adworld$Lat)
}
else{
AA<-data.frame(adworld1$Lon,adworld1$Lat)
}

names(AA)<-c("Lon", "Lat")

dev.new()

Lati<-(maxLat+minLat)/2
aspe=(1/cos(Lati*pi/180)) 



cl<-length(unique(dataf[,"Cluster"]))
clus<-unique(dataf[,"Cluster"])

div<-ceiling(cl/2)
if(div>5){
div<-5
}
ini<-ceiling(cl/div)

if(is.null(mfrowMAP)){
par(mfrow=c(ini,div), mar=c(0,0,0,0))
}
else{
par(mfrow=c(mfrowMAP[1],mfrowMAP[2]), mar=c(0,0,0,0))
}
tmp<-squishplot(xlim=c(minLon,maxLon), ylim=c(minLat,maxLat), asp=aspe)
for(z in 1:cl){
plot(0,0, xlim=c(minLon,maxLon), ylim=c(minLat, maxLat), xlab=xlab, ylab=ylab, type="n")
polygon(AA$Lon,AA$Lat,col=colcon, border="black")
dati<-subset(dataf, dataf[,"Cluster"]== clus[z])
points(x=dati$Longitude, y=dati$Latitude, cex=cexCM, pch=15, col="red")
leg<-paste("C",clus[z], sep="")
legend(x = legpos , legend = leg , bty = 'n' , text.font = 2)
}



#####Matrix for maps

matemp<-matrix(-9999,60/cell*180, 60/cell*360)
long<-seq(-180+cell/60, 180, cell/60)
lati<-c(0,seq(90-cell/60, -90, -cell/60))

matemp<-rbind(long,matemp)
matemp<-cbind(lati,matemp)



dimr<-dim(dataf)
for(z in 1:dimr[1]){
valor<-(dataf[z,"Longitude"]-long)^2
col<-which(valor==min(valor,na.rm=TRUE))
valor<-(dataf[z,"Latitude"]-lati)^2
fila<-which(valor==min(valor,na.rm=TRUE))
matemp[fila,col+1]<-dataf[z,"IgnPoorPercentage"]
}


####Plotting map

dev.new()


KnowBR::MapCell(data = matemp, Area=Area, minLon=minLon, maxLon=maxLon, minLat=minLat, maxLat=maxLat, colcon=colcon, main=main,
xlab=xlab, ylab=ylab, jpg=jpg, filejpg=filejpg, ndigits=ndigits, xl=xl, xr=xr, breaks=breaks, colscale=colscale)

}

