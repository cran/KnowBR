KnowB<-function(data, matrix=TRUE, cell=60, method="accumulation", curve= "Clench", estimator=1,
cutoff=2, cutoffCompleteness= 0, cutoffSlope= 0.1, SpR=FALSE, largematrix=FALSE,
cormatrix=FALSE, Area="World", extent=TRUE, minLon, maxLon, minLat, maxLat,
colbg="#FFFFFF", colcon="#C8C8C8", colf="black", pro = TRUE, inc = 0.005, exclude = NULL,
colexc = NULL, colfexc="black", colscale=c("#FFFFFFFF", "#C8FFFFFF","#64FFFFFF","#00FFFFFF",
"#64FF64FF","#C8FF00FF", "#FFFF00FF","#FFC800FF","#FF6400FF","#FF0000FF"), legend.pos="y",
breaks=10, xl=0, xr=0, yb=0, yt=0, asp, lab = NULL, xlab = "Longitude", ylab = "Latitude",
main1="Actual richness", main2="Predicted richness", main3="Residuals", main4="Records",
main5="Completeness", main6="Slope", cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.2, 
family = "sans", font.main = 2, font.lab = 1, font.axis = 1, lwdP=0.6, lwdC=0.1, trans=c(1,1),
log=c(0,0), ndigits=0, ini=NULL, end=NULL, ComCut=TRUE, SlopeCut=TRUE, file1 = "Actual richness.csv",
file2 = "Predicted richness.csv", file3 = "Residuals.csv", file4 = "List of species.csv",
file5 = "Species per site.csv", file6 = "Estimators.csv", file7 = "Species per record.csv",
file8 = "Records.csv", file9 = "Completeness.csv", file10 = "Slope.csv", file11 = "Standard error of the estimators.csv",
na = "NA", dec = ",", row.names = FALSE, fileEncoding = "", jpg=TRUE, jpg1="Actual richness.jpg",
jpg2="Predicted richness.jpg", jpg3="Residuals.jpg", jpg4="Correlation matrix.jpg",
jpg5="Records.jpg",  jpg6="Completeness.jpg", jpg7="Slope.jpg",cex=1.5, pch=15,
cex.labels=1.5, pchcol="red", ask=FALSE){

options(warn=-1)

if(jpg==FALSE) par(ask=ask) else yuret<-1


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

if(exists("adworld1")==FALSE){
adworld1<-1
}

if(exists("adworld2")==FALSE){
adworld2<-1
}

####End Checking

if(matrix==FALSE){
data[is.na(data)]<-0
}
x<-na.exclude(data)

if(matrix==FALSE){
matrix<-TRUE
xLo<-x[,1]
xLa<-x[,2]
dimg<-dim(x[,c(-1,-2)])
replicas<-rep(dimg[1], dimg[2])
x2<-matrix(as.matrix(x[,c(-1,-2)]), ncol=1)
headers<-names(x[,c(-1,-2)])
sps<-rep(x=headers,times=replicas)
xLo<-rep(x=xLo,times=dimg[2])
xLa<-rep(x=xLa,times=dimg[2])
x<-data.frame(sps,xLo,xLa,x2)
names(x)<-c("Species","Longitude","Latitude", "Counts")
x<-x[x$Counts>0,]
x<-x[(x$Longitude!=0 & x$Latitude!=0),]
}


if(matrix==TRUE){
if(method=="accumulation" | method=="incidence"){
b<-x[,4]
x<-x[rep(1:nrow(x[,1:3]), b), ] 
x[,4]<-1
}
}


values1<-data.frame(1,2,3,4,5,6,7,8,9,10)
values2<-data.frame(1,2,3,4,5,6,7)
values3<-data.frame(1,2,3,4,5,6,7,8,9,10,11,12,13)
sevalues1<-data.frame(1,2,3,4,5,6,7)
sevalues2<-data.frame(1,2,3,4,5,6)
sevalues3<-data.frame(1,2,3,4,5,6,7,8)
salA<-data.frame(c(2,5,6,1,0,6,5,8,7,4,9,8), nrow=4)
sal<-data.frame(c(2,5,6,1,0,6,5,8,7,4,9,8), nrow=4)
cu2<-NA;cu3<-NA;cu4<-NA;cu5<-NA
sp2<-NA;sp3<-NA;sp4<-NA;sp5<-NA
serandom<-NA;seexact<-NA;secoleman<-NA;serarefaction<-NA


if(matrix==TRUE){
re<-dim(x)

Records<-seq(1,re[1],1)

temp<-cbind(x,Records)
div<-x[,2]*cos(180*3.1416/180)+x[,3]*sin(90*3.1416/180)+(x[,2]+x[,3])*x[,2]+(x[,2]-x[,3])*x[,3]+
(x[,2]+x[,3])*x[,3]+(x[,2]-x[,3])*x[,2]+(x[,2]-x[,3])*x[,2]*x[,3]+(x[,2]+x[,3])*x[,2]*x[,3]
dt1<-cbind(x,div)
dt1<-dt1[order(dt1[,5]), ]




pp1<-subset(temp, !duplicated(temp[,1]))
dimtemp<-dim(temp)
dimpp1<-dim(pp1)

elements<-dimtemp[1]*dimpp1[1]
rm(pp1)

elements[is.na(elements)]<-0

if(elements>2000000000) elements<-0 else elements<-elements

if(SpR==FALSE) elements<-0 else elements<-elements

if(elements==0){

datosac<-x
datosf<-x


if(largematrix==FALSE){
}
else{
pp1<-subset(temp[,c(1,2,3)], !duplicated(temp[,1]))

pp1<-pp1[order(pp1[,1]), ]

pp2<-t(pp1)

coluSp<-dim(pp2)



datos3<-aggregate(x[,4],by=list(x[,1]),mean)
datos3<-datos3[order(datos3[,1]), ]
d3<-dim(datos3)

for (t in 1:d3[1]){
sele<-subset(temp,temp[,1] %in% datos3[t,1])
dimsele<-dim(sele)
matr1<-matrix(0, nrow=dimsele[1], ncol=coluSp[2]+2)
matr1[,1]<-sele[,2]
matr1[,2]<-sele[,3]
matr1[,t+2]<-1
if(t==1){
colnames(matr1)<-c("Longitude","Latitude",pp2[1,])
write.table(matr1,"Species per record.txt", row.names=FALSE)
}
else{
write.table(matr1,"Species per record.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
}
}
rm(datos3)
}




}#Big Matrix


else{
datosf<-with(dt1, table(dt1[,5],dt1[,1]))
datosac<-with(temp, table(temp[,5],temp[,1]))
}

if(elements==0){
}
else{
if(dec=="."){
write.csv(x=datosf, file = file5, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv(x=datosac, file = file7, fileEncoding = fileEncoding,
row.names=row.names,na=na)
}
else{
write.csv2(x=datosf, file = file5, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv2(x=datosac, file = file7, fileEncoding = fileEncoding,
row.names=row.names,na=na)
}



if(dec=="."){
datosf<-read.csv(file5, header=T, check.names=FALSE)
datosac<-read.csv(file7, header=T, check.names=FALSE)
}
else{
datosf<-read.csv2(file5, header=T, check.names=FALSE)
datosac<-read.csv2(file7, header=T, check.names=FALSE)
}
}


if(elements==0){
datosf<-x
datosac<-x
}
else{
datosac<-cbind(temp[,2:3],datosac)
datos2<-dt1[,2:3]
datos2<-subset(datos2, !duplicated(datos2)) 
datosf<-cbind(datos2,datosf)
datosf<-datosf[order(datosf[,1], datosf[,2]), ]
}





datosf[is.na(datosf)]<-0



species<-aggregate(x[,4],by=list(x[,1]),FUN=sum,na.rm=TRUE)
species1<-aggregate(x[,4],by=list(x[,1]),FUN=mean,na.rm=TRUE)
species<-cbind(species,species1[,2])
colnames(species)<-c("Species","Sum", "Mean")



}
else{
datosf<-x
datosac<-replace(x[,-c(1,2)], x[,-c(1,2)]>1,1)
datosac<-cbind(x[,c(1,2)],datosac)

Sum<-colSums(x[,-c(1,2)])
Mean<-colMeans(x[,-c(1,2)])
species<-rbind(Sum,Mean)
colnames(species)<-colnames(x[,-c(1,2)])
species<-cbind(c("Sum","Mean"),species)

}

if(matrix==FALSE) elements<-1 else elements<-elements

if(elements==0){
}
else{
if(dec=="."){
write.csv(x=datosf, file = file5, fileEncoding = fileEncoding,
row.names=row.names, na=na)
write.csv(x=datosac, file = file7, fileEncoding = fileEncoding,
row.names=row.names, na=na)
}
else{
write.csv2(x=datosf, file = file5, fileEncoding = fileEncoding,
row.names=row.names, na=na)
write.csv2(x=datosac, file = file7, fileEncoding = fileEncoding,
row.names=row.names, na=na)
}
}



f<-cell/60
f<-round(f,digits=18)
ff<-180/f
cc<-ff*2
matriz<-matrix(-9999, nrow=ff, ncol=cc+1)
col<-c(0,seq(from=-180+f, to=180, by=f))
row<-c(seq(from=90-f, to=-90, by=-f))
matriz<-as.data.frame(matriz)
matriz[,1]<-row
matriz<-rbind(col,matriz)
names(matriz)<-NULL
a<-dim(matriz)


options(warn=-1)
matriz1<-matriz
matriz2<-matriz
matriz3<-matriz
matriz4<-matriz
matriz5<-matriz
matriz6<-matriz

if(matrix==TRUE){
maxLat1<-ceiling(max(x[,3]))
minLat1<-floor(min(x[,3]))
maxLon1<-ceiling(max(x[,2]))
minLon1<-floor(min(x[,2]))
}
else
{
maxLat1<-ceiling(max(x[,2]))
minLat1<-floor(min(x[,2]))
maxLon1<-ceiling(max(x[,1]))
minLon1<-floor(min(x[,1]))
}

rm(x)

r1<-which(abs(row-maxLat1)==min(abs(row-maxLat1)))
r2<-which(abs(row-minLat1)==min(abs(row-minLat1)))
c1<-which(abs(col-minLon1)==min(abs(col-minLon1)))
c2<-which(abs(col-maxLon1)==min(abs(col-maxLon1)))
if(length(c1)==2) c1<-c1[2] else c1<-c1
if(length(c2)==2) c2<-c2[2] else c2<-c2
if(length(r1)==2) c1<-r1[2] else r1<-r1
if(length(r2)==2) r2<-r2[2] else r2<-r2

if(method=="abundance") datosf<-datosf else jkjllk<-0


if(matrix==TRUE){
if(method=="incidence") datosf<-datosac else jkjllk<-0
if(method=="accumulation") datosf<-datosac else jkjllk<-0
}
else{
if(method=="incidence") datosf<-datosf else jkjllk<-0
if(method=="accumulation"){
datosa<-replace(datosf[,-c(1,2)], datosf[,-c(1,2)]>1,1)
datosf<-cbind(datosf[,c(1,2)],datosa)
}
else{
jkjllk<-0
}
}



ZZ<-matrix(c("","","",""), nrow=2)

begin.time<-Sys.time() 
leng1<-length(seq(r1-f, r2+f, by = f))
uio<-1
for (z in seq(r1-f, r2+f, by = f)){

if(uio<=1){
uio<-uio+1
}
else{
if(uio>1){
end.time<-Sys.time() 
end.times <- format(end.time, "%b %d, %Y at %X")
run.time<-difftime(end.time,begin.time,units="secs")
run<-as.numeric(run.time)
run<-run/length(seq(r1-f,z, by=f))
run1<-run*(leng1-length(seq(r1-f,z, by=f)))
if(run1>=3600){
ZZ[2,2]<-"remaining hours...."
}
else{
if(run1<=60) ZZ[2,2]<-"remaining seconds...." else ZZ[2,2]<-"remaining minutes...."
}
if(run1>=3600){
minutes<-run1/3600
}
else{
if(run1<=60) minutes<-run1 else minutes<-run1/60 
}
minutes<-round(minutes, digits=1)
ZZ[1,1]<-end.times
ZZ[2,1]<-minutes
write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)
}
else{
}
}


for (h in seq(c1-f, c2+f, by = f)){

if(matrix==TRUE){
if(elements==0){

if(method=="abundance"){
datosx<-dt1[(dt1[,2]>=col[h-f])&(dt1[,2]<col[h])&(dt1[,3]>=row[z+f])&(dt1[,3]<row[z]),]
}
else{
datosx<-temp[(temp[,2]>=col[h-f])&(temp[,2]<col[h])&(temp[,3]>=row[z+f])&(temp[,3]<row[z]),]
}


dimx<-dim(datosx)

datosL<-datosx[,2:3]


if(dimx[1]==0) datosx<-datosx else datosx<-with(datosx, table(datosx[,5],datosx[,1]))


if(dimx[1]==0) datosx<-datosx else datosx<-datosx[, apply(datosx, 2, sum)!=0]#delete all columns with all values equal to zero 



}
else{
datosx<-datosf[(datosf[,1]>=col[h-f])&(datosf[,1]<col[h])&(datosf[,2]>=row[z+f])&(datosf[,2]<row[z]),]
datosx<-datosx[, apply(datosx, 2, sum)!=0]#delete all columns with all values equal to zero 
datosL<-datosx
}
}

else{
datosx<-datosf[(datosf[,1]>=col[h-f])&(datosf[,1]<col[h])&(datosf[,2]>=row[z+f])&(datosf[,2]<row[z]),]
datosx<-datosx[, apply(datosx, 2, sum)!=0]#delete all columns with all values equal to zero 
datosL<-datosx
}


if(elements==0) datosx<-datosx else datosx<-datosx[,-c(1,2)]


dimy<-dim(datosx)


dimy[is.null(dimy)] <- 0

if(dimy[1]==1){
cut<-sum(datosx, na.rm=TRUE)/dimy[2]
}
else{
cut<-dimy[1]/dimy[2]
}


if(method=="abundance") cut<-5 else cut<-cut



if(method=="abundance" & matrix==TRUE & dimy==0 & sum(as.numeric(datosx), na.rm=TRUE)>0
& length(as.numeric(datosx))>1){
salA<-vegan::estimateR(as.numeric(datosx))

if(is.na(salA[2])==FALSE){
if(salA[2]=="logical (0)") chao1<-NA else chao1<-salA[2]
}

if(is.na(salA[4])==FALSE){
if(salA[4]=="logical (0)") ACE<-NA else ACE<-salA[4]
}

Methods<-c(chao1,ACE)

if(is.na(salA[3])==FALSE){
if(salA[3]=="logical (0)") chao1.se<-NA else chao1.se<-salA[3]
}
else{
chao1.se<-NA
}

if(is.na(salA[5])==FALSE){
if(salA[5]=="logical (0)") ACE.se<-NA else ACE.se<-salA[5]
}
else{
ACE.se<-NA
}

seMethods<-c(chao1.se,ACE.se)

Methods[Methods < 0] <- NA
matriz1[z+2,h]<-dimy[2]
matriz4[z+2,h]<-dimy[1]
Lo1<-subset(datosL[,1], !duplicated(datosL[,1]))
La1<-subset(datosL[,2], !duplicated(datosL[,2]))
Longitude<-mean(Lo1)
Latitude<-mean(La1)

if (estimator==0){
pred<-mean(Methods, na.rm=T)
}
else{
pred<-mean(Methods[estimator], na.rm=T)
}
com<-(salA[1]*100/pred)

temp4<-c(Longitude, Latitude, 1,salA[1], chao1, ACE, com)

if(Longitude!=0 & Latitude!=0){
values2<-rbind(values2,temp4)
if(pred=="NaN") matriz3[z+2,h]<-(-9999) else matriz3[z+2,h]<-(dimy[2]-pred)
if(pred=="NaN"){
matriz2[z+2,h]<-(-9999)
matriz3[z+2,h]<-(-9999)
}
else{
matriz2[z+2,h]<-pred
matriz3[z+2,h]<-(dimy[2]-pred)
if(com<cutoffCompleteness){
matriz2[z+2,h]<-(-9999)
matriz3[z+2,h]<-(-9999)
}
}
}

}




if(dimy[1]==0){
hgh<-1
}
else{
if(dimy==0) {
matriz1[z+2,h]<-dimy[2]
matriz4[z+2,h]<-dimy[1]
Lo1<-subset(datosL[,1], !duplicated(datosL[,1]))
La1<-subset(datosL[,2], !duplicated(datosL[,2]))
Longitude<-mean(Lo1)
Latitude<-mean(La1)
}
else{
if(dimy[2]==1){
matriz1[z+2,h]<-dimy[2]
matriz4[z+2,h]<-dimy[1]
Lo1<-subset(datosL[,1], !duplicated(datosL[,1]))
La1<-subset(datosL[,2], !duplicated(datosL[,2]))
Longitude<-mean(Lo1)
Latitude<-mean(La1)
}
else{

if(dimy[1]==0){
hgh<-1
}
else{
if(dimy[1]==1){
if(method=="abundance"){
XX<-t(as.data.frame(colSums(datosx)))
salA<-vegan::estimateR(XX)

if(is.na(salA[2,1])==FALSE){
if(salA[2, 1]=="logical (0)") chao1<-NA else chao1<-salA[2, 1]
}

if(is.na(salA[4,1])==FALSE){
if(salA[4, 1]=="logical (0)") ACE<-NA else ACE<-salA[4, 1]
}

Methods<-c(chao1,ACE)

if(is.na(salA[3,1])==FALSE){
if(salA[3, 1]=="logical (0)") chao1.se<-NA else chao1.se<-salA[3, 1]
}
else{
chao1.se<-NA
}

if(is.na(salA[5,1])==FALSE){
if(salA[5, 1]=="logical (0)") ACE.se<-NA else ACE.se<-salA[5, 1]
}
else{
ACE.se<-NA
}

seMethods<-c(chao1.se,ACE.se)

Methods[Methods < 0] <- NA
matriz1[z+2,h]<-dimy[2]
matriz4[z+2,h]<-dimy[1]
Lo1<-subset(datosL[,1], !duplicated(datosL[,1]))
La1<-subset(datosL[,2], !duplicated(datosL[,2]))
Longitude<-mean(Lo1)
Latitude<-mean(La1)

if (estimator==0){
pred<-mean(Methods, na.rm=T)
}
else{
pred<-mean(Methods[estimator], na.rm=T)
}
com<-(dimy[2]*100/pred)

temp4<-c(Longitude, Latitude, dimy[1],dimy[2], chao1, ACE, com)
values2<-rbind(values2,temp4)
if(pred=="NaN") matriz3[z+2,h]<-(-9999) else matriz3[z+2,h]<-(dimy[2]-pred)
if(pred=="NaN"){
matriz2[z+2,h]<-(-9999)
matriz3[z+2,h]<-(-9999)
}
else{
matriz2[z+2,h]<-pred
matriz3[z+2,h]<-(dimy[2]-pred)
if(com<cutoffCompleteness){
matriz2[z+2,h]<-(-9999)
matriz3[z+2,h]<-(-9999)
}
}


if(pred=="NaN") matriz5[z+2,h]<-(-9999) else matriz5[z+2,h]<-com
}
else{
fgg<-1
}
matriz1[z+2,h]<-dimy[2]
matriz4[z+2,h]<-dimy[1]
Lo1<-subset(datosL[,1], !duplicated(datosL[,1]))
La1<-subset(datosL[,2], !duplicated(datosL[,2]))
Longitude<-mean(Lo1)
Latitude<-mean(La1)
}
else{


if(method=="incidence"){
if(cut<cutoff){
ICE<-NA
sal$chao<-NA
sal$jack1<-NA
sal$jack2<-NA
sal$boot<-NA
sal$chao.se<-NA
sal$jack1.se<-NA
sal$boot.se<-NA 
Methods<-c(NA,NA,NA,NA,NA)
seMethods<-c(NA,NA,NA)
}
else{
sal<-vegan::specpool(datosx)

ICE<-fossil::ICE(as.matrix(datosx), taxa.row=FALSE)

if(is.na(sal$chao)==FALSE){
if(sal$chao=="logical (0)") sal$chao<-NA else sal$chao<-sal$chao
}

if(is.na(sal$jack1)==FALSE){
if(sal$jack1=="logical (0)") sal$jack1<-NA else sal$jack1<-sal$jack1
}

if(sal$jack2=="logical (0)") sal$jack2<-NA else sal$jack2<-sal$jack2

if(is.na(sal$boot)==FALSE){
if(sal$boot=="logical (0)") sal$boot<-NA else sal$boot<-sal$boot
}

if(is.na(sal$chao.se)==FALSE){
if(sal$chao.se=="logical (0)") sal$chao.se<-NA else sal$chao.se<-sal$chao.se
}

if(is.na(sal$jack1.se)==FALSE){
if(sal$jack1.se=="logical (0)") sal$jack1.se<-NA else sal$jack1.se<-sal$jack1.se
}

if(is.na(sal$boot.se)==FALSE){
if(sal$boot.se=="logical (0)") sal$boot.se<-NA else sal$boot.se<-sal$boot.se
}


Methods<-c(sal$chao,ICE[1],sal$jack1,sal$jack2,sal$boot)

seMethods<-c(sal$chao.se,sal$jack1.se,sal$boot.se)
}
}
else{
fgg<-1
}


if(method=="abundance"){

if(cut<cutoff){
Methods<-c(NA,NA)
}
else{
XX<-t(as.data.frame(colSums(datosx)))
salA<-vegan::estimateR(XX)

if(salA[2, 1]=="logical (0)") chao1<-NA else chao1<-salA[2, 1]
if(salA[4, 1]=="logical (0)") ACE<-NA else ACE<-salA[4, 1]
Methods<-c(chao1,ACE)

if(is.na(salA[3,1])==FALSE){
if(salA[3, 1]=="logical (0)") chao1.se<-NA else chao1.se<-salA[3, 1]
}

if(is.na(salA[5,1])==FALSE){
if(salA[5, 1]=="logical (0)") ACE.se<-NA else ACE.se<-salA[5, 1]
}

seMethods<-c(chao1.se, ACE.se)
}
}
else{
fgg<-1
}

if(method=="accumulation"){
if(cut<cutoff){
Methods<-c(NA,NA,NA,NA)
Methodssp<-c(NA,NA,NA,NA)
sp2<-NA;sp3<-NA;sp4<-NA;sp5<-NA
seMethods<-c(NA,NA,NA,NA)
serandom<-NA;seexact<-NA;secoleman<-NA;serarefaction<-NA
}
else{



cu<- vegan::specaccum(datosx, method="random", permutations = 200)
datosc<-data.frame(cu$richness, cu$sites)
if(curve=="Clench") modelo<-try(nls(cu.richness ~ A*cu.sites/(1+B*cu.sites), data=datosc, trace=T), silent=T) else modelo<-try(nls(cu.richness ~ (A/B)*(1-exp((-B*cu.sites))), start=list(A=0.1, B=0.01),data=datosc, trace=T), silent=T)
res<-summary(modelo)
if(res[1]=="1"){
cu2<-NA
serandom<-NA
}
else{
cu2<-res$parameters[1,1]/res$parameters[2,1]
if(cu2<0){
cu2<-NA
sp2<-NA
serandom<-NA
}
else{
cu2<-cu2
leng<-length(cu$sites)
sp2<-cu$richness[leng]-cu$richness[leng-1]
serandom<-sd(c((res$parameters[1,1]/res$parameters[2,1]),
(res$parameters[1,1]+res$parameters[2,1])/(res$parameters[2,1]+res$parameters[2,2]),
(res$parameters[1,1]-res$parameters[2,1])/(res$parameters[2,1]-res$parameters[2,2])))/sqrt(3)
}
}


cu<- vegan::specaccum(datosx, method="exact")
datosc<-data.frame(cu$richness, cu$sites)
if(curve=="Clench") modelo<-try(nls(cu.richness ~ A*cu.sites/(1+B*cu.sites), data=datosc, trace=T), silent=T) else modelo<-try(nls(cu.richness ~ (A/B)*(1-exp((-B*cu.sites))), start=list(A=0.1, B=0.01),data=datosc, trace=T), silent=T)
res<-summary(modelo)
if(res[1]=="1"){
cu3<-NA
seexact<-NA
}
else{
cu3<-res$parameters[1,1]/res$parameters[2,1]
if(cu3<0){
cu3<-NA
sp3<-NA
seexact<-NA
}
else{
cu3<-cu3
leng<-length(cu$sites)
sp3<-cu$richness[leng]-cu$richness[leng-1]
seexact<-sd(c((res$parameters[1,1]/res$parameters[2,1]),
(res$parameters[1,1]+res$parameters[2,1])/(res$parameters[2,1]+res$parameters[2,2]),
(res$parameters[1,1]-res$parameters[2,1])/(res$parameters[2,1]-res$parameters[2,2])))/sqrt(3)
}
}

cu<- vegan::specaccum(datosx, method="coleman")
datosc<-data.frame(cu$richness, cu$sites)
if(curve=="Clench") modelo<-try(nls(cu.richness ~ A*cu.sites/(1+B*cu.sites), data=datosc, trace=T), silent=T) else modelo<-try(nls(cu.richness ~ (A/B)*(1-exp((-B*cu.sites))), start=list(A=0.1, B=0.01),data=datosc, trace=T), silent=T)
res<-summary(modelo)
if(res[1]=="1"){
cu4<-NA
secoleman<-NA
}
else{
cu4<-res$parameters[1,1]/res$parameters[2,1]
if(cu4<0){
cu4<-NA
sp4<-NA
secoleman<-NA
}
else{
cu4<-cu4
leng<-length(cu$sites)
sp4<-cu$richness[leng]-cu$richness[leng-1]
secoleman<-sd(c((res$parameters[1,1]/res$parameters[2,1]),
(res$parameters[1,1]+res$parameters[2,1])/(res$parameters[2,1]+res$parameters[2,2]),
(res$parameters[1,1]-res$parameters[2,1])/(res$parameters[2,1]-res$parameters[2,2])))/sqrt(3)
}
}

cu<- vegan::specaccum(datosx, method="rarefaction")
datosc<-data.frame(cu$richness, cu$sites)
if(curve=="Clench") modelo<-try(nls(cu.richness ~ A*cu.sites/(1+B*cu.sites), data=datosc, trace=T), silent=T) else modelo<-try(nls(cu.richness ~ (A/B)*(1-exp((-B*cu.sites))), start=list(A=0.1, B=0.01),data=datosc, trace=T), silent=T)
res<-summary(modelo)
if(res[1]=="1"){
cu5<-NA
serarefaction<-NA
}
else{
cu5<-res$parameters[1,1]/res$parameters[2,1]
if(cu5<0){
cu5<-NA
sp5<-NA
serarefaction<-NA
}
else{
cu5<-cu5
leng<-length(cu$sites)
sp5<-cu$richness[leng]-cu$richness[leng-1]
serarefaction<-sd(c((res$parameters[1,1]/res$parameters[2,1]),
(res$parameters[1,1]+res$parameters[2,1])/(res$parameters[2,1]+res$parameters[2,2]),
(res$parameters[1,1]-res$parameters[2,1])/(res$parameters[2,1]-res$parameters[2,2])))/sqrt(3)
}
}

Methods<-c(cu3,cu2,cu4,cu5)
if(cut<cutoff) seMethods<-c(NA,NA,NA,NA) else seMethods<-c(serandom, seexact, secoleman, serarefaction)
Methods[Methods < 0] <- NA
seMethods[seMethods < 0] <- NA
if(cut<cutoff) Methodssp<-c(NA,NA,NA,NA) else Methodssp<-c(sp3, sp2, sp4,sp5)
Methodssp[Methodssp < 0] <- NA
}
}
else{
fgg<-1
}


Methods[Methods < 0] <- NA

if(method=="abundance"){
if(cut<cutoff){
chao1<-NA
ACE<-NA
chao1.se<-NA
ACE.se<-NA
}
else{
chao1<-salA[2, 1]
ACE<-salA[4, 1]
chao1.se<-salA[3, 1]
ACE.se<-salA[5, 1]
}
}
else{
chao1<-NA
ACE<-NA
}


if(estimator==0){
pred<-mean(Methods, na.rm=T)
}
else{
pred<-mean(Methods[estimator], na.rm=T)
}


if(method=="accumulation"){
if (estimator==0){
slope<-mean(Methodssp, na.rm=T)
pred<-mean(Methods, na.rm=T)
}
else{
slope<-mean(Methodssp[estimator], na.rm=T)
pred<-mean(Methods[estimator], na.rm=T)
}
}
else{
}

com<-(dimy[2]*100/pred)

Lo1<-subset(datosL[,1], !duplicated(datosL[,1]))
La1<-subset(datosL[,2], !duplicated(datosL[,2]))
Longitude<-mean(Lo1)
Latitude<-mean(La1)

if(method=="incidence"){
temp4<-c(Longitude, Latitude, dimy[1],dimy[2], sal$chao,ICE[1], sal$jack1,sal$jack2,sal$boot,com)
values1<-rbind(values1,temp4)
setemp4<-c(Longitude, Latitude, dimy[1],dimy[2], sal$chao.se, sal$jack1.se,sal$boot.se)
sevalues1<-rbind(sevalues1,setemp4)
}
else{
wuo<-1
}

if(method=="abundance"){
temp4<-c(Longitude, Latitude, dimy[1],dimy[2], chao1, ACE, com)
values2<-rbind(values2,temp4)

setemp4<-c(Longitude, Latitude, dimy[1],dimy[2], chao1.se, ACE.se)
sevalues2<-rbind(sevalues2,setemp4)
}
else{
wuo<-1
}

if(method=="accumulation"){
temp4<-c(Longitude, Latitude, dimy[1],dimy[2],cu2,cu3,cu4,cu5,sp2,sp3,sp4,sp5,com)
cu2<-NA;cu3<-NA;cu4<-NA;cu5<-NA
sp2<-NA;sp3<-NA;sp4<-NA;sp5<-NA
values3<-rbind(values3,temp4)

setemp4<-c(Longitude, Latitude, dimy[1],dimy[2],serandom,seexact,secoleman,serarefaction)
serandom<-NA;seexact<-NA;secoleman<-NA;serarefaction<-NA
sevalues3<-rbind(sevalues3,setemp4)
}
else{
wuo<-1
}

matriz1[z+2,h]<-dimy[2]
matriz4[z+2,h]<-dimy[1]
if(pred=="NaN") matriz3[z+2,h]<-(-9999) else matriz3[z+2,h]<-(dimy[2]-pred)
if(pred=="NaN"){
matriz2[z+2,h]<-(-9999)
matriz3[z+2,h]<-(-9999)
}
else{
matriz2[z+2,h]<-pred
matriz3[z+2,h]<-(dimy[2]-pred)
if(com<cutoffCompleteness){
matriz2[z+2,h]<-(-9999)
matriz3[z+2,h]<-(-9999)
}
}

if(pred=="NaN") matriz5[z+2,h]<-(-9999) else matriz5[z+2,h]<-com

if(method=="accumulation"){
if(slope=="NaN") matriz6[z+2,h]<-(-9999) else matriz6[z+2,h]<-slope

if(pred=="NaN"){
matriz2[z+2,h]<-(-9999)
}
else{
matriz2[z+2,h]<-pred
if(slope>cutoffSlope){
matriz2[z+2,h]<-(-9999)
}
}
}
else{
matriz6[z+2,h]<-(-9999)
}

}
}
}
}
}


}
}

if(method=="incidence"){
values<-values1[-1,]
sevalues<-sevalues1[-1,]
colnames(values)<-c("Longitude", "Latitude", "Records", "Actual", "Chao", "ICE", "jackknife1",
"jackknife2", "Bootstrap","Completeness")
colnames(sevalues)<-c("Longitude", "Latitude", "Records", "Actual", "Chao.se", "jackknife1.se", "Bootstrap.se")
}
else{
}

if(method=="abundance"){
values<-values2[-1,]
sevalues<-sevalues2[-1,]
colnames(values)<-c("Longitude", "Latitude","Records","Actual","Chao (unbiased variant)", "ACE","Completeness")
colnames(sevalues)<-c("Longitude", "Latitude","Records","Actual","Chao.se", "ACE.se")
}
else{
}

if(method=="accumulation"){
values<-values3[-1,]
sevalues<-sevalues3[-1,]
colnames(values)<-c("Longitude", "Latitude", "Records","Actual","random", "exact", "coleman", "rarefaction", "Slope random", "Slope exact", "Slope coleman", "Slope rarefaction", "Completeness")
colnames(sevalues)<-c("Longitude", "Latitude", "Records","Actual","serandom", "seexact", "secoleman", "serarefaction")

}
else{
}

if (missing(jpg4)) jpg4="Correlation matrix.jpg" else jpg4=jpg4
if (missing(pch)) pch=15 else pch=pch
if (missing(cex)) cex=1.5 else cex=cex
if (missing(cex.labels)) cex.labels=1.5 else cex.labels=cex.labels
if (missing(pchcol)) pchcol="red" else pchcol=pchcol

if(cormatrix==TRUE){
if(jpg==TRUE) jpeg(filename = jpg4, width = 8000, height = 4000, units = "px", pointsize = 14, bg = "white", res = 600) else hhjhk<-1
values10<-values[,-c(1,2,3)]
pairs(values10,main="",cex.main=cex.main, cex=cex, pch=pch, cex.labels = cex.labels, col=pchcol)

if(jpg==TRUE) dev.off() else hhjhk<-1
}

ZZ[1,1]<-end.times
ZZ[2,1]<-"Saving files...."
ZZ[2,2]<-""
write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)


if(dec=="."){
write.csv(x=matriz1,file = file1, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv(x=matriz2,file = file2, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv(x=matriz3,file = file3, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv(x=matriz4,file = file8, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv(x=matriz5,file = file9, fileEncoding = fileEncoding,
row.names=row.names,na=na)
if(method=="accumulation") write.csv(x=matriz6,file = file10, fileEncoding = fileEncoding, row.names=row.names,na=na) else fdssdf<-0
write.csv(x=species, file = file4, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv(x=values, file = file6, fileEncoding = fileEncoding,
row.names=row.names,na="")
write.csv(x=sevalues, file = file11, fileEncoding = fileEncoding,row.names=row.names,na="")
}
else{
write.csv2(x=matriz1,file = file1, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv2(x=matriz2,file = file2, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv2(x=matriz3,file = file3, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv2(x=matriz4,file = file8, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv2(x=matriz5,file = file9, fileEncoding = fileEncoding,
row.names=row.names,na=na)
if(method=="accumulation") write.csv2(x=matriz6,file = file10, fileEncoding = fileEncoding, row.names=row.names,na=na) else fds<-0
write.csv2(x=species, file = file4, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv2(x=values, file = file6, fileEncoding = fileEncoding,row.names=row.names,na="")
write.csv2(x=sevalues, file = file11, fileEncoding = fileEncoding,row.names=row.names,na="")
}


if(method=="accumulation") rty<-6 else rty<-5

d<-length(Area)
AA<-Area[1]
if (AA=="World"){
datos1<-adworld[2:5,]
}
else{
datos1<-rbind(adworld1,adworld2)
}

datos1<-na.exclude(datos1)



for (qw in 1:rty){

if(qw==1){
if(dec=="."){
varscale<-read.csv(file8, header=F)
}
else{
varscale<-read.csv2(file8, header=F)
}

}
else{
}

if(qw==2){
if(dec=="."){
varscale<-read.csv(file1, header=F)
}
else{
varscale<-read.csv2(file1, header=F)
}

matrizz1<-varscale
if(!is.null(end)){
datos1<-replace(matrizz1, matrizz1>=end, end)
datos1[1,]<-varscale[1,]
datos1[,1]<-varscale[,1]
varscale<-datos1
rm(datos1)
rm(matrizz1)
codlegend<-paste(">",end)
}

}
else{
}

if(qw==3){
if(dec=="."){
varscale<-read.csv(file2, header=F)
}
else{
varscale<-read.csv2(file2, header=F)
}

matrizz1<-varscale
if(!is.null(end)){
datos1<-replace(matrizz1, matrizz1>=end, end)
datos1[1,]<-varscale[1,]
datos1[,1]<-varscale[,1]
varscale<-datos1
rm(datos1)
rm(matrizz1)
codlegend<-paste(">",end)
}

}
else{
}

if(qw==4){
if(dec=="."){
varscale<-read.csv(file3, header=F)
}
else{
varscale<-read.csv2(file3, header=F)
}

}
else{
}

if(qw==5){
if(dec=="."){
varscale<-read.csv(file9, header=F)
}
else{
varscale<-read.csv2(file9, header=F)
}

}
else{
}

if(qw==6){
if(dec=="."){
varscale<-read.csv(file10, header=F)
}
else{
varscale<-read.csv2(file10, header=F)
}

}
else{
}

ZZ[1,1]<-end.times
ZZ[2,1]<-"Printing plot"
if(qw==1) ZZ[2,2]<-main4 else hjjuy<-1
if(qw==2) ZZ[2,2]<-main1 else hjjuy<-1
if(qw==3) ZZ[2,2]<-main2 else hjjuy<-1
if(qw==4) ZZ[2,2]<-main3 else hjjuy<-1
if(qw==5) ZZ[2,2]<-main5 else hjjuy<-1
if(qw==6) ZZ[2,2]<-main6 else hjjuy<-1

write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)


if (missing(inc)) inc=0.005 else inc=inc
if (AA=="World"){
if (missing(minLat)) minLat<--90 else minLat<-minLat
if (missing(maxLat)) maxLat<-90 else maxLat<-maxLat
if (missing(minLon)) minLon<--180 else minLon<-minLon
if (missing(maxLon)) maxLon<-180 else maxLon<-maxLon
if(missing(extent)) extent<-TRUE else extent<-extent
if(extent==TRUE) minLat<-minLat1-f else minLat<-minLat
if(extent==TRUE) maxLat<-maxLat1+f else maxLat<-maxLat
if(extent==TRUE) minLon<-minLon1-f else minLon<-minLon
if(extent==TRUE) maxLon<-maxLon1+f else maxLon<-maxLon
}
else{
if (missing(maxLon)){
if(max(datos1$Lon)<0) maxLon<-(max(datos1$Lon)-max(datos1$Lon)*inc) else maxLon<-(max(datos1$Lon)+max(datos1$Lon)*inc)
}
else {
maxLon<-maxLon
}
if (missing(minLon)){
if(min(datos1$Lon)<0) minLon<-(min(datos1$Lon)+min(datos1$Lon)*inc) else minLon<-(min(datos1$Lon)-min(datos1$Lon)*inc)
}
else {
minLon<-minLon
}

if (missing(maxLat)){
if(max(datos1$Lat)<0) maxLat<-(max(datos1$Lat)-max(datos1$Lat)*inc) else maxLat<-(max(datos1$Lat)+max(datos1$Lat)*inc)
}
else {
maxLat<-maxLat
}
if (missing(minLat)){
if(min(datos1$Lat)<0) minLat<-(min(datos1$Lat)+min(datos1$Lat)*inc) else minLat<-(min(datos1$Lat)-min(datos1$Lat)*inc)
}
else {
minLat<-minLat
}
}


if (missing(trans)) trans=c(1,1) else trans<-trans
if (missing(log)) log=c(0,0) else log<-log

Lon<-as.numeric(varscale[1,-1])
varLon<-as.numeric(varscale[1,-1])
a<-length(Lon)
a
for (i in 1:a){
if(i==1) varLon[i]<-((-180+Lon[i])/2) else varLon[i]<-((Lon[i-1]+Lon[i])/2)
}

Lat<-as.numeric(varscale[-1,1])
varLat<-as.numeric(varscale[-1,1])
a<-length(Lat)
a
for (i in 1:a){
if(i==1) varLat[i]<-((90+Lat[i])/2) else varLat[i]<-((Lat[i-1]+Lat[i])/2)
}
varLat<-(-varLat)


firstrow<-varscale[1,]

ajuste<-varscale[varscale[,1]<=maxLat&varscale[,1]>=minLat,]

if(firstrow==ajuste[1,]){
}
else{
ajuste<-rbind(firstrow,ajuste)
}

ajuste<-ajuste[,ajuste[1,]<=maxLon&ajuste[1,]>=minLon]



ajuste<-ajuste[-1,-1]

ajuste<-as.matrix(ajuste)



if(trans[1]==0){
ajuste<-replace(ajuste, ajuste==-9999,NA)
ajuste<-ajuste/trans[2]
ajuste<-replace(ajuste, is.na(ajuste),-9999)
}
else{
ajuste<-replace(ajuste, ajuste==-9999,NA)
ajuste<-ajuste*trans[2]
ajuste<-replace(ajuste, is.na(ajuste),-9999)
}

if(log[1]==0){
ajuste<-ajuste
}
else{
ajuste<-replace(ajuste, ajuste==-9999,NA)
ajuste<-log(ajuste+log[2])
ajuste<-replace(ajuste, is.na(ajuste),-9999)
}


ajuste<- ajuste[nrow(ajuste):1,]

varscale<-varscale[-1,-1]

varscale<-as.matrix(varscale)

if(trans[1]==0){
varscale<-replace(varscale, varscale==-9999,NA)
varscale<-varscale/trans[2]
varscale<-replace(varscale, is.na(varscale),-9999)
}
else{
varscale<-replace(varscale, varscale==-9999,NA)
varscale<-varscale*trans[2]
varscale<-replace(varscale, is.na(varscale),-9999)
}

if(log[1]==0){
varscale<-varscale
}
else{
varscale<-replace(varscale, varscale==-9999,NA)
varscale<-log(varscale+log[2])
varscale<-replace(varscale, is.na(varscale),-9999)
}

varscale<- varscale[nrow(varscale):1,]
varscale<-t(varscale)

if (maxLon>=180) maxLon<-180 else maxLon<-maxLon
if (minLon<=-180) minLon<--180 else minLon<-minLon
if (maxLat>=90) maxLat<-90 else maxLat<-maxLat
if (minLat<=-90) minLat<--90 else minLat<-minLat

if (missing(Area)) Area="World" else Area=Area
if (missing(colbg)) colbg="transparent" else colbg=colbg
if (missing(colcon)) colcon="transparent" else colcon=colcon
if (missing(colf)) colf="black" else colf=colf
if (missing(colfexc)) colfexc="black" else colfexc=colfexc
if (missing(pro)) pro=TRUE else pro=pro
if (missing(ylab)) ylab="Latitude" else ylab=ylab
if (missing(xlab)) xlab="Longitude" else xlab=xlab
if (missing(main1)) main1="Actual richness" else main1=main1
if (missing(main2)) main2="Predicted richness" else main2=main2
if (missing(main3)) main3="Residuals" else main3=main3
if (missing(main4)) main4="Records" else main4=main4
if (missing(main5)) main5="Completeness" else main5=main5
if (missing(main6)) main6="Slope" else main6=main6

if (missing(cex.lab)) cex.lab=1.4 else cex.lab=cex.lab
if (missing(cex.axis)) cex.axis=1.2 else cex.axis=cex.axis
if (missing(lwdP)) lwdP=0.6 else lwdP=lwdP
if (missing(lwdC)) lwdC=0.1 else lwdC=lwdC
if (missing(family)) family="sans" else family=family
if (missing(font.main)) font.main=2 else font.main=font.main
if (missing(font.lab)) font.lab=1 else font.lab=font.lab
if (missing(font.axis)) font.axis=1 else font.axis=font.axis
if (missing(lab)) lab=NULL else lab=lab
if (missing(exclude)) exclude=NULL else exclude=exclude
if (missing(colexc)) colexc="white" else colexc=colexc
if (missing(varscale)) varscale=NULL else varscale=varscale
color<-c("#FFFFFFFF", "#C8FFFFFF","#64FFFFFF","#00FFFFFF", "#64FF64FF","#C8FF00FF", "#FFFF00FF","#FFC800FF","#FF6400FF","#FF0000FF")
if (missing(colscale)) colscale=color else colscale=colscale

if (missing(breaks)) breaks=10 else breaks=breaks
if (missing(xl)) xl=0 else xl=xl
if (missing(xr)) xr=0 else xr=xr
if (missing(yb)) yb=0 else yb=yb
if (missing(yt)) yt=0 else yt=yt
if (missing(ndigits)) ndigits=0 else ndigits=ndigits


legend.max=max(ajuste)

if(legend.max<=10){
legend.min=(if(min(ajuste[!ajuste==-9999])==0) min(ajuste[!ajuste==-9999])+(max(ajuste)/(length(colscale)-1)) else min(ajuste[!ajuste==-9999]))
}
else{
legend.min=min(ajuste[!ajuste==-9999])
}

if(legend.min<0) legend.min<-legend.min+legend.min*0.1/100 else legend.min<-legend.min-legend.min*0.1/100
if(legend.max<0) legend.max<-legend.max-legend.max*0.1/100 else legend.max<-legend.max+legend.max*0.1/100


Lati<-(maxLat+minLat)/2
if (pro==TRUE) aspe=(1/cos(Lati*pi/180)) else aspe=1
if (missing(asp)) asp=aspe else asp=asp

x<-0
y<-0
rm(datos1)


if (missing(jpg1)) jpg1="Actual richness.jpg" else jpg1=jpg1
if (missing(jpg2)) jpg2="Predicted richness.jpg" else jpg2=jpg2
if (missing(jpg3)) jpg3="Residuals.jpg" else jpg3=jpg3
if (missing(jpg5)) jpg5="Records.jpg" else jpg5=jpg5
if (missing(jpg6)) jpg6="Completeness.jpg" else jpg6=jpg6
if (missing(jpg7)) jpg7="Slope.jpg" else jpg7=jpg7

if(qw==1) file=jpg5 else hjjuy<-1
if(qw==2) file=jpg1 else hjjuy<-1
if(qw==3) file=jpg2 else hjjuy<-1
if(qw==4) file=jpg3 else hjjuy<-1
if(qw==5) file=jpg6 else hjjuy<-1
if(qw==6) file=jpg7 else hjjuy<-1

if(jpg==TRUE) jpeg(filename = file, width = 8000, height = 4000, units = "px", pointsize = 14, bg = "white", res = 600) else hhjhk<-1


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

if(qw==1 | qw==4){
if(min(varscale[!varscale==-9999])==0) iniF<-0 else iniF<-legend.min
}

if(qw==2 | qw==3){
if(!is.null(ini)){
iniF<-ini
}
else{
if(min(varscale[!varscale==-9999])==0) iniF<-0 else iniF<-legend.min
}
}


if(qw==5){
if (ComCut==TRUE){
iniF<-cutoffCompleteness
}
else{
if(min(varscale[!varscale==-9999])==0) iniF<-(-0.00001) else iniF<-legend.min
}
}

if(qw==6){
if (SlopeCut==TRUE){
legend.max<-cutoffSlope
if(min(varscale[!varscale==-9999])==0) iniF<-(-0.00001) else iniF<-legend.min
}
else{
if(min(varscale[!varscale==-9999])==0) iniF<-(-0.00001) else iniF<-legend.min
}
}




if (maxLon==180 & minLon==-180 & minLat==-90 & maxLat==90){
xl<-185
xr<-195
}


if(qw==1) main=main4 else hjjuy<-1
if(qw==2) main=main1 else hjjuy<-1
if(qw==3) main=main2 else hjjuy<-1
if(qw==4) main=main3 else hjjuy<-1
if(qw==5) main=main5 else hjjuy<-1
if(qw==6){
main=main6
ndigits=2
color7<-rev(colscale[-1])
colscale<-append(colscale[1],color7)
}


par(lwd=lwdP,fg="black",family=family)

tmp<-squishplot(xlim=c(minLon,maxLon), ylim=c(minLat,maxLat), asp=aspe)

legend.freq1=abs((legend.max-iniF)/(length(colscale)-1))

legend.freq=abs((legend.max-iniF)/(breaks-1))

if (missing(legend.pos)){
if((maxLon-minLon)>260 & (maxLon-minLon)/(maxLat-minLat)>2.265) legend.pos="x" else legend.pos="y"
}

if(legend.pos=="x"){
if (missing(cex.main)) cex.main=1.3 else cex.main=cex.main
}

if(legend.pos=="y"){
if (missing(cex.main)) cex.main=1.6 else cex.main=cex.main
}


if (legend.pos=="y") par(oma=c(0,0,0,1)) else  par(oma=c(0,0,2,0))
image(varLon, varLat,varscale,xlim=c(minLon,maxLon),ylim=c(minLat,maxLat), axes=F, xaxs="i", yaxs="i", xlab="",ylab="", col=colscale, breaks=c(iniF,seq(iniF,legend.max,by=legend.freq1)))

par(new=T,lwd=lwdP)

plot(x,y,xlim=c(minLon,maxLon),ylim=c(minLat,maxLat),xlab=xlab, main="", axes=TRUE,
ylab = ylab, cex.lab=cex.lab, cex.axis= cex.axis,type="n",bty="l",
font.lab=font.lab, font.axis=font.axis,lab=lab,yaxs="i",xaxs="i",yaxt="n",xaxt="n")
mtext(text=main,side=3, line=0.3, cex=cex.main, font=font.main)


axis(side=1,xlim=c(minLon,maxLon),lwd=lwdP)
axis(side=2,ylim=c(minLat,maxLat),lwd=lwdP)

if (colbg=="#FFFFFF") rect(0, 0, 0, 0, col = colbg) else rect(minLon, minLat, maxLon, maxLat, col = colbg)

if (legend.pos=="y"){
if (xl==0){
x1<-(maxLon-minLon)*(-0.00106495)+0.747382095+maxLon
x2<-(maxLon-minLon)*(-0.003194851)+2.060146284+maxLon
}
else{
x1<-xl
x2<-xr
}

if(legend.max<=10){
sequ<-(seq(iniF,legend.max,by=legend.freq))
sequ<-round(sequ, digits=ndigits)

}
else{
if(iniF==0){
legend.freq=abs((legend.max-iniF)/(breaks-1))
sequ<-(seq(iniF,legend.max,by=legend.freq))
sequ<-round(sequ, digits=ndigits)

}
else{
sequ<-(seq(iniF,legend.max,by=legend.freq))
sequ<-round(sequ, digits=ndigits)


}

}

if(qw==2 | qw==3){ 
if(!is.null(end)){
lensequ<-length(sequ)
sequ[lensequ]<-codlegend
}
}


plotrix::color.legend(xl=x1, yb=minLat, xr= x2,
yt=maxLat, sequ, gradient="y", align="rb", cex=1.2, rect.col=colscale[-1])
}
else{
if (yb==0){
if(!is.null(main)){
y1<-maxLat+(maxLat-minLat)*(0.101851852)-1.333333333
y2<-maxLat+(maxLat-minLat)*(0.157407407)-1.333333333
}
else{
y1<-maxLat+(maxLat-minLat)*(0.027777778)
y2<-maxLat+(maxLat-minLat)*(0.083333333)
}
}
else{
y1<-yb
y2<-yt
}

if(legend.max<=10){
sequ<-(seq(iniF,legend.max,by=legend.freq))
sequ<-round(sequ, digits=ndigits)
}
else{
sequ<-(seq(iniF,legend.max,by=legend.freq))
sequ<-round(sequ, digits=ndigits)
}

if(qw==2 | qw==3){ 
if(!is.null(end)){
lensequ<-length(sequ)
sequ[lensequ]<-codlegend
}
}

plotrix::color.legend(xl=minLon, yb=y1, xr=maxLon, yt=y2, sequ,
gradient="x", align="lt", cex=1.2, rect.col=colscale[-1])
}

if (AA=="World") {
polygon(adworld$Lon,adworld$Lat,col=colcon, border=colf)
if(!is.null(exclude)){
polygon(adworld2$Lon,adworld2$Lat,col=colexc, border=colfexc)
}
}
else {
polygon(adworld1$Lon,adworld1$Lat,col=colcon, border=colf)
polygon(adworld2$Lon,adworld2$Lat,col=colexc, border=colfexc)
}

par(tmp)

if(jpg==TRUE) dev.off() else hhjk<-1
}
ZZ[1,1]<-"END"
ZZ[1,2]<-""
ZZ[2,1]<-""
ZZ[2,2]<-""
write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)
rm(datos2)
rm(datosac)
rm(datosf)
rm(datosL)
rm(datosx)
rm(matriz)
rm(matriz1)
rm(matriz2)
rm(matriz3)
rm(matriz4)
rm(matriz5)
rm(matriz6)
rm(values)
rm(values1)
rm(values2)
rm(values3)
rm(varscale)
rm(temp)
rm(temp4)
}
