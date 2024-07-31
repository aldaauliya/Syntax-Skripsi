library(maptools)
library(spdep)
library(sp)
library(spData)
library(INLA)
library(sf)
library(terra)
library(RColorBrewer)
library(gstat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(sarima)
library(extrafont)
library(shapefiles)
library(lattice)
library(splines2)
library(Matrix)
library(CAMAN)
library(CARBayes)
library(CARBayesST)
library(spatialreg)
library(brinla)
library(gplots)
library(tmap)
library(rgeos)
library(rgdal)
library(spData)
library(classInt)
library(raster)
library(ggplot2)
library(ggsn)
library(gtools)
library(lmtest)
library(splm) 
library(plm) 
library(tseries)
library(mvtnorm)
library(orcutt)
library(HoRM)
library(ctv)
library(splines)
library(splancs)


#Input MAP
Indonesia <- getData('GADM', country='IDN', level=3)
Bandung <- Indonesia[Indonesia$NAME_2 == "Kota Bandung",]
plot(Bandung)
Kecamatan <- as.data.frame(Bandung$NAME_3);Kecamatan

#Input Data
setwd("D:/SKRIPSI")
DBD <- read.csv("DATA SKRIPSI FIX.csv", sep=";", dec=",")
head(DBD)
summary(DBD)
str(DBD)

#Analisis Deskriptif
summary_per_tahun <- DBD %>%
  group_by(Tahun) %>%
  summarise(
    Total_Penduduk = sum(Jumlah.Penduduk),
    Rata_Rata_Y = mean(Y),
    Maxmimum_Y = max(Y),
    Manimum_Y = min(Y),
    STDEV_Y = sd(Y)
  )

print(summary_per_tahun)

# - Buat ID untuk menghubungkan Peta dengan Data
ID <- c(1:30)                                                                                                                                                                                                                                                                                                                                                                                                                                              
Bandung$id <- ID
CoordK <- coordinates(Bandung)
XC <- coordinates(Bandung)[,1]
YC <- coordinates(Bandung)[,2]
DBD$id <- rep(ID, 60) 
DBD$XC <- rep(XC, 60)
DBD$YC <- rep(YC, 60)

# - Menghubungkan Peta dengan Data
BandungData <- fortify(Bandung)
id <- as.numeric(unique(BandungData$id))
DBD$id <-rep(id,60)

data.shp <- merge(BandungData, DBD, by="id", all.X=TRUE)
dataMap <- data.shp[order(data.shp$order), ]

## MATRIKS PEMBOBOT SPASIAL ##
# Queen
W <- poly2nb(Bandung, row.names=ID, queen=TRUE) #Mendapatkan W 
WB <- nb2mat(W, style='B', zero.policy = TRUE) #Menyajikan dalam bentuk matrix biner "B" 
Wls<-nb2listw(W, zero.policy = TRUE)
Ws <- as(as_dgRMatrix_listw(Wls), "CsparseMatrix")


#atau ini gatau betul apa tidak
#Ws <- as(spatialreg::as_dgRMatrix_listw(Wls),"CsparseMatrix")

plot(Bandung, col="pink",main="Queen")
plot(Wls,CoordK,add=T, col="red", lwd=1)


## PENGUJIAN AUTOKORELASI SPATIALTEMPORAL ##
Wsmatrix <- as.matrix(Ws)
KasusDBD <- DBD$Y
Y1 <- as.matrix(KasusDBD)
MoranST <- MoranST.MC(Y1,Wsmatrix,100) #function MoranST.MC nya di beda sheet


##PEMODELAN
dbd <- DBD$Y
KP <- DBD$X1
PHBS <- DBD$X2
CHUJAN <- DBD$X3
SUHU <- DBD$X4
LEMBAB <- DBD$X5
Ei <- DBD$Ei

#INLA SETUP
control<- list(
  predictor = list(compute = TRUE,link=1),
  results = list(return.marginals.random = TRUE,
                 return.marginals.predictor=TRUE),
  compute = list(hyperpar=TRUE, return.marginals=TRUE, dic=TRUE, mlik =
                   TRUE, cpo = TRUE,
                 po = TRUE, waic=TRUE, graph=TRUE, openmp.strategy="huge"))

#MODEL YG HANYA MEMPERHATIKAN EFEK FIXED
Model1 <- dbd~KP+PHBS+CHUJAN+SUHU+LEMBAB
RModelRun1<- inla(Model1, family="poisson",data=DBD,E=Ei,
                  control.compute = control$compute,
                  control.predictor = control$predictor,
                  control.inla = list(tolerance=1e-20, h=1e-8,int.strategy = "eb", strategy ="simplified.laplace")
)
summary(RModelRun1)


#MODEL MEMPERHATIKAN EFEK SPASIAL = LEROUX
#PAKE besagproper2 
Model2 <- dbd~f(ID1,model="besagproper2",graph = Ws, constr = T)+KP+PHBS+CHUJAN+SUHU+LEMBAB
RModelRun2 <- inla(Model2,family="poisson",data=DBD,E=Ei,
                   control.compute = control$compute,
                   control.predictor = control$predictor,
                   control.inla = list(tolerance=1e-20, h=1e-8,int.strategy = "eb", strategy = "simplified.laplace")
)
summary(RModelRun2)

#MODEL MEMPERHATIKAN EFEK TEMPORAL = RW2
Model3<- dbd~f(IT1,model="rw2", constr=T)+KP+PHBS+CHUJAN+SUHU+LEMBAB
RModelRun3<- inla(Model3, family="poisson",data=DBD,E=Ei,
                  control.compute = control$compute,
                  control.predictor = control$predictor,
                  control.inla = list(tolerance=1e-20, h=1e-8,int.strategy = "eb", strategy ="simplified.laplace")
)
summary(RModelRun3)

#MODEL MEMPERHATIKAN EFEK SPASIAL DAN TEMPORAL
DBD$IT2 <- DBD$IT1
DBD$IT3<-DBD$IT1
DBD$ID2<-DBD$ID1

Model4 <- dbd ~ KP+PHBS+CHUJAN+SUHU+LEMBAB+
  f(IT1,model="rw2", constr=T)+
  f(ID1,model="besagproper2",graph=Ws, constr=T)

RModelRun4 <- inla(Model4, family="poisson",data=DBD,E=Ei,
                  control.compute = control$compute,
                  control.predictor = control$predictor,
                  control.inla = list(tolerance=1e-20, h=1e-8,int.strategy = "eb", strategy ="simplified.laplace")
)
summary(RModelRun4)


#MODEL INTERAKSI
Model5<- dbd ~ KP+PHBS+CHUJAN+SUHU+LEMBAB+
  f(ID1,model="besagproper2",graph=Ws, constr=T)+ 
  f(IT1,model="rw2", constr=T)+
  f(ID2,model="besagproper2",graph=Ws,initial=1,constr=T, group=IT2,control.group=list(model="rw2"))

RModelRun5<- inla(Model5, family="poisson",data=DBD,E=Ei,
                  control.compute = control$compute,
                  control.predictor = control$predictor,
                  control.inla = list(tolerance=1e-20, h=1e-8,int.strategy = "eb", strategy = "simplified.laplace")
)
summary(RModelRun5)


## - R Korelasi 
pred1=RModelRun1$summary.fitted.values$mean*Ei
cor1=cor(pred1,dbd,use="complete.obs")

pred2=RModelRun2$summary.fitted.values$mean*Ei
cor2=cor(pred2,dbd,use="complete.obs")

pred3=RModelRun3$summary.fitted.values$mean*Ei
cor3=cor(pred3,dbd,use="complete.obs")

pred4=RModelRun4$summary.fitted.values$mean*Ei
cor4=cor(pred4,dbd,use="complete.obs")

pred5=RModelRun5$summary.fitted.values$mean*Ei
cor5=cor(pred5,dbd,use="complete.obs")

DIC=c(RModelRun1$dic$dic,RModelRun2$dic$dic,
      RModelRun3$dic$dic,RModelRun4$dic$dic,
      RModelRun5$dic$dic)
cor=c(cor1,cor2,cor3,cor4,cor5)
data.frame(DIC,cor) 


#EVALUASI MODEL
plot(RModelRun5)
bri.hyperpar.summary(RModelRun5)

## MENGHITUNG TAKSIRAN RISIKO RELATIF ##

#Model fix buat RR nya
Modelfix<- dbd ~ KP+PHBS+SUHU+CHUJAN+LEMBAB+
  f(IT1,model="rw2", constr=T)+
  f(ID2,model="besagproper2",graph=Ws,initial=1,constr=T, group=IT2,control.group=list(model="rw2"))

RModelfix<- inla(Modelfix, family="poisson",data=DBD,E=Ei,
                 control.compute = control$compute,
                 control.predictor = control$predictor,
                 control.inla = list(tolerance=1e-20, h=1e-8,int.strategy = "eb", strategy = "simplified.laplace")
)
summary(RModelfix)

RR <- RModelfix$summary.fitted.values$mean
write.table(RR, "Taksiran_RR_BARU.csv", sep = ";", row.names = FALSE, col.names = TRUE)
head(RR)

#PLOT RR PER TAHUN
# Data yang diberikan
tahun <- c(2019, 2020, 2021, 2022, 2023)
rata_rata <- c(1.320266463, 0.848486289, 1.066935986, 1.521728747, 0.546966254)

# Membuat data frame 
data <- data.frame(Tahun = tahun, Rata_rata = rata_rata)

# Membuat plot garis
plot(data$Tahun, data$Rata_rata,type = "l", col = "black", xlab = "Tahun", ylab = "Rata-rata", main = "Plot Rata-rata per Tahun")

#PEMETAAN
#RR 2019
STB1_19=RModelfix$summary.fitted.values$mean[1:30]
STB2_19=RModelfix$summary.fitted.values$mean[31:60]
STB3_19=RModelfix$summary.fitted.values$mean[61:90]
STB4_19=RModelfix$summary.fitted.values$mean[91:120]
STB5_19=RModelfix$summary.fitted.values$mean[121:150]
STB6_19=RModelfix$summary.fitted.values$mean[151:180]
STB7_19=RModelfix$summary.fitted.values$mean[181:210]
STB8_19=RModelfix$summary.fitted.values$mean[211:240]
STB9_19=RModelfix$summary.fitted.values$mean[241:270]
STB10_19=RModelfix$summary.fitted.values$mean[271:300]
STB11_19=RModelfix$summary.fitted.values$mean[301:330]
STB12_19=RModelfix$summary.fitted.values$mean[331:360]
RR_2019 = data.frame(STB1_19,STB2_19,STB3_19,STB4_19,
STB5_19,STB6_19,STB7_19,STB8_19,STB9_19,STB10_19,
STB11_19,STB12_19)



#RR2020
STB1_20=RModelfix$summary.fitted.values$mean[361:390]
STB2_20=RModelfix$summary.fitted.values$mean[391:420]
STB3_20=RModelfix$summary.fitted.values$mean[421:450]
STB4_20=RModelfix$summary.fitted.values$mean[451:480]
STB5_20=RModelfix$summary.fitted.values$mean[481:510]
STB6_20=RModelfix$summary.fitted.values$mean[511:540]
STB7_20=RModelfix$summary.fitted.values$mean[541:570]
STB8_20=RModelfix$summary.fitted.values$mean[571:600]
STB9_20=RModelfix$summary.fitted.values$mean[601:630]
STB10_20=RModelfix$summary.fitted.values$mean[631:660]
STB11_20=RModelfix$summary.fitted.values$mean[661:690]
STB12_20=RModelfix$summary.fitted.values$mean[691:720]

#RR2021
STB1_21=RModelfix$summary.fitted.values$mean[721:750]
STB2_21=RModelfix$summary.fitted.values$mean[751:780]
STB3_21=RModelfix$summary.fitted.values$mean[781:810]
STB4_21=RModelfix$summary.fitted.values$mean[811:840]
STB5_21=RModelfix$summary.fitted.values$mean[841:870]
STB6_21=RModelfix$summary.fitted.values$mean[871:900]
STB7_21=RModelfix$summary.fitted.values$mean[901:930]
STB8_21=RModelfix$summary.fitted.values$mean[931:960]
STB9_21=RModelfix$summary.fitted.values$mean[961:990]
STB10_21=RModelfix$summary.fitted.values$mean[991:1020]
STB11_21=RModelfix$summary.fitted.values$mean[1021:1050]
STB12_21=RModelfix$summary.fitted.values$mean[1051:1080]


#RR2022
STB1_22=RModelfix$summary.fitted.values$mean[1081:1110]
STB2_22=RModelfix$summary.fitted.values$mean[1111:1140]
STB3_22=RModelfix$summary.fitted.values$mean[1141:1170]
STB4_22=RModelfix$summary.fitted.values$mean[1171:1200]
STB5_22=RModelfix$summary.fitted.values$mean[1201:1230]
STB6_22=RModelfix$summary.fitted.values$mean[1231:1260]
STB7_22=RModelfix$summary.fitted.values$mean[1261:1290]
STB8_22=RModelfix$summary.fitted.values$mean[1291:1320]
STB9_22=RModelfix$summary.fitted.values$mean[1321:1350]
STB10_22=RModelfix$summary.fitted.values$mean[1351:1380]
STB11_22=RModelfix$summary.fitted.values$mean[1381:1410]
STB12_22=RModelfix$summary.fitted.values$mean[1411:1440]


#RR2023
STB1_23=RModelfix$summary.fitted.values$mean[1441:1470]
STB2_23=RModelfix$summary.fitted.values$mean[1471:1500]
STB3_23=RModelfix$summary.fitted.values$mean[1501:1530]
STB4_23=RModelfix$summary.fitted.values$mean[1531:1560]
STB5_23=RModelfix$summary.fitted.values$mean[1561:1590]
STB6_23=RModelfix$summary.fitted.values$mean[1591:1620]
STB7_23=RModelfix$summary.fitted.values$mean[1621:1650]
STB8_23=RModelfix$summary.fitted.values$mean[1651:1680]
STB9_23=RModelfix$summary.fitted.values$mean[1681:1710]
STB10_23=RModelfix$summary.fitted.values$mean[1711:1740]
STB11_23=RModelfix$summary.fitted.values$mean[1741:1770]
STB12_23=RModelfix$summary.fitted.values$mean[1771:1800]

#BIKIN PETA 
#Input MAP
Indonesia <- getData('GADM', country='IDN', level=3)
Bandung <- Indonesia[Indonesia$NAME_2 == "Kota Bandung",]
plot(Bandung)
Kecamatan <- as.data.frame(Bandung$NAME_3)

# - Buat ID untuk menghubungkan Peta dengan Data
ID <- c(1:30)                                                                                                                                                                                                                                                                                                                                                                                                                                              
Bandung$id <- ID
Coordinate <- coordinates(Bandung)
Koordinate <- data.frame(long=Coordinate[,1],
lat=Coordinate[,2],id=ID)
XC <- coordinates(Bandung)[,1]
YC <- coordinates(Bandung)[,2]

# - Menghubungkan Peta dengan Data
BandungData <- fortify(Bandung)
head(BandungData)
BandungData$id <-ID
DBD$id <-rep(id,60)

data.shp <- merge(BandungData, DBD, by="id", all.X=TRUE)
dataMap <- data.shp[order(data.shp$order), ]

#PETA 2019
#Setting data
#2019 Jan
dataku119=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB1_19,
Kecamatan=Bandung$NAME_3)
merge.shp119=merge(BandungData,dataku119,by="id",all.x=T)
FinalData119=merge.shp119[order(merge.shp123$order),]

#2019 Feb
dataku219=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB2_19,
Kecamatan=Bandung$NAME_3)
merge.shp219=merge(BandungData,dataku219,by="id",all.x=T)
FinalData219=merge.shp219[order(merge.shp219$order),]

#2019 Mar
dataku319=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB3_19,
Kecamatan=Bandung$NAME_3)
merge.shp319=merge(BandungData,dataku319,by="id",all.x=T)
FinalData319=merge.shp319[order(merge.shp319$order),]

#2019 April
dataku419=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB4_19,
Kecamatan=Bandung$NAME_3)
merge.shp419=merge(BandungData,dataku419,by="id",all.x=T)
FinalData419=merge.shp419[order(merge.shp419$order),]

#2019 Mei
dataku519=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB5_19,
Kecamatan=Bandung$NAME_3)
merge.shp519=merge(BandungData,dataku519,by="id",all.x=T)
FinalData519=merge.shp519[order(merge.shp519$order),]

#2019 Juni
dataku619=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB6_19,
Kecamatan=Bandung$NAME_3)
merge.shp619=merge(BandungData,dataku619,by="id",all.x=T)
FinalData619=merge.shp619[order(merge.shp619$order),]

#2019 Juli
dataku719=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB7_19,
Kecamatan=Bandung$NAME_3)
merge.shp719=merge(BandungData,dataku719,by="id",all.x=T)
FinalData719=merge.shp719[order(merge.shp719$order),]

#2019 Agustus
dataku819=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB8_19,
Kecamatan=Bandung$NAME_3)
merge.shp819=merge(BandungData,dataku819,by="id",all.x=T)
FinalData819=merge.shp819[order(merge.shp819$order),]

#2019 September
dataku919=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB9_19,
Kecamatan=Bandung$NAME_3)
merge.shp919=merge(BandungData,dataku919,by="id",all.x=T)
FinalData919=merge.shp919[order(merge.shp919$order),]

#2019 Okt
dataku1019=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB10_19,
Kecamatan=Bandung$NAME_3)
merge.shp1019=merge(BandungData,dataku1019,by="id",all.x=T)
FinalData1019=merge.shp1019[order(merge.shp1019$order),]

#2019 Nov
dataku1119=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB11_19,
Kecamatan=Bandung$NAME_3)
merge.shp1119=merge(BandungData,dataku1119,by="id",all.x=T)
FinalData1119=merge.shp1119[order(merge.shp1119$order),]

#2019 Dec
dataku1219=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB12_19,
Kecamatan=Bandung$NAME_3)
merge.shp1219=merge(BandungData,dataku1219,by="id",all.x=T)
FinalData1219=merge.shp1019[order(merge.shp1219$order),]


#Final
bulan=rep(month.name,each=18299)
FinalBDGData19=rbind(FinalData119,FinalData219,FinalData319,
		FinalData419,FinalData519,FinalData619,
		FinalData719,FinalData819,FinalData919,
		FinalData1019,FinalData1119,FinalData1219)
FinalBDGData19$Bulan=bulan
str(FinalBDGData19)


FinalBDGData19$Bulan <- factor(FinalBDGData19$Bulan, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
#Ggplot Map
nlcr=5
plotclr=brewer.pal(nlcr,"Reds")

p= ggplot(data=FinalBDGData19)+
geom_polygon(aes(x=long,y=lat,group=group,fill=ST),color="black",size=0.25)+
geom_text(aes(x=long,y=lat,label = id),data = Koordinate, size = 3, hjust = 0.5)+
facet_wrap(.~Bulan,scales="free",ncol=3)+
scale_fill_gradientn(colours=plotclr)+
theme_bw(base_size=15)+
ylab("")+
xlab("")+
labs(fill = "ST")+
theme(legend.position="right")+
theme(text =element_text(size=12)) +
guides(fill = guide_colorbar(title.position = "top",title.vjust = 1, frame.colour = "black",barwidth = 1,barheight = 10))

#SIMPAN
ggsave("plot2019fix.png", plot = p, width = 12, height = 8, dpi = 300)

#PETA 2020
#Setting data
#2020 Jan
dataku120=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB1_20,
Kecamatan=Bandung$NAME_3)
merge.shp120=merge(BandungData,dataku120,by="id",all.x=T)
FinalData120=merge.shp120[order(merge.shp123$order),]

#2020 Feb
dataku220=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB2_20,
Kecamatan=Bandung$NAME_3)
merge.shp220=merge(BandungData,dataku220,by="id",all.x=T)
FinalData220=merge.shp220[order(merge.shp220$order),]

#2020 Mar
dataku320=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB3_20,
Kecamatan=Bandung$NAME_3)
merge.shp320=merge(BandungData,dataku320,by="id",all.x=T)
FinalData320=merge.shp320[order(merge.shp320$order),]

#2020 April
dataku420=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB4_20,
Kecamatan=Bandung$NAME_3)
merge.shp420=merge(BandungData,dataku420,by="id",all.x=T)
FinalData420=merge.shp420[order(merge.shp420$order),]

#2020 Mei
dataku520=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB5_20,
Kecamatan=Bandung$NAME_3)
merge.shp520=merge(BandungData,dataku520,by="id",all.x=T)
FinalData520=merge.shp520[order(merge.shp520$order),]

#2020 Juni
dataku620=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB6_20,
Kecamatan=Bandung$NAME_3)
merge.shp620=merge(BandungData,dataku620,by="id",all.x=T)
FinalData620=merge.shp620[order(merge.shp620$order),]

#2020 Juli
dataku720=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB7_20,
Kecamatan=Bandung$NAME_3)
merge.shp720=merge(BandungData,dataku720,by="id",all.x=T)
FinalData720=merge.shp720[order(merge.shp720$order),]

#2020 Agustus
dataku820=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB8_20,
Kecamatan=Bandung$NAME_3)
merge.shp820=merge(BandungData,dataku820,by="id",all.x=T)
FinalData820=merge.shp820[order(merge.shp820$order),]

#2020 September
dataku920=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB9_20,
Kecamatan=Bandung$NAME_3)
merge.shp920=merge(BandungData,dataku920,by="id",all.x=T)
FinalData920=merge.shp920[order(merge.shp920$order),]

#2020 Okt
dataku1020=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB10_20,
Kecamatan=Bandung$NAME_3)
merge.shp1020=merge(BandungData,dataku1020,by="id",all.x=T)
FinalData1020=merge.shp1020[order(merge.shp1020$order),]

#2020 Nov
dataku1120=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB11_20,
Kecamatan=Bandung$NAME_3)
merge.shp1120=merge(BandungData,dataku1120,by="id",all.x=T)
FinalData1120=merge.shp1120[order(merge.shp1120$order),]

#2020 Dec
dataku1220=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB12_20,
Kecamatan=Bandung$NAME_3)
merge.shp1220=merge(BandungData,dataku1220,by="id",all.x=T)
FinalData1220=merge.shp1020[order(merge.shp1220$order),]


#Final
bulan=rep(month.name,each=18299)
FinalBDGData20=rbind(FinalData120,FinalData220,FinalData320,
		FinalData420,FinalData520,FinalData620,
		FinalData720,FinalData820,FinalData920,
		FinalData1020,FinalData1120,FinalData1220)
FinalBDGData20$Bulan=bulan
str(FinalBDGData20)


FinalBDGData20$Bulan <- factor(FinalBDGData20$Bulan, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
#Ggplot Map
nlcr=5
plotclr=brewer.pal(nlcr,"Reds")

p= ggplot(data=FinalBDGData20)+
geom_polygon(aes(x=long,y=lat,group=group,fill=ST),color="black",size=0.25)+
geom_text(aes(x=long,y=lat,label = id),data = Koordinate, size = 3, hjust = 0.5)+
facet_wrap(.~Bulan,scales="free",ncol=3)+
scale_fill_gradientn(colours=plotclr)+
theme_bw(base_size=15)+
ylab("")+
xlab("")+
labs(fill = "ST")+
theme(legend.position="right")+
theme(text =element_text(size=12)) +
guides(fill = guide_colorbar(title.position = "top",title.vjust = 1, frame.colour = "black",barwidth = 1,barheight = 10))

#SIMPAN
ggsave("plot2020fix.png", plot = p, width = 12, height = 8, dpi = 300)

#PETA 2021
#Setting data
#2021 Jan
dataku121=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB1_21,
Kecamatan=Bandung$NAME_3)
merge.shp121=merge(BandungData,dataku121,by="id",all.x=T)
FinalData121=merge.shp121[order(merge.shp123$order),]

#2021 Feb
dataku221=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB2_21,
Kecamatan=Bandung$NAME_3)
merge.shp221=merge(BandungData,dataku221,by="id",all.x=T)
FinalData221=merge.shp221[order(merge.shp221$order),]

#2021 Mar
dataku321=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB3_21,
Kecamatan=Bandung$NAME_3)
merge.shp321=merge(BandungData,dataku321,by="id",all.x=T)
FinalData321=merge.shp321[order(merge.shp321$order),]

#2021 April
dataku421=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB4_21,
Kecamatan=Bandung$NAME_3)
merge.shp421=merge(BandungData,dataku421,by="id",all.x=T)
FinalData421=merge.shp421[order(merge.shp421$order),]

#2021 Mei
dataku521=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB5_21,
Kecamatan=Bandung$NAME_3)
merge.shp521=merge(BandungData,dataku521,by="id",all.x=T)
FinalData521=merge.shp521[order(merge.shp521$order),]

#2021 Juni
dataku621=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB6_21,
Kecamatan=Bandung$NAME_3)
merge.shp621=merge(BandungData,dataku621,by="id",all.x=T)
FinalData621=merge.shp621[order(merge.shp621$order),]

#2021 Juli
dataku721=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB7_21,
Kecamatan=Bandung$NAME_3)
merge.shp721=merge(BandungData,dataku721,by="id",all.x=T)
FinalData721=merge.shp721[order(merge.shp721$order),]

#2021 Agustus
dataku821=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB8_21,
Kecamatan=Bandung$NAME_3)
merge.shp821=merge(BandungData,dataku821,by="id",all.x=T)
FinalData821=merge.shp821[order(merge.shp821$order),]

#2021 September
dataku921=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB9_21,
Kecamatan=Bandung$NAME_3)
merge.shp921=merge(BandungData,dataku921,by="id",all.x=T)
FinalData921=merge.shp921[order(merge.shp921$order),]

#2021 Okt
dataku1021=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB10_21,
Kecamatan=Bandung$NAME_3)
merge.shp1021=merge(BandungData,dataku1021,by="id",all.x=T)
FinalData1021=merge.shp1021[order(merge.shp1021$order),]

#2021 Nov
dataku1121=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB11_21,
Kecamatan=Bandung$NAME_3)
merge.shp1121=merge(BandungData,dataku1121,by="id",all.x=T)
FinalData1121=merge.shp1121[order(merge.shp1121$order),]

#2021 Dec
dataku1221=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB12_21,
Kecamatan=Bandung$NAME_3)
merge.shp1221=merge(BandungData,dataku1221,by="id",all.x=T)
FinalData1221=merge.shp1021[order(merge.shp1221$order),]


#Final
bulan=rep(month.name,each=18299)
FinalBDGData21=rbind(FinalData121,FinalData221,FinalData321,
		FinalData421,FinalData521,FinalData621,
		FinalData721,FinalData821,FinalData921,
		FinalData1021,FinalData1121,FinalData1221)
FinalBDGData21$Bulan=bulan
str(FinalBDGData21)


FinalBDGData21$Bulan <- factor(FinalBDGData21$Bulan, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
#Ggplot Map
nlcr=5
plotclr=brewer.pal(nlcr,"Reds")

p= ggplot(data=FinalBDGData21)+
geom_polygon(aes(x=long,y=lat,group=group,fill=ST),color="black",size=0.25)+
geom_text(aes(x=long,y=lat,label = id),data = Koordinate, size = 3, hjust = 0.5)+
facet_wrap(.~Bulan,scales="free",ncol=3)+
scale_fill_gradientn(colours=plotclr)+
theme_bw(base_size=15)+
ylab("")+
xlab("")+
labs(fill = "ST")+
theme(legend.position="right")+
theme(text =element_text(size=12)) +
guides(fill = guide_colorbar(title.position = "top",title.vjust = 1, frame.colour = "black",barwidth = 1,barheight = 10))

#SIMPAN
ggsave("plot2021fix.png", plot = p, width = 12, height = 8, dpi = 300)


#PETA 2022
# January
dataku122=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB1_22,
Kecamatan=Bandung$NAME_3)
merge.shp122=merge(BandungData,dataku122,by="id",all.x=T)
FinalData122=merge.shp122[order(merge.shp122$order),]

# February
dataku222=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB2_22,
Kecamatan=Bandung$NAME_3)
merge.shp222=merge(BandungData,dataku222,by="id",all.x=T)
FinalData222=merge.shp222[order(merge.shp222$order),]

# March
dataku322=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB3_22,
Kecamatan=Bandung$NAME_3)
merge.shp322=merge(BandungData,dataku322,by="id",all.x=T)
FinalData322=merge.shp322[order(merge.shp322$order),]

# April
dataku422=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB4_22,
Kecamatan=Bandung$NAME_3)
merge.shp422=merge(BandungData,dataku422,by="id",all.x=T)
FinalData422=merge.shp422[order(merge.shp422$order),]

# May
dataku522=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB5_22,
Kecamatan=Bandung$NAME_3)
merge.shp522=merge(BandungData,dataku522,by="id",all.x=T)
FinalData522=merge.shp522[order(merge.shp522$order),]

# June
dataku622=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB6_22,
Kecamatan=Bandung$NAME_3)
merge.shp622=merge(BandungData,dataku622,by="id",all.x=T)
FinalData622=merge.shp622[order(merge.shp622$order),]

# July
dataku722=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB7_22,
Kecamatan=Bandung$NAME_3)
merge.shp722=merge(BandungData,dataku722,by="id",all.x=T)
FinalData722=merge.shp722[order(merge.shp722$order),]

# August
dataku822=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB8_22,
Kecamatan=Bandung$NAME_3)
merge.shp822=merge(BandungData,dataku822,by="id",all.x=T)
FinalData822=merge.shp822[order(merge.shp822$order),]

# September
dataku922=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB9_22,
Kecamatan=Bandung$NAME_3)
merge.shp922=merge(BandungData,dataku922,by="id",all.x=T)
FinalData922=merge.shp922[order(merge.shp922$order),]

# October
dataku1022=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB10_22,
Kecamatan=Bandung$NAME_3)
merge.shp1022=merge(BandungData,dataku1022,by="id",all.x=T)
FinalData1022=merge.shp1022[order(merge.shp1022$order),]

# November
dataku1122=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB11_22,
Kecamatan=Bandung$NAME_3)
merge.shp1122=merge(BandungData,dataku1122,by="id",all.x=T)
FinalData1122=merge.shp1122[order(merge.shp1122$order),]

# December
dataku1222=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB12_22,
Kecamatan=Bandung$NAME_3)
merge.shp1222=merge(BandungData,dataku1222,by="id",all.x=T)
FinalData1222=merge.shp1222[order(merge.shp1222$order),]


# Final
bulan=rep(month.name,each=18299)
FinalBDGData22=rbind(FinalData122,FinalData222,FinalData322,
		FinalData422,FinalData522,FinalData622,
		FinalData722,FinalData822,FinalData922,
		FinalData1022,FinalData1122,FinalData1222)
FinalBDGData22$Bulan=bulan
str(FinalBDGData22)


FinalBDGData22$Bulan <- factor(FinalBDGData22$Bulan, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
# Ggplot Map
nlcr=5
plotclr=brewer.pal(nlcr,"Reds")

p= ggplot(data=FinalBDGData22)+
geom_polygon(aes(x=long,y=lat,group=group,fill=ST),color="black",size=0.25)+
geom_text(aes(x=long,y=lat,label = id),data = Koordinate, size = 3, hjust = 0.5)+
facet_wrap(.~Bulan,scales="free",ncol=3)+
scale_fill_gradientn(colours=plotclr)+
theme_bw(base_size=15)+
ylab("")+
xlab("")+
labs(fill = "ST")+
theme(legend.position="right")+
theme(text =element_text(size=12)) +
guides(fill = guide_colorbar(title.position = "top",title.vjust = 1, frame.colour = "black",barwidth = 1,barheight = 10))

# Save the plot
ggsave("plot2022fix.png", plot = p, width = 12, height = 8, dpi = 300)

#PETA 2023
#Setting data
#2023 Jan
dataku123=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB1_23,
Kecamatan=Bandung$NAME_3)
merge.shp123=merge(BandungData,dataku123,by="id",all.x=T)
FinalData123=merge.shp123[order(merge.shp123$order),]

#2023 Feb
dataku223=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB2_23,
Kecamatan=Bandung$NAME_3)
merge.shp223=merge(BandungData,dataku223,by="id",all.x=T)
FinalData223=merge.shp223[order(merge.shp223$order),]

#2023 Mar
dataku323=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB3_23,
Kecamatan=Bandung$NAME_3)
merge.shp323=merge(BandungData,dataku323,by="id",all.x=T)
FinalData323=merge.shp323[order(merge.shp323$order),]


#2024 April
dataku423=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB4_23,
Kecamatan=Bandung$NAME_3)
merge.shp423=merge(BandungData,dataku423,by="id",all.x=T)
FinalData423=merge.shp423[order(merge.shp423$order),]

#2023 Mei
dataku523=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB5_23,
Kecamatan=Bandung$NAME_3)
merge.shp523=merge(BandungData,dataku523,by="id",all.x=T)
FinalData523=merge.shp523[order(merge.shp523$order),]

#2023 Juni
dataku623=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB6_23,
Kecamatan=Bandung$NAME_3)
merge.shp623=merge(BandungData,dataku623,by="id",all.x=T)
FinalData623=merge.shp623[order(merge.shp623$order),]

#2023 Juli
dataku723=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB7_23,
Kecamatan=Bandung$NAME_3)
merge.shp723=merge(BandungData,dataku723,by="id",all.x=T)
FinalData723=merge.shp723[order(merge.shp723$order),]

#2023 Agustus
dataku823=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB8_23,
Kecamatan=Bandung$NAME_3)
merge.shp823=merge(BandungData,dataku823,by="id",all.x=T)
FinalData823=merge.shp823[order(merge.shp823$order),]

#2023 September
dataku923=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB9_23,
Kecamatan=Bandung$NAME_3)
merge.shp923=merge(BandungData,dataku923,by="id",all.x=T)
FinalData923=merge.shp923[order(merge.shp923$order),]

#2023 Okt
dataku1023=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB10_23,
Kecamatan=Bandung$NAME_3)
merge.shp1023=merge(BandungData,dataku1023,by="id",all.x=T)
FinalData1023=merge.shp1023[order(merge.shp1023$order),]

#2023 Nov
dataku1123=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB11_23,
Kecamatan=Bandung$NAME_3)
merge.shp1123=merge(BandungData,dataku1123,by="id",all.x=T)
FinalData1123=merge.shp1123[order(merge.shp1123$order),]

#2023 Dec
dataku1223=data.frame(id=id,x=Coordinate[,1],
y=Coordinate[,2],ST=STB12_23,
Kecamatan=Bandung$NAME_3)
merge.shp1223=merge(BandungData,dataku1223,by="id",all.x=T)
FinalData1223=merge.shp1023[order(merge.shp1223$order),]


#Final
bulan=rep(month.name,each=18299)
FinalBDGData=rbind(FinalData123,FinalData223,FinalData323,
		FinalData423,FinalData523,FinalData623,
		FinalData723,FinalData823,FinalData923,
		FinalData1023,FinalData1123,FinalData1223)
FinalBDGData$Bulan=bulan
str(FinalBDGData)
head(FinalBDGData)

FinalBDGData$Bulan <- factor(FinalBDGData$Bulan, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
#Ggplot Map
nlcr=5
plotclr=brewer.pal(nlcr,"Reds")

p= ggplot(data=FinalBDGData)+
geom_polygon(aes(x=long,y=lat,group=group,fill=ST),color="black",size=0.25)+
geom_text(aes(x=long,y=lat,label = id),data = Koordinate, size = 3, hjust = 0.5)+
facet_wrap(.~Bulan,scales="free",ncol=3)+
scale_fill_gradientn(colours=plotclr)+
theme_bw(base_size=15)+
ylab("")+
xlab("")+
labs(fill = "ST")+
theme(legend.position="right")+
theme(text =element_text(size=12)) +
guides(fill = guide_colorbar(title.position = "top",title.vjust = 1, frame.colour = "black",barwidth = 1,barheight = 10))

#SIMPAN
ggsave("plot.png", plot = p, width = 12, height = 8, dpi = 300)
