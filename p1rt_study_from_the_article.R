library(latex2exp); library(maps); library(sf); library(magrittr); library(rnaturalearth)
library(chronosphere); library(grDevices); library(dplyr);library(RColorBrewer)

######################### Source functions ##############################################

# The file with the needed functions is available in the github project folder
source("p1rt_functions_for_the_study")

##########################################################################################
############################# Get data  ##################################################

# pr_cmip6_yearmax and df_lonlat files are available in the github project folder
pr_cmip6_yearmax <- readRDS("df_pr_cmip56_yearmax.rds")
df_lonlat <- readRDS("df_lonlat.rds")

##########################################################################################
#                           EXEMPLE 1: Study of 1 gridpoint                              #
##########################################################################################

################# Gridpoint studied in the article    ####################################
gridpoint <- 1821
model <- "IPSL-CM6A-LR(CM6)"
run <- "r1i1p1f1"
varcmip <- "cmip6"
lon <- df_lonlat[gridpoint ,1]
lat <- df_lonlat[gridpoint ,2]
years <- seq(1850, 2100,by=1)
mat_X_Z <- traj_from_data(variable.df=pr_cmip6_yearmax,
                          grid_points=gridpoint,
                          model.choice=model,
                          run.choice=run,
                          var=varcmip)
matx <- mat_X_Z$matX
matz <- mat_X_Z$matZ
tt <- mat_X_Z$time.vec

########################### P1rt vs t  ####################################################
h <- 60.5
GZ<-matGZ_func(matx,matz)
param <- weibullGMM_NonStationaire(GZ, tt, tt, h = h, kern=dEpan, truevalues=NULL)
lambda.t <- param$lambdahat
k.t <- param$khat
############################# decadal records ############################################
r <-10 # decadal
zalpha <- qnorm(0.95)
covna <- matcovNA_alone(matx,matz,tt,tt,h,graphiques=FALSE)
covnb <- matcovNB_alone(matx,matz,tt,tt,h,method="MC",m=2000,graphiques=FALSE)
covn <- covna + covnb
G_emp <- ecdf(matx)
GmZ <- G_emp(matz)
varestimpr<-varp1rfar_t(r,covn, matx, matz, GmZ,tt,tt,h)
varpr <-varestimpr$varp1r_t
p1rthat <- varestimpr$p1r_t[,1]
std <- sqrt(varpr)
lowerbn <- p1rthat -(zalpha*std)
upperbn<- p1rthat + (zalpha*std)
par(mfrow = c(1, 1),mar=c(5,5,5,5))
plot(years,p1rthat,type="l", ylim=c(0,1),xlab = "years",
     ylab="probability",
     main = "(a) Decadal record probability",
     col="red",cex.main=1.3, cex.lab=1.5)
polygon(c(years,rev(years)),
        c(lowerbn,rev(upperbn)),col = "lightcyan", border = FALSE)
lines(years,p1rthat,col="red")
v <- c(expression(widehat(p["1,10"])~"(t)"),expression(p["0,10"]~"(t)"))
legend(locator(1),legend=v,col=c("red","blue"),lty=c(1,2),cex=1.2,y.intersp=1.3)
l1 <- c(parse(text=sprintf('paste(Lat,\' = %s\')',lat)),
        parse(text=sprintf('paste(Lon,\' = %s\')',lon)),
        parse(text=sprintf('paste(Model,\' = %s\')',"IPSL-CM5A-MR")))
legend(locator(1),legend=l1,bty='n',cex=1.2) # Click on plot 
which(lowerbn > 0.1)[1]+1849
abline(h=0.1, col="blue", lty=2)
abline(v=2002, col="gray", lty=2)

############################# centennial records ############################################
r <-100 # centennial
varestimpr<-varp1rfar_t(r,covn, matx, matz, GmZ,tt,tt,h)
varpr <-varestimpr$varp1r_t
p1rthat <- varestimpr$p1r_t[,1]
std <- sqrt(varpr)
lowerbn <- p1rthat -(zalpha*std)
upperbn<- p1rthat + (zalpha*std)
plot(years,p1rthat,type="l", ylim=c(0,1),xlab = "years",
     ylab="probability",
     main = "(b) Centennial record probability",
     col="red",cex.main=1.3, cex.lab=1.5)
polygon(c(years,rev(years)),
        c(lowerbn,rev(upperbn)),col = "lightcyan", border = FALSE)
lines(years,p1rthat,col="red")
v <- c(expression(widehat(p["1,100"])~"(t)"),expression(p["0,100"]~"(t)"))
legend(locator(1),legend=v,col=c("red","blue"),lty=c(1,2),cex=1.2,y.intersp=1.3)
l1 <- c(parse(text=sprintf('paste(Lat,\' = %s\')',lat)),
        parse(text=sprintf('paste(Lon,\' = %s\')',lon)),
        parse(text=sprintf('paste(Model,\' = %s\')',"IPSL-CM5A-MR")))
legend(locator(1),legend=l1,bty='n',cex=1.2) # click on plot
abline(h=0.01, col="blue", lty=2)
which(lowerbn > 0.01)[1]+1849 
abline(v=2046, col="gray", lty=2)


##########################################################################################
#               EXEMPLE 2 :  Maps  using cmip6  IPSL-CM6A-LR model                       #                             #
##########################################################################################

############################# Option (1) : Estimation of p_{1,r}(t) #######################
pr_cmip6_yearmax <- readRDS("/Users/pgonzale/Documents/datasets/df_pr_cmip56_yearmax.rds")
df_lonlat <- readRDS("/Users/pgonzale/Documents/datasets/df_lonlat.rds")

model.choice <- "IPSL-CM6A-LR(CM6)"
run.choice <- "r1i1p1f1"
var <- "cmip6"
trajmap <- traj_map_from_data(variable.df = pr_cmip6_yearmax,model.choice = model.choice,run.choice = run.choice,var = var )
matx <- trajmap$matrix_x
matz <- trajmap$matrix_z
tt <- c(1:dim(matz)[1])
h <- 60.5
r <- 10
p1r.mat <- matrix(0,nrow = 251 ,ncol= 2592)
var.mat <- matrix(0,nrow = 251 ,ncol= 2592)
for(i in 1:dim(matz)[2]){
  traj_map_p1r <- matrix(0,nrow = 251 ,ncol= 2592)
  traj_map_var <- matrix(0,nrow = 251 ,ncol= 2592)
  covna <- matcovNA_alone(matx[,i],matz[,i],tt,tt,h,graphiques=FALSE)
  covnb <- matcovNB_alone(matx[,i],matz[,i],tt,tt,h,method="MC",graphiques=FALSE)
  covn<-covna+covnb
  G_emp <- ecdf(matx[,i])
  GmZ <- G_emp(matz[,i])
  varestimpr<-varp1rfar_t(r,covn, matx[,i], matz[,i], GmZ,tt,tt,h)
  varpr <-varestimpr$varp1r_t
  p1rthat <- varestimpr$p1r_t[,1]
  p1r.mat[,i] <- as.numeric(p1rthat)
  var.mat[,i] <- as.numeric(varpr)
}
# p1r.mat : \hat{p}_{1,r}(t)
# var.mat : variance of the estimator \hat{p}_{1,r}(t)

##########################################################################################
#               EXEMPLE 2 : Maps  using cmip6  IPSL-CM6A-LR model                        #                             #
##########################################################################################

#################### Option (2) : Visualise results from the article  #####################
#                            [ import files from folder xxxxx ]
zalpha <- qnorm(0.95)
p1r.mat <- matrix(0,nrow = 251 ,ncol= 2592)
var.mat <- matrix(0,nrow = 251 ,ncol= 2592)
lowerbound.mat <- matrix(0,nrow = 251 ,ncol= 2592)
upperbound.mat <- matrix(0,nrow = 251 ,ncol= 2592)
for (i in 1:2592){
  # read files
  filename_p1r <- paste("p1rvarcmip6/p1r_",i,"_cmip6_r1i1p1f1.rds",sep ="")
  filename_var <- paste("p1rvarcmip6/var_",i,"_cmip6_r1i1p1f1.rds",sep ="")
  p <- readRDS(filename_p1r)
  v <- readRDS(filename_var)
  p1r.mat[,i] <- p
  var.mat[,i] <- v
  std.mat <- sqrt(var.mat)
  lowerbound.mat[,i] <- p1r.mat[,i] - (zalpha*std.mat[,i])
  upperbound.mat[,i] <- p1r.mat[,i] + (zalpha*std.mat[,i])
}

############################# Identify emergence times ###########################################
r <- 10 
index_signal <- rep(0, dim(p1r.mat)[2])
for (i in 1:dim(p1r.mat)[2] ){
  if ( isTRUE(p1r.mat[251,i]>1/r) ){
    year <- tail(which(lowerbound.mat[,i] < 1/r),n=1) + 1 + 1849
    if (identical(year, numeric(0))){year <- 1850}
    if (year == 2101){year <- NA}
  }
  if (isTRUE(p1r.mat[251,i]<1/r)){
    year <- tail(which(upperbound.mat[,i] > 1/r),n=1) + 1 + 1849
    if (identical(year, numeric(0))){year <- 1850}
    if (year == 2101){year <- NA}
  }
  index_signal[i] <- year
}

####################### Probability ratio of decadal records in 2050 ####################################

target_crs <- st_crs("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
sf_use_s2(FALSE)
world1 <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()
offset <- 180 - 0
polygon <- st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)
world2 <- world1 %>% st_difference(polygon)
world3 <- sf::st_transform(world2, crs = target_crs)
par(mfrow=c(2,1))
#nf <- layout( matrix(c(1,2), ncol=1) )
p1r_vec <- p1r.mat[201,]
y <- as.vector(matrix(p1r_vec))
unique.lat <- seq(-87.5, 87.5, by = 5)
unique.lon <- seq(-177.5, 177.5, by = 5)
lat <- rep(unique.lat, each = length(unique.lon))
lon <- rep(unique.lon, length(unique.lat))
for (j in 1: length(y)){
  index_signal[is.na(index_signal)] <- 3000
  if (isTRUE(2050  < index_signal[j]))
  {
    y[j] <- NA
  }
}
y <- matrix(y, nrow = length(unique.lon))
dat <- data.frame(lat = lat, lon = lon, y = y)
data.mat <- y
df <- st_as_sf(
  x = dat,
  coords = c("lon", "lat"),
  crs = "+proj=lonlat"
)
df_proj <- st_transform(df, crs = target_crs)
xy_proj <- st_coordinates(df_proj)
xmat <- matrix(xy_proj[, 1], nrow(data.mat))
ymat <- matrix(xy_proj[, 2], nrow(data.mat))
colsp <- brewer.pal(11, "BrBG")
colsp <- colorRampPalette(colsp)(64)
r_br <- c(1/7,1/6,1/5,1/4,1/3,1/2,1,2,3,4,5,6,7)
r_breaks1 <- seq(1/7,1,length.out=33)[1:32]
r_breaks2 <- seq(1,7,length.out=33)
r_breaks <- log(c(r_breaks1,r_breaks2))
r_brlab <- c("1/7","1/6","1/5","1/4","1/3","1/2",1,2,3,4,5,6,7)
p1rpzerp_vec <-p1r_vec*10
dat <- data.frame(lat = lat, lon = lon, y = y)
data.mat <- y
dataplot <- data.mat * 10
dataplot <- replace(dataplot, dataplot<1/7, 7)
dataplot <- replace(dataplot, dataplot>7, 7)
breaks <- c(-max(p1rpzerp_vec), -6, -5,  -4, -3, -1.5, 1.5, 3, 4, 5, 6, max(p1rpzerp_vec))
fields::image.plot(
  xmat, ymat, log(dataplot),
  xlab = "", ylab = "",
  asp = 1, xaxt = "n", yaxt = "n", bty = "n",
  col = colsp,
  breaks = r_breaks,
  axis.args = list(at=log(r_br),label=r_brlab),
    smallplot = c(.7, .72, .25, .8),
  legend.args = list( text = TeX("$\\widehat{p_{1,10}}\\,(2050)\\, x\\, 1/10$"),
                      cex = 0.7,
                      side = 3,
                      line = .5)
)
plot(st_geometry(world3), add = TRUE, graticule = st_crs(4326), axes = TRUE, border = "white", lwd = 3)
plot(st_geometry(world3), add = TRUE, graticule = st_crs(4326), axes = TRUE)
title("a) Probability ratio of decadal records in 20250",cex.main=1)

############################# Time of emergence of decadak records #########################################
dataplot <- index_signal
data.mat <- matrix(data=dataplot,nrow=72,ncol=36)
colsp <- rev(brewer.pal(8, "Spectral"))
colsp <- c(colsp,"gray","gray" )
r_breaks <- seq(1850,2100,length.out=11)
fields::image.plot(
  xmat, ymat, data.mat,
  xlab = "", ylab = "",
  asp = 1, xaxt = "n", yaxt = "n", bty = "n",
  col = rev(colsp),
  breaks = r_breaks,
  smallplot = c(.7, .72, .25, .8),
  legend.args = list( text = "year",
                      cex = 1,
                      side = 3,
                      line = .5)
)
plot(st_geometry(world3), add = TRUE, graticule = st_crs(4326), axes = TRUE, border = "white", lwd = 3)
plot(st_geometry(world3), add = TRUE, graticule = st_crs(4326), axes = TRUE)
title("b) Time of emergence of decadal records",cex.main=1)


