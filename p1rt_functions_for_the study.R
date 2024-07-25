library(dplyr)
library(stringr)
library(stats4); library(gmm); library(stats); library(np); library(EWGoF)

##########################################################################################
#        Functions for the estimation of p_{1,r}(t) as in exemple file                  #
##########################################################################################
# See file p1rt_study_from_article.R

matGZ_func<- function(matX,matZ){
  
  ## This function of computes G(Ztj) for a given trajectory X and Z, where the imput variables are two matrix
  ## matX, matZ, conteining each one multiple trajectories X and Z
  ## This allow us to compute G(Ztj) multiple times for differents trajectories.
  ## matX and matZ must have the same number of columns
  ##
  ## Input :
  ## - matrices matX and matZ of trajectories X and Z, where each column is a trayectory
  ##
  ## Output :
  ## - matrix of dimentions ( dim(matX)[1] x J ) containong G(Ztj) values
  ##

  matX <- as.matrix(matX)
  matZ <- as.matrix(matZ)
  dimnsZ=dim(matZ)
  matGm <- matrix(nrow=dimnsZ[1],ncol=dimnsZ[2])
  for (j in 1:dimnsZ[2]){
    X <- matX[,j]
    Z <- matZ[,j]
    G_empirique<-ecdf(X)
    matGm[,j] <- G_empirique(Z)
  }
  return(matGm)
}

weibullGMM_NonStationaire <- function(matGm, tt, t_eval, h, kern=dEpan, truevalues=NULL){
  
  ## This function computes the estimates of multiple lambda_t and k_t from the output of the
  ## function matGZ_func(), via M-estimation
  ##
  ## Input :
  ## - matGm : a matrix issued from matGZ_func(), containing the values \hat{G}m(Z_tj) , where each column is associated
  ##   to a different trayectory X_t and Z_t.
  ## - tt: vector tt of length(Z) containing the Z's time steps
  ## - t_eval: evaluation vector (in practice t_eval = tt)
  ## - truevalues : an optional vector, which can contain appropriate starting values for the GMM estimation process of
  ##   lambda_1 and k_1, i.e, the estimations for the first time step.
  ##   if the starting values ar not given, they will be chosen as [1,1]
  ##   for the next times steps t=2, .... the starting value at t will be the optimal found in t-1.
  ## - kernel bandwith value h
  ##
  ## Output :a list of 2 elements
  ##   the first one ($lambdahat) containing the GMM estimate of the scale lambda_t parameter for each column of matGm
  ##   (ie of matX and matZ) and for each time step t.
  ##   the first one ($khat) containing the GMM estimate of the scale k_t parameter for each column of matGm
  ##   (ie of matX and matZ) and for each time step t.
  ##
  ## Requires : the gmm() routine from the gmm package and the utilitary function function_gmm_noyaux()
  ##
  ## In gmm(), the "nlminb" optimization routine is used.
  matGm <- as.matrix(matGm)
  Nligne <- dim(matGm)[1]
  Ncol <- dim(matGm)[2]
  lambdahat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  khat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  if ( is.null(truevalues) ){
    startvalueslambd <- 1
    startvaluesk <- 1
  } else {
    startvalueslambd <- truevalues[1] 
    startvaluesk <- truevalues[2]
  }
  for (j in 1:Ncol){
    Gnhat = matGm[,j]
    lambdahat.vec = lambdahat.mat[,j]
    khat.vec = khat.mat[,j]
    startvalueslambd_t <- startvalueslambd
    startvaluesk_t<- startvaluesk
    for (i in 1:Nligne){
      startval_t=as.numeric(c(startvalueslambd_t,startvaluesk_t))
      fg_noyaux <- function(theta,vecx){function_gmm_noyaux(theta,vecx,index=i,tt.vec=tt,t_eval.vec=t_eval,bandwidth=h)}
      EGMMweibull_NonStationary <- gmm(g=fg_noyaux,
                                       x=Gnhat,
                                       t0=startval_t,
                                       optfct="nlminb",
                                       lower=c(10^(-5),10^(-5)),upper=c(Inf,Inf),
                                       onlyCoefficients = TRUE
      )
      lambdahat.vec[i] <- EGMMweibull_NonStationary$coefficients[1]
      khat.vec[i] <- EGMMweibull_NonStationary$coefficients[2]
      startvalueslambd_t <- lambdahat.vec[i]
      startvaluesk_t <- khat.vec[i]
    }
    lambdahat.mat[,j] <- lambdahat.vec
    khat.mat[,j] <- khat.vec
  }
  return( list("lambdahat"=lambdahat.mat,"khat"=khat.mat) )
  
}

matcovNA_alone <- function(X,Z,tt,t_eval,h,graphiques=TRUE){
  # weighed variance of NA
  Z<-as.numeric(Z)
  X <- as.numeric(X)
  J <- length(Z)
  p12_hat.vec <- as.vector(p12_NonPar(X,Z,tt,t_eval,h)$p12.mat)
  p13_hat.vec<- p13_NonPar(X,Z,tt,tt,h)
  G_emp <- ecdf(X);matGm<- G_emp(Z)
  tetha<-weibullGMM_NonStationaire(matGm, tt, t_eval, h, kern=dEpan, truevalues=NULL)
  lambda_h <- tetha$lambdahat
  k_h <- tetha$khat
  p14W_t <- p1rfarW_temps(lambda_h,k_h,matrix(4,ncol=1,nrow=J))$p1r
  p15W_t <- p1rfarW_temps(lambda_h,k_h,matrix(5,ncol=1,nrow=J))$p1r
  A2<-matcovNA_A2_jit(p12_hat.vec,p13_hat.vec)
  B2<-matcovNA_B2_jit(p13_hat.vec,p15W_t)
  C2<-matcovNA_C2_jit(p12_hat.vec,p13_hat.vec,p14W_t)
  Kh <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  list_unweighed_NA_t<- array(NA,c(2,2,J))
  list_unweighed_NA_t[1,1,] <- as.numeric(A2)
  list_unweighed_NA_t[1,2,] <- as.numeric(C2)
  list_unweighed_NA_t[2,1,] <- as.numeric(C2)
  list_unweighed_NA_t[2,2,] <- as.numeric(B2)
  list_weighed_matcov_N_t<- array(NA,c(2,2,J)) 
  for (index_t in 1:J){
    denom<- 1/J * sum(Kh[,index_t])
    list_khj_t <- rep(0,J)
    for (j in 1:dim(Kh)[1]){
      Khj<-Kh[j,index_t]
      list_khj_t[j]<- (Khj/denom)^2
    }
    list_weighed_matcov_N_t[1,1,index_t] <- 1/(J^2) * sum(list_khj_t*2*A2)
    list_weighed_matcov_N_t[1,2,index_t] <-  1/(J^2) *sum(list_khj_t*2*C2)
    list_weighed_matcov_N_t[2,1,index_t] <- 1/(J^2) *sum(list_khj_t*2*C2)
    list_weighed_matcov_N_t[2,2,index_t] <- 1/(J^2) *sum(list_khj_t*2*B2)
  }
  if (graphiques==TRUE){
    plot(tt,list_weighed_matcov_N_t[1,1,])
    plot(tt,list_weighed_matcov_N_t[1,2,])
    plot(tt,list_weighed_matcov_N_t[2,2,])
    plot(lambda_h)
    plot(k_h)
  }
  return(list_weighed_matcov_N_t)
}

matcovNB_alone <- function(X,Z,tt,t_eval,h,method="integration",m=2000,graphiques=TRUE){
  
  method <-  method # integration, MC, CVMC
  
  # weighed variance of NB
  
  X <- as.numeric(X)
  Z <- as.numeric(Z)
  J<-length(Z)
  I <- length(X)
  p12_hat.vec <- as.vector(p12_NonPar(X,Z,tt,t_eval,h)$p12.mat)
  p13_hat.vec<- p13_NonPar(X,Z,tt,tt,h)
  G_emp <- ecdf(X);matGm<- G_emp(Z)
  tetha<-weibullGMM_NonStationaire(matGm, tt, t_eval, h, kern=dEpan, truevalues=NULL)
  lambda_h <- tetha$lambdahat
  k_h <- tetha$khat
  uw_matcovNB_a<-matcovNB_aji(p12_hat.vec,lambda_h,k_h,method=method,m=m)
  uw_matcovNB_b <- matcovNB_bji(p13_hat.vec,lambda_h,k_h,method=method,m=m)
  uw_matcovNB_c <- matcovNB_cji(p12_hat.vec,p13_hat.vec,lambda_h,k_h,method=method,m=m)
  Kh <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  list_weighed_matcov_NB_t<- array(NA,c(2,2,J)) 
  for (index_t in 1:J){
    denom <- 1/J * sum(Kh[,index_t])
    inx <- 1
    list_Wji_t <- rep(0,31626)
    for (j in 1:dim(Kh)[1]){
      Khj<-Kh[j,index_t]
      for (i in 1:j){
        Khi <- Kh[i,index_t]
        KhjKhi <- (Khj * Khi)/(denom)^2
        list_Wji_t[inx]<- KhjKhi
        inx <- inx +1
      }
    }
    list_weighed_matcov_NB_t[1,1,index_t] <- (1/I) * 1/(J^2) * sum(list_Wji_t * 2*uw_matcovNB_a )
    list_weighed_matcov_NB_t[1,2,index_t] <- (1/I) * 1/(J^2) * sum(list_Wji_t * uw_matcovNB_c )
    list_weighed_matcov_NB_t[2,1,index_t] <- (1/I) * 1/(J^2) * sum(list_Wji_t * uw_matcovNB_c)
    list_weighed_matcov_NB_t[2,2,index_t] <- (1/I) * 1/(J^2) * sum(list_Wji_t * 2*uw_matcovNB_b)
  }
  if (graphiques==TRUE){
    plot(tt,list_weighed_matcov_NB_t[1,1,])
    plot(tt,list_weighed_matcov_NB_t[1,2,])
    plot(tt,list_weighed_matcov_NB_t[2,2,])
  }
  return(list_weighed_matcov_NB_t)
}

varp1rfar_t <- function(r,matcovN_t, X.vec, Z.vec, GmZ,tt,t_eval,h){
  # Var{p1r_t} and Var{far_t}
  J <- length(Z.vec)
  theta_t <-weibullGMM_NonStationaire (GmZ, tt, t_eval, h, kern=dEpan, truevalues=NULL)
  lam_t <- theta_t[[1]]
  k_t <- theta_t[[2]]
  p1rW_t <- p1rfarW_temps(as.matrix(lam_t),as.matrix(k_t),matrix(r,ncol=1,nrow=J))
  
  Jacov12 <- jacobianFunctiong12(lam_t,k_t)
  Jacovrminus1 <- jacobianFunctiongrminus1(lam_t,k_t,rep(r,J))
  
  list_variancep1r_t <- rep(0,J)
  list_variancefar_t <- rep(0,J)
  for (i in 1:J){
    Jacov12_inv <- solve(Jacov12[[1]][[i]])
    Jacov12T_inv <- solve(t(Jacov12[[1]][[i]]))
    list_variancep1r_t[i] <- Jacovrminus1[[i]]%*% Jacov12_inv %*% matcovN_t[,,i] %*% Jacov12T_inv %*% t(Jacovrminus1[[i]])
  }
  return (list("varp1r_t"=list_variancep1r_t,"p1r_t"=p1rW_t$p1r))
}


##########################################################################################
#      Custom-made functions for the specific study  (usend in the application file)     #
##########################################################################################

traj_from_data <- function(variable.df, grid_points, model.choice, run.choice,var="cmip6"){ # tmax or pr
  
  # just extract trayectories: matx, matz and tt
  
  if (var=="cmip6"){
    hist_fin <- 2014
    rcp85_0 <- 2015
  }
  
  if (var=="cmip5"){
    hist_fin <- 2005
    rcp85_0 <- 2006
  }
  
  Z_historical <- variable.df %>%
    dplyr::select (institute,model, experiment, run, year, one_of(str_c(grid_points))) %>%
    filter(experiment == "historical" & model == model.choice & run == run.choice & between(year, 1850, hist_fin))  %>%
    arrange(year) %>%
    dplyr::select(!c(institute,model,experiment,run, year))
  Z_historical<-as.data.frame(Z_historical)
  Z_rcp85<- variable.df %>%
    dplyr::select (institute, model, experiment, run, year, one_of(str_c(grid_points))) %>%
    filter(experiment == "rcp85" & model == model.choice & run == run.choice & between(year, rcp85_0, 2100))  %>%
    arrange(year) %>%
    dplyr::select(!c(institute,model,experiment,run, year))
  Z_rcp85<-as.data.frame(Z_rcp85)
  matz <- bind_rows(Z_historical,Z_rcp85)
  matz <-as.matrix(matz)
  
  matx<- variable.df %>%
    dplyr::select (institute, model, experiment, run, year, one_of(str_c(grid_points))) %>%
    filter(experiment == "historicalNat" & model == model.choice & run == run.choice)  %>%
    arrange(year) %>% dplyr::select(!c(institute,model,experiment,run, year))
  matx<-as.matrix(matx)
  
  tt<-c(1:dim(matz)[1])
  
  return(list("matZ"=matz,"matX"=matx,"time.vec"=tt))
}

traj_map_from_data <-function(variable.df, model.choice, run.choice, varname){
  
  ## for a given variable it creates the matrix matlam and matk,
  ## which contains the the values of the variables for each gridpoint
  ## and for each t from 1850 to 2100
  ##
  ## Input
  ## variable.df: the data base containing Tmax trayectories for differents models, scenarios and runs,
  ##          e.g. tmax_cmip6_yearmax.rds"
  ## model.choice: the name of the model in between "", e.g. IPSL-CM6A-LR"
  ## run.choice:  the name of the tun in between "", e.g. "r1i1p1f1"
  ## varname: name of the varible of study in between "", e.g. "tmax" or "pr"
  ##
  ## Output
  ## two files.rds, one with de values of lambda and the other with the values of k. the file's names are
  ## "matlam_",model.choice,"_",run.choice,".rds"
  ##
  ##
  ##
  
  stringYear<-toString(model.choice)
  stringYear<-toString(run.choice)
  
  Z_historical <- variable.df %>%
    filter(experiment == "historical" & model == model.choice & run == run.choice & between(year, 1850, 2014))  %>%
    arrange(year) %>%dplyr::select(!c(institute,model,experiment,run,year))
  Z_historical<-as.data.frame(Z_historical)
  Z_rcp85<- variable.df %>%
    filter(experiment == "rcp85" & model == model.choice & run == run.choice & between(year, 2015, 2100))  %>%
    arrange(year) %>%
    dplyr::select(!c(institute,model,experiment,run,year))
  Z_rcp85<-as.data.frame(Z_rcp85)
  matz <- bind_rows(Z_historical,Z_rcp85)
  matz <-as.matrix(matz)
  
  matx<- variable.df %>%
    filter(experiment == "historicalNat" & model == model.choice & run == run.choice)  %>%
    arrange(year) %>% dplyr::select(!c(institute,model,experiment,run, year))
  matx<-as.matrix(matx)
  
  tt<-c(1:dim(matz)[1])
  
  return(list("matrix_x"=matx,"matrix_z"=matz))
  
}


##########################################################################################
#         Utilitary Functions                                                            #  
##########################################################################################

dEpan <- function(x){
  ## Function of Epanechnikov density distribution
  k <- (3/4)* (1-x^2)
  k[-1>x] <- 0
  k[x>1] <- 0
  return (k)
}

function_gmm_noyaux<-function(theta,vecx,index,tt.vec,t_eval.vec,bandwidth,kern=dEpan){
  
  ## Utilitary function required for the M-estimation of (lambda_t,k_t)
  ## in function weibullGMM_NonStationaire() and weibullGMM_NonStationaire_startval_1()
  ##
  ## Input :
  ## - theta : matrix of dimesion 2 x J containing the tentative values of lambda_t and k_t
  ## - vecx  : a vector of size n supposed to contain the values \hat{G}_m(Z_tj)
  ## - index : value i that tell us the time step tt[i] in wich we are doing the estimation
  ## - tt: vector tt of length(Z) containing the Z's time steps
  ## - t_eval: evaluation vector (in practice t_eval = tt)
  ## - bandwidth: kernel bandwith value
  ##
  ## Output :
  ## - a J x 2 matrix which first and second columns respectively contain the values
  ##      (  J * Kij_tj/sum(Kij_tj)*\hat{G}_m(Z_tj) - p_12t(lambda_tj,k_tj)  ) and
  ##        (  J * Kij_tj/sum(Kij_tj)*(\hat{G}_m(Z_tj)^2) - p_13t(lambda_tj,k_tj)  )
  ##
  ## Requires : laplaceWeibull()
  
  point_eval=t_eval.vec[index]
  lambdaval=theta[1]; kval=theta[2]
  
  p12val <- laplaceWeibull(j=1,lambda=lambdaval,k=kval,lowerbnd=10^(-5))
  p13val <- laplaceWeibull(j=2,lambda=lambdaval,k=kval,lowerbnd=10^(-5))
  
  Kij_ti <- outer(point_eval,tt.vec,function(zz,z)dEpan((zz-z)/bandwidth))
  Kij_ti <- t(Kij_ti)
  W <- Kij_ti/sum(Kij_ti)
  
  p12_nonMoyenne <- W*vecx
  p12_nonMoyenne <- as.vector (p12_nonMoyenne)
  p13_nonMoyenne <- W*(vecx)^2
  p13_nonMoyenne <- as.vector (p13_nonMoyenne)
  
  M <- cbind( length(vecx)*p12_nonMoyenne - p12val ,  length(vecx)*p13_nonMoyenne - p13val )
  return(M)
}

laplaceWeibull <- function(j,lambda,k,lowerbnd=10^(-5),upperbnd=1,fac=1,tol=10^(-5)){
  
  ## This function computes E(G(Z)^j)
  ## for any given integer j, where W is a Weibull(lambda,k) variable.
  ##
  ## The computation is based on integration rather than partial power series,
  ##
  ## Input :
  ## - single values of lambda and k (in practice issued from weibullGMM_NonStationaire
  ##   or other estimation method)
  ## - j value , in practice j = r-1
  ## - other optional parameters controll the way the numerical integration is conducted
  ##
  ##
  ## Output :
  ## - the numerical evaluation of E(exp(-j*W)) when W~Weibull(lambda,k)
  ##
  ## Requires : funcLaplace()
  ##
  vala=fac*(j*lambda)^k
  upperbndmodif=upperbnd^vala
  lowerbndmodif=upperbndmodif*10^(-5)
  I <- integrate(f=funcLaplace,
                 lower=lowerbndmodif,
                 upper=upperbndmodif,
                 subdivisions=1000L,
                 rel.tol=tol,
                 m=j,lam=lambda,k=k,a=vala,
                 stop.on.error = FALSE)
  resultat <- I$value
  return(resultat)
}

funcLaplace <- function(x,m,lam,k,a){
  # Utilitary function used in the integrate() statement in function laplaceWeibull()
  (1/a) * exp( -(m*lam/a^(1/k))*(-log(x))^(1/k) ) * x^(1/a - 1)
}

p12_NonPar <- function(X.mat,Z.mat,tt,t_eval,h,kern= dEpan){
  
  ## Computation of p12_t using kernels

  X.mat <- as.matrix(X.mat)
  Z.mat <- as.matrix(Z.mat)
  
  GmZ.mat <- matGZ_func(X.mat,Z.mat)
  N <- dim(GmZ.mat)[2]
  J <- dim(GmZ.mat)[1]
  Kij <- outer(t_eval,tt,function(zz,z) kern((zz - z) / h))
  W <- Kij / rowSums(Kij)
  p12.mat <- matrix(NA, ncol = N, nrow = length(t_eval))
  for(i in 1:N){
    GZ_2<-GmZ.mat[,i]
    p12.mat[,i] <- W %*% GZ_2
  }
  return( list("p12.mat"=p12.mat,"GmZ.mat"=GmZ.mat) )
}

p13_NonPar <- function(X.mat,Z.mat,tt,t_eval,h,kern= dEpan){
  
  ## Computation of p13_t using kernels

  X.mat <- as.matrix(X.mat)
  Z.mat <- as.matrix(Z.mat)
  
  GmZ.mat <- matGZ_func(X.mat,Z.mat)
  N <- dim(GmZ.mat)[2]
  J <- dim(GmZ.mat)[1]
  Kij <- outer(t_eval,tt,function(zz,z) kern((zz - z) / h))
  W <- Kij / rowSums(Kij)
  p13.mat <- matrix(NA, ncol = N, nrow = length(t_eval))
  for(i in 1:N){
    GZ_2<-(GmZ.mat[,i])^2
    p13.mat[,i] <- W %*% GZ_2
  }
  return(p13.mat)
}

p1rfarW_temps<- function(lam.mat,k.mat,r.mat,lowerbnd=10^(-5)){
  
  r.mat <- as.matrix(r.mat)
  lam.mat<-as.matrix(lam.mat)
  k.mat<-as.matrix(k.mat)
  n = dim(lam.mat)[1]
  N = dim(lam.mat)[2]
  p1r.mat = matrix(0,n,N)
  for (j in 1:N){
    lam.vec = lam.mat[,j]
    k.vec = k.mat[,j]
    r.vec = r.mat[,j]
    p1r.vec = p1r.mat[,j]
    for (i in 1:n){ 
      p1r_t <- laplaceWeibull(r.vec[i]-1,lam.vec[i],k.vec[i],lowerbnd=lowerbnd)
      p1r.vec[i] <- p1r_t
    }
    p1r.mat[,j] <- p1r.vec
  }
  far.mat <- 1 - 1 / (r.mat * p1r.mat)
  return( list("p1r"=p1r.mat,"far"=far.mat) )
}

matcovNA_A2_jit <- function(p12_hat_t,p13_hat_t){
  ## This fucntion computes the unweighted list A_2ij of NA's variance-covariance matrix
  ##    A2_ji = p13_hat.ti-(p12_hat.ti)^2
  list_unweighed_matcov_NA_aji <- p13_hat_t - (p12_hat_t)^2

  return (list_unweighed_matcov_NA_aji)
}

matcovNA_B2_jit <- function(p13_hat_t,p15W_t){
  ## This fucntion computes the unweighted list B2_ij of NA's variance-covariance matrix
  ##    B2_ji = p135W.ti - (p13_hat.ti)^2
  list_unweighed_matcov_NA_bji <- p15W_t -  (p13_hat_t)^2

  return (list_unweighed_matcov_NA_bji)
}

matcovNA_C2_jit <- function(p12_hat_t,p13_hat_t,p14W_t){
  ## This fucntion computes the unweighted list C2_ji of NA's variance-covariance matrix
  ##    C2_ji = p14W_ti - (p12_hat.ti * p13_hat.ti)
  ##
  list_unweighed_matcov_NA_cji <- p14W_t - (p12_hat_t * p13_hat_t)

  return (list_unweighed_matcov_NA_cji)
}

matcovNB_aji <- function(p12_hat_t,lambda_t,k_t,upperbnd=Inf,lowerbnd=10^(-5),tol=10^(-5),method="MC",m=2000){
  
  ## This function computes the unweighted term Aij of NB_n's asymptotic variance-covariance matrix
  ##
  ##    Aij = E( min ( G(Z_tj),G(Z_ti) ) - G(Z_tj )* G(Z_ti) ) =  E ( M_2ji + ( p12_hat.j * p12_hat.i ) )
  ##
  ## Input :
  ## - \hat{p}_12.t and \hat{p}_13.t vectors from kernel estimation
  ## - two vectors of length J containing Weibull's estimated parameters \hat{lambda}_t,\hat{k}_t
  ##
  ## Output :
  ## - Vector containing A1_ji terms
  ##
  ## Used in : matcovNB
  ##
  ## Requires : Mrfuncji()
  
  # cat("start: computation of A1ji","\n")
  method <- method
  
  list_unweighed_matcov_NB_aji <- rep(NA,31626)
  #list_unweighed_matcov_NB_aji_2 <- rep(NA,31626)
  inx <- 1
  for (j in 1:length(p12_hat_t)){
    # cat("start: components for j= ",j,"\n")
    lambda.j <- lambda_t[j]
    k.j <- k_t[j]
    p12_hat.j <- p12_hat_t[j]
    
    for (i in 1:j){
      # cat("start: components for i= ",i,"\n")
      lambda.i <- lambda_t[i]
      k.i <- k_t[i]
      p12_hat.i <- p12_hat_t[i]
      
      ##
      if (method == "integration"){
        matcov_NB_aji <-Mrfuncji(2,lambda.j,k.j,lambda.i,k.i,lowerbnd=lowerbnd,upperbnd=upperbnd,tol=tol) - ( p12_hat.j * p12_hat.i)
      }
      if (method == "CVMC"){
        matcov_NB_aji <- Mrji_CVMC(m = m,r=2,lambda.j,k.j,lambda.i,k.i) - ( p12_hat.j * p12_hat.i)
      }
      if (method == "MC"){
        matcov_NB_aji <- Mrji_MC(m = m,r=2,lambda.j,k.j,lambda.i,k.i) - ( p12_hat.j * p12_hat.i)
      }
      ##
      
      
      list_unweighed_matcov_NB_aji[inx] <- matcov_NB_aji
      #list_unweighed_matcov_NB_aji_2[inx] <- matcov_NB_aji_2
      inx <- inx +1
      #cat("done: unweighed list of Aji components for ",j,", ",i,"\n")
    }
  }
  return (list_unweighed_matcov_NB_aji)
}

Mrfuncji <- function(r,lam.j,k.j,lam.i,k.i,lowerbnd=10^(-5),upperbnd=Inf,tol=10^(-5)){
  
  I <- integrate(f=funcMrji,
                 lower=lowerbnd,
                 upper=upperbnd,
                 subdivisions=100L,
                 rel.tol=tol,
                 m=r-1,
                 lam.j=lam.j,
                 k.j=k.j,
                 lam.i=lam.i,
                 k.i=k.i,
                 lowerbnd=lowerbnd,
                 stop.on.error = FALSE)
  return(I$value)
}

funcMrji <- function(x,m,lam.j,k.j,lam.i,k.i,lowerbnd=10^(-5),tol=10^(-5)){
  # Utilitary function inside Mrfuncji
  nx=length(x)
  facteur1.j = exp( -(m) * x) * k.j / ((lam.j)^k.j) * x^(k.j-1) * exp( -(x/lam.j)^k.j )
  facteur1.i = exp( -(m) * x) * k.i / ((lam.i)^k.i) * x^(k.i-1) * exp( -(x/lam.i)^k.i )
  facteur2.j=rep(0,nx)
  facteur2.i=rep(0,nx)
  for (i in 1:nx){
    facteur2.j[i] = integral_Mrji_Eji(m=m,lam=lam.j,k=k.j,lowerbnd=(10^-6)*x[i],upperbnd=x[i],tol=tol)
  }
  for (i in 1:nx){
    facteur2.i[i] = integral_Mrji_Eji(m=m,lam=lam.i,k=k.i,lowerbnd=(10^-6)*x[i],upperbnd=x[i],tol=tol)
  }
  return( facteur1.j * facteur2.i + facteur1.i * facteur2.j)
}

integral_Mrji_Eji <- function(m,lam,k,lowerbnd=10^(-5),upperbnd=Inf,tol=10^(-5)){
  I <- integrate(f=integrand_Mrji_Eji,
                 lower=lowerbnd,
                 upper=upperbnd,
                 subdivisions=1000L,
                 rel.tol=tol,
                 m=m,lam=lam, k=k,
                 stop.on.error = FALSE)
  resultat <- I$value
  return(resultat)
}

integrand_Mrji_Eji <- function(u,m,lam,k,a){
  # m est r-1 ou r-2
  exp( -(m-1) * u) * k / ((lam)^k) * u^(k-1) * exp( -(u/lam)^k )
}

# Mrji by control variate Monte Carlo
Mrji_CVMC <- function(m =2000,r,lam.j,k.j,lam.i,k.i){
  # 3. Mrji by control variate Monte Carlo
  
  # B) Option 1: Our control variate Y  is
  # Y = exp( -(r-2) rWeibull.j - (r-1) rWeibull.i)  where i<j
  # i.e we remplace min( - rWeibull.j, - rWeibull.i) by - rWeibull.i
  # In the case r = 2: EY = p,2.i, in the case r = 3: EY = p,2.j * p1,3.i
  # We will first estimate p,2.i, p,2.j, p1,3.i using control variate MC
  # and then we use EY as control variate of a second CVMC
  
  # setting the parameters alpha.j, beta.j, alpha.i, beta.i of two
  # Gamma(alpha,beta) distributed r.v. such that its mean and variance
  # are equal to the one of the Weibull(k,l)
  wj.mean <- lam.j * gamma(1+1/k.j)
  wj.var <- lam.j *lam.j* (gamma(1+2/k.j) - (gamma(1+1/k.j))^2)
  beta.j <- wj.mean/wj.var
  alpha.j<- beta.j * wj.mean
  
  wi.mean <- lam.i * gamma(1+1/k.i)
  wi.var <- lam.i *lam.i* (gamma(1+2/k.i) - (gamma(1+1/k.i))^2)
  beta.i <- wi.mean/wi.var
  alpha.i<- beta.i * wi.mean
  
  # sampling  correlated  Weibull(k,l) and Gamma(alpha,beta) samples
  rUnif.i <- runif(m)
  rGamma.i <- qgamma(rUnif.i,shape = alpha.i,scale = 1/beta.i)
  rWeibull.i <-qweibull(rUnif.i,shape = k.i,scale = lam.i)
  
  rUnif.j <- runif(m)
  rGamma.j <- qgamma(rUnif.j,shape = alpha.j,scale = 1/beta.j)
  rWeibull.j <-qweibull(rUnif.j,shape = k.j,scale = lam.j)
  
  
  if(r==2){
    # FIRST STEP : control variate for the estimation of p_{1,r}
    # EY = p12.i
    EY.i<- (1+(2-1)/beta.i)^(-alpha.i)
    VY.i <- (1+2*(r-1)/beta.i)^(-alpha.i) - (EY.i)^2
    X.i <-rep(NA,m); Y.i <- X.i
    X.i <- exp(-(2-1)*rWeibull.i)
    Y.i <- exp(-(2-1)*rGamma.i)
    var.Y.i <- (sd(Y.i))^2
    cov.XY.i <- cov(X.i,Y.i)
    c.i <- cov.XY.i/var.Y.i
    EY.i<- c.i * EY.i
    VY.i <- c.i^2 * VY.i
    Y.i <- c.i * exp(-(2-1)*rGamma.i)
    mean.diff.i <- mean(X.i - Y.i)
    var.diff.i <- (sd(X.i-Y.i))^2
    p12.i <- mean.diff.i + EY.i
    
    EY <- p12.i
    
    # SECOND STEP : control variate for the estimation of Mrji
    X <- exp(-(r-2)*rWeibull.j - (r-2)*rWeibull.i - pmax(rWeibull.j,rWeibull.i) )
    Y <- exp(-(r-1)*rWeibull.i)
    var.Y<-(sd(Y))^2
    cov.XY<-cov(X,Y)
    c <- cov.XY/var.Y
    EY <- c * EY
    Y <- c * Y
    
    mean.diff<-mean(X-Y)
    var.diff<-(sd(X-Y))^2
    mean.Y<-mean(Y)
    var.Y<-(sd(Y))^2
    mean.X<-mean(X)
    var.X<-(sd(X))^2
    cov.XY<-cov(X,Y)
    
    Mr_CVMC<- mean.diff+EY
    
  }
  
  
  if(r==3){
    # FIRST STEP : control variate for the estimation of p_{1,r}
    # EY = p12.j * p13.i
    EY.i<- (1+(3-1)/beta.i)^(-alpha.i)
    VY.i <- (1+2*(3-1)/beta.i)^(-alpha.i) - (EY.i)^2
    X.i <-rep(NA,m); Y.i <- X.i
    X.i <- exp(-(3-1)*rWeibull.i)
    Y.i <- exp(-(3-1)*rGamma.i)
    var.Y.i <- (sd(Y.i))^2
    cov.XY.i <- cov(X.i,Y.i)
    c.i <- cov.XY.i/var.Y.i
    EY.i<- c.i * EY.i
    #VY.i <- c.i^2 * VY.i
    Y.i <- c.i * Y.i
    mean.diff.i <- mean(X.i - Y.i)
    var.diff.i <- (sd(X.i-Y.i))^2
    p13.i <- mean.diff.i + EY.i
    
    EY.j<- (1+(2-1)/beta.j)^(-alpha.j)
    VY.j <- (1+2*(2-1)/beta.j)^(-alpha.j) - (EY.j)^2
    X.j <-rep(NA,m); Y.j <- X.j
    X.j <- exp(-(2-1)*rWeibull.j)
    Y.j <- exp(-(2-1)*rGamma.j)
    var.Y.j <- (sd(Y.j))^2
    cov.XY.j <- cov(X.j,Y.j)
    c.j <- cov.XY.j/var.Y.j
    EY.j<- c.j * EY.j
    VY.j <- c.j^2 * VY.j
    Y.j <- c.j * Y.j
    mean.diff.j <- mean(X.j - Y.j)
    var.diff.j <- (sd(X.j-Y.j))^2
    p12.j <- mean.diff.j + EY.j
    
    EY <- p12.j * p13.i
    
    # SECOND STEP : control variate for the estimation of Mrji
    X <- exp(-(r-2)*rWeibull.j - (r-2)*rWeibull.i - pmax(rWeibull.j,rWeibull.i) )
    Y <- exp(-(r-2)*rWeibull.j-(r-1)*rWeibull.i)
    var.Y<-(sd(Y))^2
    cov.XY<-cov(X,Y)
    c <- cov.XY/var.Y
    EY <- c * EY
    Y <- c * Y
    
    mean.diff<-mean(X-Y)
    var.diff<-(sd(X-Y))^2
    mean.Y<-mean(Y)
    var.Y<-(sd(Y))^2
    mean.X<-mean(X)
    var.X<-(sd(X))^2
    cov.XY<-cov(X,Y)
    
    Mr_CVMC<- mean.diff+EY
  }
  
  return(Mr_CVMC)
}

# Mrji by Monte Carlo simulation
invW<- function(p, lambda,k){lambda * (-(log(1-p)))^(1/k)}

Mrji_MC <- function(m = 2000,r,lam.j,k.j,lam.i,k.i){
  Wj <- invW(runif(m), lambda = lam.j, k = k.j)
  eWj <- exp(-(r-2)*Wj)
  
  Wi <- invW(runif(m), lambda = lam.i, k = k.i)
  eWi<- exp(-(r-2)*Wi)
  
  min_eW <- pmin(exp(-Wj),exp(-Wi))
  
  Mrji_MC <- mean(eWj*eWi*min_eW)
  
  return(Mrji_MC)
}

matcovNB_bji <- function(p13_hat_t,lambda_t,k_t,lowerbnd=10^(-5),tol=10^(-5),method="MC",m=2000){
  
  ## This fucntion computes the unweighted term Cij of NB_n's asymptotic variance-covariance matrix
  ##
  ##    Bij = 4* E( (G(Z_tj) * G(Z_ti) * [ min( G(Z_tj),G(Z_ti) ) - G(Z_tj) * G(Z_ti) ] ) = 4( M_3ji - p13_hat.j * p13_hat.i )
  ##
  ## Input :
  ## - \hat{p}_12.t and \hat{p}_13.t vectors from kernel estimation
  ## - two vectors of length t containing Weibull's estimated parameters \hat{lambda}_t,\hat{k}_t
  ##
  ## Output :
  ## - Vector containing B1_ji terms
  ##
  ## Used in : matcovNB
  ##
  ## Requires : Mrfuncji()
  
  # cat("start: computation of B1ji","\n")
  method <- method
  
  list_unweighed_matcov_NB_bji <- rep(NA,31626)
  inx <- 1
  for (j in 1:length(p13_hat_t)){
    # cat("start: components for j= ",j,"\n")
    lambda.j <- lambda_t[j]
    k.j <- k_t[j]
    p13_hat.j <- p13_hat_t[j]
    
    for (i in 1:j){
      # cat("start: components for i= ",i,"\n")
      lambda.i <- lambda_t[i]
      k.i <- k_t[i]
      p13_hat.i <- p13_hat_t[i]
      
      ##
      if (method == "integration"){
        NB_bji <- 4*(Mrfuncji(3,lambda.j,k.j,lambda.i,k.i,lowerbnd=lowerbnd,upperbnd=Inf,tol=tol) - (p13_hat.j * p13_hat.i))
      }
      if (method == "CVMC"){
        NB_bji <- 4*(Mrji_CVMC(m = m,r=3,lambda.j,k.j,lambda.i,k.i) - (p13_hat.j * p13_hat.i))
      }
      if (method == "MC"){
        NB_bji <- 4*(Mrji_MC(m = m,r=3,lambda.j,k.j,lambda.i,k.i) - (p13_hat.j * p13_hat.i))
      }
      
      list_unweighed_matcov_NB_bji[inx] <- NB_bji
      inx <- inx +1
    }
  }
  return (list_unweighed_matcov_NB_bji)
}

matcovNB_cji <- function(p12_hat_t,p13_hat_t,lambda_t,k_t,lowerbnd=10^(-5),tol=10^(-5),method="MC",m=2000){
  
  ## This fucntion computes the unweighted term Cij of NB_n's asymptotic variance-covariance matrix
  ##
  ##    Cji =  E ( 2* G(Z_tj) * [min(G(Z_tj),G(Z_ti)) - G(Z_tj)*G(Z_ti)] )  =  2* ( E_ji - p13_hat.j  * p12_hat.i )
  ##
  ## Input :
  ## - \hat{p}_12.t and \hat{p}_13.t vectors from kernel estimation
  ## - two vectors of length J containing Weibull's estimated parameters \hat{lambda}_t,\hat{k}_t
  ##
  ## Output :
  ## - Vector containing C1_ji terms
  ##
  ## Used in : matcovNB
  ##
  ## Requires : calculEGzjminGzjGzi_partieA1() calculEGzjminGzjGzi_partieB()
  
  # cat("start: computation of C1ji","\n")
  method <- method
  
  list_unweighed_matcov_NB_cji <- rep(NA,31626)
  inx <- 1
  for (j in 1:length(p12_hat_t)){
    # cat("start: components for j= ",j,"\n")
    lambda.j <- lambda_t[j]
    k.j <- k_t[j]
    p13_hat.j <- p13_hat_t[j]
    p12_hat.j <- p12_hat_t[j]
    
    for (i in 1:j){
      # cat("start: components for i= ",i,"\n")
      lambda.i <- lambda_t[i]
      k.i <- k_t[i]
      p12_hat.i <- p12_hat_t[i]
      p13_hat.i <- p13_hat_t[i]
      
      # Adding all parts
      
      ##
      if (method == "integration"){
        Eji <- Eji_func(lambda.j,k.j,lambda.i,k.i,lowerbnd=10^(-5),upperbnd=Inf,tol=10^(-5))
        Eij <- Eji_func(lambda.i,k.i,lambda.j,k.j,lowerbnd=10^(-5),upperbnd=Inf,tol=10^(-5))
      }
      if (method == "CVMC"){
        Eji <- Eji_CVMC(m = m,lambda.j,k.j,lambda.i,k.i)
        Eij <- Eji_CVMC(m = m,lambda.i,k.i,lambda.j,k.j)
      }
      if (method == "MC"){
        Eji <- Eji_MC(m = m,lambda.j,k.j,lambda.i,k.i)
        Eij <- Eji_MC(m = m,lambda.i,k.i,lambda.j,k.j)
      }
      
      matcov_NB_cji <- 2*(Eji+Eij) - 2*((p13_hat.j * p12_hat.i)+(p13_hat.i * p12_hat.j))
      
      list_unweighed_matcov_NB_cji[inx] <- matcov_NB_cji
      inx <- inx +1
      
    }
  }
  return (list_unweighed_matcov_NB_cji)
}

Eji_func <- function(lam.j,k.j,lam.i,k.i,lowerbnd=10^(-5),upperbnd=Inf,tol=10^(-5)){
  
  I <- integrate(f=funcEji,
                 lower=lowerbnd,
                 upper=upperbnd,
                 subdivisions=100L,
                 rel.tol=tol,
                 lam.j=lam.j,
                 k.j=k.j,
                 lam.i=lam.i,
                 k.i=k.i,
                 lowerbnd=lowerbnd,
                 stop.on.error = FALSE)
  return(I$value)
}

funcEji <- function(x,lam.j,k.j,lam.i,k.i,lowerbnd=10^(-5),tol=10^(-5)){
  # Utilitary function inside Mrfuncji
  nx=length(x)
  facteur1.j = exp( -(2) * x) * k.j / ((lam.j)^k.j) * x^(k.j-1) * exp( -(x/lam.j)^k.j )
  facteur1.i = exp( - x) * exp( -(x/lam.i)^k.i ) * k.i / ((lam.i)^k.i) * x^(k.i-1)
  facteur2.j=rep(0,nx)
  facteur2.i=rep(0,nx)
  for (i in 1:nx){
    facteur2.j[i] = integral_Mrji_Eji(m=2,lam=lam.j,k=k.j,lowerbnd=(10^-6)*x[i],upperbnd=x[i],tol=tol)
  }
  for (i in 1:nx){
    facteur2.i[i] = integral_Mrji_Eji(m=1,lam=lam.i,k=k.i,lowerbnd=(10^-6)*x[i],upperbnd=x[i],tol=tol)
  }
  return( facteur1.i * facteur2.j + facteur1.j * facteur2.i)
}

integral_Mrji_Eji <- function(m,lam,k,lowerbnd=10^(-5),upperbnd=Inf,tol=10^(-5)){
  I <- integrate(f=integrand_Mrji_Eji,
                 lower=lowerbnd,
                 upper=upperbnd,
                 subdivisions=1000L,
                 rel.tol=tol,
                 m=m,lam=lam, k=k,
                 stop.on.error = FALSE)
  resultat <- I$value
  return(resultat)
}

integrand_Mrji_Eji <- function(u,m,lam,k,a){
  # m est r-1 ou r-2
  exp( -(m-1) * u) * k / ((lam)^k) * u^(k-1) * exp( -(u/lam)^k )
}

Eji_CVMC <- function(m = 2000,lam.j,k.j,lam.i,k.i){
  
  # setting the parameters alpha.j, beta.j, alpha.i, beta.i of two
  # Gamma(alpha,beta) distributed r.v. such that its mean and variance
  # are equal to the one of the Weibull(k,l)
  wj.mean <- lam.j * gamma(1+1/k.j)
  wj.var <- lam.j *lam.j* (gamma(1+2/k.j) - (gamma(1+1/k.j))^2)
  beta.j <- wj.mean/wj.var
  alpha.j<- beta.j * wj.mean
  
  wi.mean <- lam.i * gamma(1+1/k.i)
  wi.var <- lam.i *lam.i* (gamma(1+2/k.i) - (gamma(1+1/k.i))^2)
  beta.i <- wi.mean/wi.var
  alpha.i<- beta.i * wi.mean
  
  # sampling  correlated  Weibull(k,l) and Gamma(alpha,beta) samples
  rUnif.i <- runif(m)
  rGamma.i <- qgamma(rUnif.i,shape = alpha.i,scale = 1/beta.i)
  rWeibull.i <-qweibull(rUnif.i,shape = k.i,scale = lam.i)
  
  rUnif.j <- runif(m)
  rGamma.j <- qgamma(rUnif.j,shape = alpha.j,scale = 1/beta.j)
  rWeibull.j <-qweibull(rUnif.j,shape = k.j,scale = lam.j)
  
  # FIRST STEP : control variate for the estimation of p_{1,r}
  # EY = p12.j * p12.i
  EY.i<- (1+(2-1)/beta.i)^(-alpha.i)
  VY.i <- (1+2*(2-1)/beta.i)^(-alpha.i) - (EY.i)^2
  X.i <-rep(NA,m); Y.i <- X.i
  X.i <- exp(-(2-1)*rWeibull.i)
  Y.i <- exp(-(2-1)*rGamma.i)
  var.Y.i <- (sd(Y.i))^2
  cov.XY.i <- cov(X.i,Y.i)
  c.i <- cov.XY.i/var.Y.i
  EY.i<- c.i * EY.i
  VY.i <- c.i^2 * VY.i
  Y.i <- c.i * Y.i
  mean.diff.i <- mean(X.i - Y.i)
  var.diff.i <- (sd(X.i-Y.i))^2
  p12.i <- mean.diff.i + EY.i
  
  EY.j<- (1+(2-1)/beta.j)^(-alpha.j)
  VY.j <- (1+2*(2-1)/beta.j)^(-alpha.j) - (EY.j)^2
  X.j <-rep(NA,m); Y.j <- X.j
  X.j <- exp(-(2-1)*rWeibull.j)
  Y.j <- exp(-(2-1)*rGamma.j)
  var.Y.j <- (sd(Y.j))^2
  cov.XY.j <- cov(X.j,Y.j)
  c.j <- cov.XY.i/var.Y.i
  EY.j<- c.j * EY.j
  VY.j <- c.j^2 * VY.j
  Y.j <- c.j * Y.j
  mean.diff.j <- mean(X.j - Y.j)
  var.diff.j <- (sd(X.j-Y.j))^2
  p12.j <- mean.diff.j + EY.j
  
  EY <- p12.j * p12.i
  
  # SECOND STEP : control variate for the estimation of Mrji
  X <- exp(-rWeibull.j - pmax(rWeibull.j,rWeibull.i) )
  Y <- exp(-rWeibull.j-rWeibull.i)
  var.Y<-(sd(Y))^2
  cov.XY<-cov(X,Y)
  c <- cov.XY/var.Y
  EY <- c * EY
  Y <- c * Y
  
  mean.diff<-mean(X-Y)
  var.diff<-(sd(X-Y))^2
  mean.Y<-mean(Y)
  var.Y<-(sd(Y))^2
  mean.X<-mean(X)
  var.X<-(sd(X))^2
  cov.XY<-cov(X,Y)
  
  Eji_CVMC<- mean.diff+EY
  
  
  
  return(Eji_CVMC)
}

invW<- function(p, lambda,k){lambda * (-(log(1-p)))^(1/k)}

Eji_MC <- function(m = 2000,lam.j,k.j,lam.i,k.i){
  Wj <- invW(runif(m), lambda = lam.j, k = k.j)
  eWj <- exp(-Wj)
  
  Wi <- invW(runif(m), lambda = lam.i, k = k.i)
  eWi <- exp(-Wi)
  
  min_eW <- pmin(eWj,eWi)
  
  Eji_MC <- mean(eWj * min_eW)
  
  return(Eji_MC)
}

jacobianFunctiong12 <- function(lam.vec,k.vec,debugg=FALSE){
  
  ## This function computes the J (2x2) jacobian matrices of function
  ## g : (lambda_tj,k_tj) ->
  ##     ( g1(lambda_tj,k_tj), g2(lambda,k) ) = ( E(G(Z_tj)) , E(G^2(Z_tj)) )
  ##
  ## at the values (lambda_t,k_t) given as inputs.
  ##
  ## Input :
  ## - vectors of same sizes containing values of lambda and k
  ##   (in practice, these values are estimates of lambda and k issued from weibullGMMestim())
  ##
  ## Output :
  ## - a list containing the J 2x2 jacobian matrices described above
  ##
  ## Used : weibullFsolve() (as J12), matcovtheta_t(), varp1rfar_t()
  ##
  ## Requires : dgjoverdlambdafunc(), dgjoverdkfunc()
  lam.vec <- as.numeric(lam.vec)
  k.vec <- as.numeric(k.vec)
  lvec <- length(lam.vec)
  listejacobiennes <- list()
  for (i in 1:lvec){
    lambda <- lam.vec[i] ; k <- k.vec[i]
    dg1surdlambda <- dgjoverdlambdafunc(1,lambda,k)
    if (debugg){ cat("i=",i,": dg1dlam,") }
    dg1surdk <- dgjoverdkfunc(1,lambda,k)
    if (debugg){ cat("dg1dk,")}
    dg2surdlambda <- dgjoverdlambdafunc(2,lambda,k)
    if (debugg){ cat("dg2dlam,")}
    dg2surdk <- dgjoverdkfunc(2,lambda,k)
    if (debugg){ cat("dg2dk\n")}
    
    listejacobiennes[[i]] <- matrix(c(dg1surdlambda,dg1surdk,dg2surdlambda,dg2surdk),
                                    2,2,byrow=TRUE)
  }
  return(list(listejacobiennes))
}

foncdgjoverdlambda <- function(u,j,lam,k,a){
  ## Utilitary function used inside dgjoverdlambdafunc()
  (-j/a^((1/k)+1)) * (-log(u))^(1/k) * exp( -(j*lam/a^(1/k))*(-log(u))^(1/k) ) * u^(1/a - 1)
}

dgjoverdlambdafunc <- function(j,lambda,k,lowerbnd=10^(-8),fac=0.5){
  
  ## This utilitary function computes the partial derivative lambda
  ## of the expectation E( exp(-j*W) ) where W is Weibull(lambda,k).
  ##
  ## Used in : jacobianFunctiong12(), jacobianFunctiongrminus1()
  ##
  ## Requires : foncdgjoverdlambda()
  
  vala=fac*(j*lambda)^k
  I <- integrate(f=foncdgjoverdlambda,lower=lowerbnd,upper=1,subdivisions=1000L,
                 j=j,lam=lambda,k=k,a=vala,
                 stop.on.error = FALSE)
  return(I$value)
}

foncdgjoverdk <- function(u,j,lam,k,a){
  ## Utilitary function used inside dgjoverdkfunc()
  (-lam/k^2) * log( (1/a)*(-log(u)) ) * foncdgjoverdlambda(u,j,lam,k,a)
}

dgjoverdkfunc <- function(j,lambda,k,lowerbnd=10^(-6),fac=1){
  
  ## This utilitary function computes the partial derivative k
  ## of the expectation E( exp(-j*W) ) where W is Weibull(lambda,k).
  ##
  ## Used in : jacobianFunctiong12(), jacobianFunctiongrminus1()
  ##
  ## Requires : foncdgjoverdk()
  
  vala=fac*(j*lambda)^k
  I <- integrate(f=foncdgjoverdk,lower=lowerbnd,upper=1,subdivisions=1000L,
                 j=j,lam=lambda,k=k,a=vala,stop.on.error = FALSE)
  return(I$value)
}

jacobianFunctiongrminus1 <- function(lam.vec,k.vec,r.vec){
  lam.vec <- as.numeric(lam.vec)
  k.vec <- as.numeric(k.vec)
  lvec <- length(lam.vec)
  listejacobiennes <- list()
  for (i in 1:lvec){
    lambda <- lam.vec[i] ; k <- k.vec[i] ; r <- r.vec[i]
    dgoverdlambda <- dgjoverdlambdafunc(r-1,lambda,k)
    dgoverdk <- dgjoverdkfunc(r-1,lambda,k)
    listejacobiennes[[i]] <- matrix(c(dgoverdlambda,dgoverdk),nrow=1,ncol=2)
  }
  return(listejacobiennes)
}

