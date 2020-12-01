VAR = function (x, p = 1, include.mean = T, fixed = NULL) {
   Tn = dim(x)[1]
   k = dim(x)[2]
   idm = k * p
   ne = Tn - p
   ist = p + 1
   y = x[ist:Tn, ]
   if (include.mean) {
      idm = idm + 1
      xmtx = cbind(rep(1, ne), x[p:(Tn - 1), ])
   } else {
      xmtx = x[p:(Tn - 1), ]
   }
   if (p > 1) {
      for (i in 2:p) {
         xmtx = cbind(xmtx, x[(ist - i):(Tn - i), ])
      }
   }
   ndim = ncol(xmtx)
   if (length(fixed) == 0) {
      paridx = matrix(1, ndim, k)
   } else {
      paridx = fixed
   }
   res = NULL
   beta = matrix(0, ndim, k)
   sdbeta = matrix(0, ndim, k)
   npar = 0
   for (i in 1:k) {
      idx = c(1:ndim)[paridx[, i] == 1]
      resi = y[, i]
      if (length(idx) > 0) {
         xm = as.matrix(xmtx[, idx])
         npar = npar + dim(xm)[2]
         xpx = t(xm) %*% xm
         xpxinv = solve(xpx)
         xpy = t(xm) %*% as.matrix(y[, i], ne, 1)
         betai = xpxinv %*% xpy
         beta[idx, i] = betai
         resi = y[, i] - xm %*% betai
         nee = dim(xm)[2]
         sse = sum(resi * resi)/(Tn - p - nee)
         dd = diag(xpxinv)
         sdbeta[idx, i] = sqrt(dd * sse)
      }
      res = cbind(res, resi)
   }
   sse = t(res) %*% res/(Tn - p)
   Phi = NULL
   Ph0 = NULL
   jst = 0
   if (include.mean) {
      Ph0 = beta[1, ]
      se = sdbeta[1, ]
      jst = 1
   }
   if (include.mean) {
      for (i in 1:k) {
         if (abs(Ph0[i]) > 1e-08) 
            npar = npar - 1
      }
   }
   for (i in 1:p) {
      phi = t(beta[(jst + 1):(jst + k), ])
      se = t(sdbeta[(jst + 1):(jst + k), ])
      jst = jst + k
      Phi = cbind(Phi, phi)
   }
   VAR <- list(data = x, residuals = res, secoef = sdbeta, 
               Sigma = sse, Phi = Phi, Ph0 = Ph0)
}
extract = function(data,k) {
   t = nrow(data)
   n = ncol(data)
   xx = crossprod(data);
   ee = eigen(xx)
   eval = ee$values
   evc = ee$vectors
   
   lam = sqrt(n)*evc[,1:k];
   fac = data%*%lam/n;
   
   return = list(lam=lam,fac=fac)
}

### BAYES_DFM
library("MASS")
library("matlab")
library("mvtnorm")
library("R.matlab")
options(warn=-1)

dfmdgp = function(L=NULL,Sigma=NULL,T=NULL,N=NULL,K=NULL,lag_f=NULL,PHI=NULL,PSI=NULL) {
   if (is.null(T)) {
      T=200
   }
   if (is.null(N)) {
      N=9
   }
   if (is.null(K)) {
      K=3
   }
   if (is.null(lag_f)) {
      lag_f=1;
   }
   if (is.null(L)) {
      L = matrix(c(1,0,0,0,1,0,0,0,1,0.99,0.60,0.34,
                0.99,0,0,0.78,0,0.54,
                0.35,0.85,0.78,0,0.33,0.90,
                0.78,0,0.90), ncol=3)
   }
   if (is.null(Sigma)) {
      Sigma = diag(c(0.2,0.19,0.36,0.2,0.2,0.19,0.19,0.36,0.36));
   }
   if (is.null(PHI)) {
      PHI = matrix(c(0.5,0,0,0,0.5,0,0,0,0.5),ncol=3);
   }
   if (is.null(PSI)) {
      PSI = matrix(c(1,0.5,0.5,0,1,0.5,0,0,1),ncol=3);
   }
   
   #BFMDGP - Bayesian Factor Model DGP
   Q = solve(PSI%*%t(PSI));
   
   f = rbind(matrix(rnorm(lag_f*K),lag_f,K), zeros(T,K));
   # Now generate f from VAR (L,PHI,PSI)
   for (nn in lag_f:(T+lag_f)){
      u = t(chol(Q))%*%rnorm(K);
      flag = embed(f,(lag_f));
      f[nn,] = flag[nn,]%*%PHI + t(u);
   }
   
   f = f[(lag_f+1):(T+lag_f),];
   X = f%*%t(L) + matrix(rnorm(T*N),ncol=N)%*%t(chol(Sigma))
   X
}

#BAYES_DFM This implements the state-space estimation of the dynamic
#factors, their loadings and autoregressive parameters. The program was
#written by Piotr Eliasz (2003) Princeton University for his PhD thesis and
#and used in Bernake, Boivin and Eliasz (2005), published in QJE.
#Modified by Dimitris Korobilis on Monday, June 11, 2007
#--------------------------------------------------------------------------
# Estimate dynamic factor model:
#                X_t = Lam x F_t + e_t
#                F_t = B(L) x F_t-1 + u_t
# The previous (Eliasz, 2005) version was:
#         [X_t ; Y_t] = [Lf Ly ; 0 I] x [F_t ; Y_t] + e_t
#         [F_t ; Y_t] = [B(L) ; 0 I] x [F_t-1 ; Y_t-1] + u_t
#
#--------------------------------------------------------------------------
      
#---------------------------LOAD DATA--------------------------------------
# simulated dataset

#X = dfmdgp(L,Sigma,T,N,K,lag_f,PHI,PSI)
X = as.matrix(read.table("/Users/user/Dropbox/18. Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/R/Yraw.dat"))
T=dim(X)[1]
N=dim(X)[2]

#Demean xraw
for (i in 1:N) {
   X[,i] = X[,i]-mean(X[,i])
}

# Number of factors & lags in B(L):
K=2;
lags=2;
#----------------------------PRELIMINARIES---------------------------------
# Set some Gibbs - related preliminaries
nods = 500;   # Number of draws
bid = 200;    # Number of burn-in-draws
thin = 1;      # Consider every thin-th draw (thin value)
    
# store draws in:
Ldraw=zeros(nods-bid,N,K);
Bdraw=zeros(nods-bid,K,K,lags);
Qdraw=zeros(nods-bid,K,K);
Fdraw=zeros(nods-bid,T,K);
Rdraw=zeros(nods-bid,N);
#********************************************************

# STANDARDIZE for PC only
X_st = X
for (i in 1:N) {
   X_st[,i] = X[,i]/sd(X[,i])
}
# First step - extract PC from X
data = X_st
EX=psych::principal(X_st,nfactors=K)
F0=EX$scores
Lf=EX$loadings
Lf=Lf[1:ncol(X_st),1:K]

#EX=extract(X_st,K)
#F0=EX$fac
#Lf=EX$lam

# Transform factors and loadings for LE normalization
Q=qr(t(Lf))
ql=qr.Q(Q)
rl=qr.R(Q)
Lf=rl;  # do not transpose yet, is upper triangular
F0=F0%*%ql;

# Need identity in the first K columns of Lf, call them A for now
A=Lf[,1:K];
if (ncol(Lf)==K){
  Lf=diag(K)
} else {
  Lf=t(cbind(diag(K),solve(A)%*%Lf[,(K+1):N]));
}
F0=F0%*%A;

# Obtain R:
e=X_st-F0%*%t(Lf);
R=(t(e)%*%e)/T;
R=diag(diag(R));
L=Lf;
# Run a VAR in F, obtain initial B and Q

ydata=F0
lag=1
xdata=NULL

var0 = VAR(F0,lags)
B = var0$Phi
v = var0$residuals
Q = var0$Sigma

# Put it all in state-space representation, write obs equ as XY=FY*L+e
XY=X;   #Tx(N+M)
FY=F0;

# adjust for lags in state equation, Q is KxK
Q=rbind(cbind(Q, matrix(0,K,K*(lags-1))),matrix(0,K*(lags-1),K*lags));
B=rbind(B,cbind(diag(K*(lags-1)),matrix(0,K*(lags-1),K)));

# start with
Sp=matrix(0,K*lags,1);
Pp=diag(K*lags);  
               
# Proper Priors:-----------------------------------------------------------
# on VAR -- Normal-Wishart, after Kadiyala, Karlsson, 1997
# on Q -- si
# on B -- increasing tightness
# on observable equation:
# N(0,I)-iG(3,0.001)
               
# prior distributions for VAR part, need B and Q
vo=K+2;
s0=3;
alpha=0.001;
L_var_prior=diag(K);
Qi=zeros(K,1);
               
# singles out latent factors
indexnM=ones(K,lags);
indexnM=as.matrix(which(indexnM==1));
#***************End of Preliminaries & PriorSpecification******************
                 
#==========================================================================
#========================== Start Sampling ================================
#==========================================================================

#************************** Start the Gibbs "loop" ************************
print('Number of iterations');
for (rep in 1:nods) {
  if (rep%%200==0) {
    print(paste0(round(rep/nods*100,2),"%"));
  }
               
  # STEP 1. =========|DRAW FACTORS
  # generate Gibbs draws of the factors
  H=L;
  F=B;
  t=nrow(XY)
  n=ncol(XY)
  kml=nrow(Sp);
  km=ncol(L);
  S=zeros(t,kml);
  P=zeros(kml^2,t);
  Sdraw=zeros(t,kml);
  for (i in 1:t) {
    y = XY[i,];
    nu = y - H%*%Sp[1:km];   # conditional forecast error
    f = H%*%Pp[1:km,1:km]%*%t(H) + R;    # variance of the conditional forecast error
    finv = solve(f);
               
    Stt = Sp + Pp[,1:km]%*%t(H)%*%finv%*%nu;
    Ptt = Pp - Pp[,1:km]%*%t(H)%*%finv%*%H%*%Pp[1:km,];
               
    if (i < t) {
      Sp = F%*%Stt;
      Pp = F%*%Ptt%*%t(F) + Q;
    }
               
    S[i,] = t(Stt);
    P[,i] = matrix(Ptt,kml^2,1);
  }
               
  # draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
  Sdraw[t,]=S[t,];
  Sdraw[t,indexnM]=rmvnorm(1,mean=c(Sdraw[t,]),sigma=Ptt);

  # iterate 'down', drawing at each step, use modification for singular Q
  Qstar=Q[1:km,1:km];
  Fstar=F[1:km,];
  for (i in (t-1):1) {
    Sf = matrix(Sdraw[(i+1),1:km],ncol=1);
    Stt = S[i,];
    Ptt = matrix(P[,i],kml,kml);
    f = Fstar%*%Ptt%*%t(Fstar) + Qstar;
    finv = solve(f);
    nu = Sf - Fstar%*%Stt;
                                       
    Smean = Stt + Ptt%*%t(Fstar)%*%finv%*%nu;
    Svar = Ptt - Ptt%*%t(Fstar)%*%finv%*%Fstar%*%Ptt;
                                       
    Sdraw[i,] = Smean;
    Sdraw[i,indexnM] = rmvnorm(1,t(Sdraw[i,indexnM]),Svar[indexnM,indexnM]);
  }
  FY=Sdraw[,1:km];
  # Demean
  FY=FY-matrix(rep(colMeans(FY),T),T,byrow=T);
  
  # STEP 2. ========|DRAW COEFFICIENTS
  # -----------------------2.1. STATE EQUATION---------------------------
  # first univ AR for scale in priors
  for (i in 1:km) {
    xxx = embed(FY[,i],lags+1)
    lm1 = lm(xxx[,1]~xxx[,-1])
    Bi = lm1$coefficients[-1]
    vi = lm1$residuals
    Qi[i]=sum(vi^2)/length(vi-lags-1)
  }
  Q_prior=diag(c(Qi));
  B_var_prior=diag(c(kronecker(1/t(Qi),1/(1:lags))));
  var3=VAR(FY,lags);
  Bd=var3$Phi
  v=var3$residuals
  Qd=var3$Sigma
  B_hat=t(Bd);
  Z=zeros(c(T,km,lags));
  for (i in 1:lags){
    Z[(lags+1):T,,i]=FY[(lags+1-i):(T-i),];
  }
  Z = matrix(Z,ncol=lags*2)
  Z = Z[(lags+1):T,];
  iB_var_prior=solve(B_var_prior);
  B_var_post=solve(iB_var_prior+t(Z)%*%Z);
  B_post=B_var_post%*%crossprod(Z)%*%B_hat;
  Q_post=t(B_hat)%*%crossprod(Z)%*%B_hat+Q_prior+(T-lags)*Qd-t(B_post)%*%(iB_var_prior+crossprod(Z))%*%B_post;
                                
  # Draw Q from inverse Wishart
  iQd=matrix(rnorm((T+vo)*km),(T+vo),km)%*%chol(ginv(Q_post));
  Qd=solve(t(iQd)%*%iQd);
  Q[1:km,1:km]=Qd;
                                       
  # Draw B conditional on Q
  vecB_post=matrix(B_post,km*km*lags,1);
  vecBd = vecB_post+t(chol(kronecker(Qd,B_var_post)))%*%rnorm(km*km*lags);
  Bd = t(matrix(vecBd,km*lags,km));
  B[1:km,]=Bd;
                                       
  # Truncate to ensure stationarity
  while (max(abs(eigen(B)$values))>0.999) {
    vecBd = vecB_post+t(chol(kronecker(Qd,B_var_post)))%*%matrix(rnorm(km*km*lags),ncol=1);
    Bd = t(matrix(vecBd,km*lags,km))
    B[1:km,]=Bd;
  }
                                       
  # ----------------------2.2. OBSERVATION EQUATION----------------------
  # OLS quantities
  L_OLS = solve(crossprod(FY))%*%(t(FY)%*%X[,c((K+1):N)]);
  R_OLS = t(X - FY%*%t(L))%*%(X - FY%*%t(L))/(T-N);
  
  if (K==ncol(X)) {
    L=diag(K)
  } else {
    L=t(cbind(diag(K), L_OLS));
  }             
  
  for (n in 1:N) {
    ed=matrix(X[,n],ncol=1)-FY%*%t(matrix(L[n,],nrow=1));
               
    # draw R(n,n)
    R_bar=s0+t(ed)%*%ed+matrix(L[n,],nrow=1)%*%ginv(L_var_prior+ginv(t(FY)%*%FY))%*%matrix(t(L[n,]),ncol=1);        
    Rd=rchisq(1,T+alpha);
    Rd=R_bar/Rd;
    R[n,n]=Rd;
  }

  # Save draws
  if (rep > bid) {
    Ldraw[(rep-bid),,]=L[1:N,1:km];
    Bdraw[(rep-bid),,,]=array(Bd,c(km,km,lags));
    Qdraw[(rep-bid),,]=Qd;
    Fdraw[(rep-bid),,]=FY;
    Rdraw[(rep-bid),]=diag(R);
  }
}

# =========================================================================
# ==========================Finished Sampling==============================
# =========================================================================

# Do thining in case of high correlation
thin_val = seq(1,((nods-bid)/thin),thin);
Ldraw = Ldraw[thin_val,,];
Bdraw = Bdraw[thin_val,,,];
Qdraw = Qdraw[thin_val,,];
Fdraw = Fdraw[thin_val,,];
Rdraw = Rdraw[thin_val,];

# Average over Gibbs draws
Fdraw2=apply(Fdraw,2:3,mean);
Fdraw2
Ldraw2=apply(Ldraw,2:3,mean);
Ldraw2
Qdraw2=apply(Qdraw,2:3,mean);
Qdraw2
Bdraw2=apply(Bdraw,2:4,mean);
Bdraw2
Rdraw2=apply(Rdraw,2,mean);
Rdraw2

# Get matrix of autoregressive parameters B
Betas = NULL;
for (dd in 1:lags) {
   B = apply(Bdraw,2:4,mean)
   beta = B[,,dd];
   beta_new = zeros(K,K);
   for (jj in 1:K) {
      beta_new[,jj] = beta[,jj];
   }
   Betas = cbind(Betas, beta_new); #ok<AGROW>
}
Betas

par(mfcol = c(K,1), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:K) {
   plot(Fdraw2[,i],type="l",xaxs="i",las=1,xlab="",ylab="",tck=.02)
}

### END
