ts_prior=function(rawdat,tau,p,plag) {
   yt = t(rawdat[(plag+1):(tau+plag),]);
   # m is the number of elements in the state vector
   m = p + plag*(p^2);
   Zt=NULL;
   
   for (i in (plag+1):(tau+plag)) {
      ztemp = diag(p);
      for (j in 1:plag) {
         xlag = rawdat[(i-j),1:p];
         xtemp = matrix(0,p,p*p);
         for (jj in 1:p) {
            xtemp[jj,((jj-1)*p+1):(jj*p)] = xlag;
         }
         ztemp = cbind(ztemp, xtemp);
      }
      Zt = rbind(Zt, ztemp);
   }
   
   vbar = matrix(0,m,m);
   xhy = matrix(0,m,1);
   
   for (i in 1:tau) {
      zhat1 = Zt[((i-1)*p+1):(i*p),];
      vbar = vbar + t(zhat1)%*%zhat1;
      xhy = xhy + t(zhat1)%*%as.matrix(yt[,i]);
   }
   vbar = solve(vbar);
   aols = vbar%*%xhy;
   
   sse2 = matrix(0,p,p);
   for (i in 1:tau) {
      zhat1 = Zt[((i-1)*p+1):(i*p),];
      sse2 = sse2 + (yt[,i] - zhat1%*%aols)%*%t(yt[,i] - zhat1%*%aols);
   }
   hbar = sse2/tau;
   vbar = matrix(0,m,m);
   for (i in 1:tau) {
      zhat1 = Zt[((i-1)*p+1):(i*p),];
      vbar = vbar + t(zhat1)%*%solve(hbar)%*%zhat1;
   }
   vbar = solve(vbar);
   achol = t(chol(hbar));
   ssig = matrix(0,p,p);
   for (i in 1:p) {
      ssig[i,i] = achol[i,i]; 
      for (j in 1:p) {
         achol[j,i] = achol[j,i]/ssig[i,i];
      }
   }
   achol = solve(achol);
   numa = p*(p-1)/2;
   a0 = matrix(0,numa,1);
   ic = 1;
   for (i in 2:p) {
      for (j in 1:(i-1)) {
         a0[ic,1] = achol[i,j];
         ic = ic + 1;
      }
   }
   ssig1 = matrix(0,p,1);
   for (i in 1:p) {
      ssig1[i,1] = log(ssig[i,i]^2);
   }
   
   hbar1 = solve(tau*hbar);
   hdraw = matrix(0,p,p);
   a02mo = matrix(0,numa,numa);
   a0mean = matrix(0,numa,1);
   
   for (irep in 1:4000) {
      hdraw = rWishart(1,tau,hbar1)[,,1];
      hdraw = solve(hdraw);
      achol = t(chol(hdraw));
      ssig = matrix(0,p,p);
      for (i in 1:p) {
         ssig[i,i] = achol[i,i]; 
         for (j in 1:p) {
            achol[j,i] = achol[j,i]/ssig[i,i];
         }
      }
      
      achol = solve(achol);
      a0draw = matrix(0,numa,1);
      ic = 1;
      for (i in 2:p) {
         for (j in 1:(i-1)) {
            a0draw[ic,1] = achol[i,j];
            ic = ic + 1;
         }
      }
      a02mo = a02mo + a0draw%*%t(a0draw);
      a0mean = a0mean + a0draw; 
   }
   a02mo = a02mo/4000;
   a0mean = a0mean/4000;
   a02mo = a02mo - a0mean%*%t(a0mean);
   return(list(aols=aols,vbar=vbar,a0=a0,ssig1=ssig1,a02mo=a02mo))
}
carter_kohn=function(y,Z,Ht,Qt,m,p,t,B0,V0,kdraw=NULL) {
   # Carter and Kohn (1994), On Gibbs sampling for state space models.
   # Kalman Filter
   if (is.null(kdraw)) {
      kdraw=c(ones(t,1))
   }
   bp = B0;
   Vp = V0;
   bt = zeros(t,m);
   Vt = zeros(m^2,t);
   R = Ht;
   log_lik = 0;
   for (i in 1:t) {
      H = Z[((i-1)*p+1):(i*p),,drop=F];
      cfe = y[,i] - H%*%bp;   # conditional forecast error
      f = H%*%Vp%*%t(H) + R;    # variance of the conditional forecast error
      inv_f = t(H)%*%solve(f);
      btt = bp + Vp%*%inv_f%*%cfe;
      Vtt = Vp - Vp%*%inv_f%*%H%*%Vp;
      if (i < t) {
         bp = btt;
         Vp = Vtt + kdraw[i]*Qt;
      }
      bt[i,] = btt;
      Vt[,i] = matrix(Vtt,m^2,1);
   }
   
   # draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
   bdraw = zeros(t,m);
   bdraw[t,] = t(btt) + rnorm(m)%*%chol(Vtt,pivot=T)
   
   # Backward recurssions
   for (i in 1:(t-1)) {
      bf = bdraw[(t-i+1),];
      btt = bt[(t-i),];
      Vtt = matrix(c(Vt[,(t-i)]),m,m,byrow=T);
      f = Vtt + kdraw[(t-i)]*Qt;
      inv_f = Vtt%*%solve(f);
      cfe = bf - btt;
      bmean = btt + inv_f%*%cfe;
      bvar = Vtt - inv_f%*%Vtt;
      bdraw[(t-i),] = t(bmean) + rnorm(m)%*%chol(bvar,pivot=T)
   }
   t(bdraw)
}

library("MASS")
library("matlab")
library("mvtnorm")
library("R.matlab")

# TVP-VAR Time varying structural VAR with constant covariance matrix
# ------------------------------------------------------------------------------------
# This code implements the Homoskedastic TVP-VAR using the Carter and Kohn (1994)
# algorithm for state-space models.
#  ************************************************************************************
# The model is:
#
#     Y(t) = B0(t) + B1(t)xY(t-1) + B2(t)xY(t-2) + u(t) 
# 
#  with u(t)~N(0,H).
# The state equation is
#
#            B(t) = B(t-1) + error
#
# where B(t) = [B0(t),B1(t),B2(t)]'.
#
# ------------------------------------------------------------------------------------

#----------------------------------LOAD DATA----------------------------------------
# Load Korobilis (2008) quarterly data
ydata = as.matrix(read.table("/Users/user/Dropbox/18. Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/R/Yraw.dat"))
# Demean and standardize data
#for (i in 1:ncol(ydata)) {
#   ydata[,i] = scale(ydata[,i],T,T)
#}
Y = ydata

# Number of observations and dimension of X and Y
t=dim(Y)[1]; # t is the time-series observations of Y
M=dim(Y)[2]; # M is the dimensionality of Y

# Number of factors & lags:
tau = 40; # tau is the size of the training sample
p = 2;    # p is number of lags in the VAR part

# ===================================| VAR EQUATION |==============================
# Generate lagged Y matrix. This will be part of the X matrix
# Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
ylag = embed(Y,p+1)[,-c(1:M)]; # Y is [T x M]. ylag is [T x (Mp)]
ylag = ylag[(tau+1):nrow(ylag),];

K = M + p*(M^2); # K is the number of elements in the state vector
# Create Z_t matrix.
Z = zeros((t-tau-p)*M,K);
for (i in 1:(t-tau-p)) {
   ztemp = eye(M);
   for (j in 1:p) {  
      xtemp = t(matrix(ylag[i,((j-1)*M+1):(j*M)]));
      xtemp = kronecker(diag(M),xtemp);
      ztemp = cbind(ztemp, xtemp);  #ok<AGROW>
   }
   Z[((i-1)*M+1):(i*M),] = ztemp;
}

# Redefine FAVAR variables y
y = t(Y[(tau+p+1):t,]);
# Time series observations
t=ncol(y);# t is now 215 - p - tau = 173

#----------------------------PRELIMINARIES---------------------------------
# Set some Gibbs - related preliminaries
nrep = 1000;     # Number of replications
nburn = 1000;      # Number of burn-in-draws
it_print = 100; # Print in the screen every "it_print"-th iteration

#========= PRIORS:
# To set up training sample prior a-la Primiceri, use the following subroutine
PRIOR = ts_prior(Y,tau,M,p);
B_OLS = PRIOR$aols
VB_OLS = PRIOR$vbar
A_OLS = PRIOR$a0
VA_OLS = PRIOR$a02mo
sigma_OLS = PRIOR$ssig1

#-------- Now set prior means and variances (_prmean / _prvar)
# This is the Kalman filter initial condition for the time-varying
# parameters B(t)
# B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean = B_OLS;
B_0_prvar = 4*VB_OLS;

# Note that for IW distribution I keep the _prmean/_prvar notation...
# Q is the covariance of B(t)
# Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean = (0.01^2)*tau*VB_OLS;
Q_prvar = tau;

# Sigma is the covariance of the VAR covariance, SIGMA
# Sigma ~ IW(I,M+1)
Sigma_prmean = diag(M);
Sigma_prvar = M+1;

#========= INITIALIZE MATRICES:
# Specify covariance matrices for measurement and state equations
consQ = 0.0001;
Qdraw = consQ*diag(K);
Qchol = sqrt(consQ)*diag(K);
Btdraw = zeros(K,t);
Sigmadraw = 0.1*diag(M);

# Storage matrices for posteriors and stuff
Bt_postmean = zeros(nrep,K,t);
Qmean = zeros(nrep,K,K);
Sigmamean = zeros(nrep,M,M);
#----------------------------- END OF PRELIMINARIES ---------------------------

#====================================== START SAMPLING ========================================
#==============================================================================================
print('Number of iterations');
for (irep in 1:(nrep + nburn)){ # GIBBS iterations starts here
  # Print iterations
  if (irep%%it_print == 0) {
    print(irep);
  }
  # -----------------------------------------------------------------------------------------
  #   STEP I: Sample B_t from p(B_t|y,Sigma) (Drawing coefficient states, pp. 844-845)
  # -----------------------------------------------------------------------------------------
  Btdraw = carter_kohn(y,Z,Sigmadraw,Qdraw,K,M,t,B_0_prmean,B_0_prvar);
  tail(Z)
  tail(Btdraw)
  
  Btemp = t(Btdraw[,2:t]) - t(Btdraw[,1:(t-1)]);
  sse_2Q = matrix(0,K,K);
  for (i in 1:(t-1)) {
    sse_2Q = sse_2Q + crossprod(t(Btemp[i,]));
  }
  
  Qdraw=Qinv = solve(sse_2Q + Q_prmean);
  Qinvdraw = rWishart(1,t+Q_prvar,Qinv)[,,1];
  Qdraw = solve(Qinvdraw);
  Qchol = chol(Qdraw);
  
  # -----------------------------------------------------------------------------------------
  #   STEP I: Sample Sigma from p(Sigma|y,B_t) which is i-Wishart
  # ----------------------------------------------------------------------------------------
  yhat = zeros(M,t);
  for (i in 1:t) {
    yhat[,i] = y[,i] - Z[((i-1)*M+1):(i*M),]%*%Btdraw[,i];
  }
             
  sse_2S = zeros(M,M);
  for (i in 1:t) {
    sse_2S = sse_2S + crossprod(t(yhat[,i]));
  }
             
  Sigmadraw=Sigmainv = solve(sse_2S + Sigma_prmean);
  Sigmainvdraw = rWishart(1,t+Sigma_prvar,Sigmainv)[,,1];
  Sigmadraw = solve(Sigmainvdraw);
  Sigmachol = chol(Sigmadraw);
             
  #----------------------------SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES -----------------
  if (irep > nburn) {
    # Save only the means of B(t), Q and SIGMA. Not memory efficient to
    # store all draws (at least for B(t) which is large). If you want to
    # store all draws, it is better to save them in a file at each iteration.
    # Use the MATLAB command 'save' (type 'help save' in the command window
    # for more info)
    Bt_postmean[(irep-nburn),,] = Btdraw;
    Qmean[(irep-nburn),,] = Qdraw;
    Sigmamean[(irep-nburn),,] = Sigmadraw;
  } # END saving after burn-in results 
} #END main Gibbs loop (for irep = 1:nrep+nburn)

#=============================GIBBS SAMPLER ENDS HERE==================================
Bt_postmean1 = apply(Bt_postmean,2:3,mean);   # Posterior mean of B(t) (VAR regression coeff.)
B_05 = apply(B_postmean,2:3,quantile,0.05);
B_50 = apply(B_postmean,2:3,quantile,0.50);
B_95 = apply(B_postmean,2:3,quantile,0.95);
Sigmamean1 = apply(Sigmamean,2:3,mean)        # Posterior mean of SIGMA (VAR covariance matrix)
Sigt_05 = apply(Sigt_postmean,2:3,quantile,0.05);
Sigt_50 = apply(Sigt_postmean,2:3,quantile,0.50);
Sigt_95 = apply(Sigt_postmean,2:3,quantile,0.95);
Qmean1 = apply(Qmean,2:3,mean);               # Posterior mean of Q (covariance of B(t))

par(mfcol = c(7,3), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:dim(B_postmean1)[1]) {
   plot(B_postmean1[i,],type="l",xaxs="i",las=1,xlab="",ylab="",tck=0.02)#,ylim=c(min(B_05[i,]),max(B_95[i,])))
   #lines(B_50[i,],col="steelblue4",lty=2)
   #lines(B_05[i,],col="steelblue4")
   #lines(B_95[i,],col="steelblue4")
}

par(mfcol = c(M,1), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:M) {
   plot(Sigt_postmean1[i,],type="l",xaxs="i",las=1,ylim=c(min(Sigt_05[i,]),max(Sigt_95[i,])),xlab="",ylab="",tck=0.02)
   lines(Sigt_50[i,],col="steelblue4",lty=2)
   lines(Sigt_05[i,],col="steelblue4")
   lines(Sigt_95[i,],col="steelblue4")
}

### END
