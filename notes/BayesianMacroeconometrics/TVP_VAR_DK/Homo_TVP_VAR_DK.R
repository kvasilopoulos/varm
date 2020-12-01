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
kalfilt=function(y1,Z,Ht,Qt,m,p,t) {
   # Run the Kalman filter and then return what -- mean of w 
   #Kalman filter code
   Kkeep = zeros(m*t,p);
   Lkeep = zeros(m*t,m);
   Fkeep = zeros(p*t,p);
   a = zeros(m,t+1);
   v = zeros(p,t);
   Pt = zeros(m,m);
   llik = 0;
   for (i in 1:t) {
      htemp = Ht[((i-1)*p+1):(i*p),];
      ztemp = Z[((i-1)*p+1):(i*p),];
      v[,i] = y1[,i] - ztemp%*%a[,i];
      Ft = ztemp%*%Pt%*%t(ztemp) + htemp;
      Ftinv = solve(Ft);
      llik = llik + log(det(Ft)) + v[,i]%*%Ftinv%*%v[,i];
      Fkeep[((i-1)*p+1):(i*p),] = Ftinv;
      Kt = Pt%*%t(ztemp)%*%Ftinv ;
      Kkeep[((i-1)*m+1):(i*m),] = Kt;
      Ltt = eye(m) - Kt%*%ztemp;
      Lkeep[((i-1)*m+1):(i*m),] = Ltt;
      a[,(i+1)] = a[,i] + Kt%*%v[,i];
      Pt = Pt%*%t(Ltt) + Qt;
   }
   llik = -.5*llik;
   
   #Backward recursion to evaluate rt and, thus, whatt*/
   rt = zeros(m,t+1);
   pm = p+m;
   what = zeros(pm*t,1);
   for (i in t:1) {
      htemp = Ht[((i-1)*p+1):(i*p),];
      ztemp = Z[((i-1)*p+1):(i*p),];
      lterm = Lkeep[((i-1)*m+1):(i*m),];
      fterm = Fkeep[((i-1)*p+1):(i*p),];
      kterm = Kkeep[((i-1)*m+1):(i*m),];
      what[((i-1)*pm+1):((i-1)*pm+p),1] = htemp%*%fterm%*%v[,i] - htemp%*%t(kterm)%*%rt[,(i+1)];
      what[((i-1)*pm+p+1):(i*pm),1] = Qt%*%rt[,(i+1)];
      rt[,i] = t(ztemp)%*%fterm%*%v[,i] + t(lterm)%*%rt[,(i+1)];
   }

   alph = zeros(m,t+1);
   for (i in 1:t) {
      alph[,(i+1)] = alph[,i] + Qt%*%rt[,(i+1)];
   }
   return = list(what=what,alph=alph,llik=llik)
}
dk=function(y,p,m,t,Qchol,Ht,Qt,Z) {
   #First draw w as in page 605 of DK
   pm = p + m;
   wplus = zeros(pm*t,1);
   for (i in 1:t) {
      Hchol = chol(Ht[((i-1)*p+1):(i*p),]);
      wplus[((i-1)*pm+1):((i-1)*pm+p),1] = t(Hchol)%*%rnorm(p);
      wplus[((i-1)*pm+p+1):((i-1)*pm+pm),1] = t(Qchol)%*%rnorm(m);
   }
   
   #Now get implied draw of y
   aplus = zeros(m,t+1);
   pm = p + m;
   yplus = zeros(p,t);
   for (i in 1:t) {
      ztemp = Z[((i-1)*p+1):(i*p),];
      yplus[,i] = ztemp%*%aplus[,i] + wplus[((i-1)*pm+1):((i-1)*pm+p),1];
      aplus[,(i+1)] = aplus[,i] + wplus[((i-1)*pm+p+1):(i*pm),1] ;
   }

   kf1 = kalfilt(y,Z,Ht,Qt,m,p,t);
   #what=kf1$what
   ahat=kf1$alph
   
   kf2 = kalfilt(yplus,Z,Ht,Qt,m,p,t);
   #whatp=kf2$what
   ahatp=kf2$alph
   tail(ahat)
   tail(ahatp)
   atilda = ahat - ahatp + aplus;
   atilda
}

library("MASS")
library("matlab")
library("mvtnorm")
library("R.matlab")

# TVP-VAR Time varying structural VAR with constant covariance matrix
# ------------------------------------------------------------------------------------
# This code implements the Homoskedastic TVP-VAR using the Durbin and Koopman (2002)
# algorithm for state-space models.
# ************************************************************************************
# The model is:
#
#  Y(t) = B0(t) + B1(t)xY(t-1) + B2(t)xY(t-2) + u(t) 
# 
#  with u(t)~N(0,H).
# The state equation is
#
#B(t) = B(t-1) + error
#
# where B(t) = [B0(t),B1(t),B2(t)]'.
#
# ************************************************************************************
#NOTE: 
#There are references to equations of Primiceri, "Time Varying Structural Vector
#Autoregressions & Monetary Policy",(2005),Review of Economic Studies 72,821-852
#for your convenience. The definition of vectors/matrices is also based on thispaper. 
#VERSION:
#This version: June 15, 2008
#AUTHOR: 
#Dimitris Korobilis
#Department of Economics, University of Strathclyde
#Enjoy! (...and extend/modify)
# ------------------------------------------------------------------------------------

#----------------------------------LOAD DATA----------------------------------------
# Load Korobilis (2008) quarterly data
#load ydata.dat;
#load yearlab.dat;

# # Demean and standardize data
# t2 = size(ydata,1);
# stdffr = std(ydata(:,3));
# ydata = (ydata- repma2t(mean(ydata,1),t2,1))./repmat(std(ydata,1),t2,1);

#Y=ydata;

Y = as.matrix(read.table("/Users/user/Dropbox/18. Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/R/Yraw.dat"))
head(Y)
tail(Y)

# Number of observations and dimension of X and Y
t=nrow(Y);  # t is the time-series observations of Y
M=ncol(Y);  # M is the dimensionality of Y

# Number of factors & lags:
tau = 40; # tau is the size of the training sample
p = 2; # p is number of lags in the VAR part
numa = M*(M-1)/2; # Number of lower triangular elements of A_t (other than 0's and 1's)
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
#yearlab = yearlab[(tau+p+1):t];
# Time series observations
t=ncol(y);# t is now 215 - p - tau = 173

#----------------------------PRELIMINARIES---------------------------------
# Set some Gibbs - related preliminaries
nburn = 500;# Number of burn-in-draws
nrep = 500;  # Number of replications
it_print = 100;  #Print in the screen every "it_print"-th iteration

#========= PRIORS:
# To set up training sample prior a-la Primiceri, use the following subroutine
PRIOR = ts_prior(Y,tau,M,p)
B_OLS=PRIOR$aols
VB_OLS=PRIOR$vbar
A_OLS=PRIOR$a0
sigma_OLS=PRIOR$ssig1
VA_OLS=PRIOR$a02mo

# Or use uninformative values
#B_OLS = zeros(K,1);
#VB_OLS = eye(K);

#-------- Now set prior means and variances (_prmean / _prvar)
# This is the Kalman filter initial condition for the time-varying
# parameters B(t)
# B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean = B_OLS;
B_0_prvar = 4*VB_OLS;
inv_B_0_prvar = solve(B_0_prvar);

# Note that for IW distribution I keep the _prmean/_prvar notation...
# Q is the covariance of B(t)
# Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean = 0.0001*tau*VB_OLS;
Q_prvar = tau;

# Sigma is the covariance of the VAR covariance, SIGMA
# Sigma ~ IW(I,M+1)
Sigma_prmean = eye(M);
Sigma_prvar = M+1;

#========= INITIALIZE MATRICES:
# Specify covariance matrices for measurement and state equations
consQ = 0.0001;
Qdraw = consQ*eye(K);
Qchol = sqrt(consQ)*eye(K);
Btdraw = zeros(K,t+1);
Btdrawc = zeros(K,t+1);
Sigmadraw = 0.1*eye(M);
Sigmainvdraw = solve(Sigmadraw);
Ht = zeros(M*t,M);  
for (i in 1:t) {
   Ht[((i-1)*M+1):(i*M),] = Sigmadraw;
}

# Storage matrices for posteriors and stuff
Bt_postmean = zeros(nrep,K,t+1);
Qmean = zeros(nrep,K,K);
Sigmamean = zeros(nrep,M,M);
#----------------------------- END OF PRELIMINARIES ---------------------------
 
#====================================== START SAMPLING ========================================
#==============================================================================================
tic(); # This is just a timer
print('Number of iterations');
for (irep in 1:(nrep+nburn)) { # GIBBS iterations starts here
   # Print iterations
   if (irep%%it_print==0) {
      print(paste0(round(irep/(nrep+nburn-1)*100,2),"%"));
      toc();
   }
   # -----------------------------------------------------------------------------------------
   #STEP I: Sample B from p(B|y,A,Sigma,V) (Drawing coefficient states, pp. 844-845)
   # -----------------------------------------------------------------------------------------
   # Draw initial condition for B(t)
   vbar=zeros(K,K);
   xhy=zeros(K,1);
   for (i in 1:t) {
      zhat1 = Z[((i-1)*M+1):(i*M),];
      yhat1 = y[,i] - zhat1%*%Btdrawc[,i];
      HHat1 = Sigmainvdraw;  #solve(Ht[((i-1)*M+1):(i*M),]);  #use the commented code if you are estimating the Heteroskedastic TVP-VAR
      vbar = vbar + t(zhat1)%*%HHat1%*%zhat1;
      xhy = xhy + t(zhat1)%*%HHat1%*%yhat1;
   }
   vbar = solve(vbar + inv_B_0_prvar);
   B0hat = vbar%*%(inv_B_0_prvar%*%B_0_prmean + xhy);
   B0draw = B0hat + t(chol(vbar))%*%rnorm(K); # Draw from the initial condition B(0)
 
   ya = zeros(M,t);
   for (i in 1:t) {
      ya[,i] = y[,i] - Z[((i-1)*M+1):(i*M),]%*%B0draw;
   }
   
   # Now get a draw of B(t) using the Durbin and Koopman (2002) smoother
   Btdrawc = dk(ya,M,K,t,Qchol,Ht,Qdraw,Z);
   
   # Add on the initial condition B(0)
   for (i in 1:(t+1)) {
      Btdraw[,i] = Btdrawc[,i] + B0draw;
   }
   
   # Now get the SSE in the state equation (to estimate the covariance Q)
   Btemp = t(Btdraw[,2:(t+1)]) - t(Btdraw[,1:t]);
   sse_2 = zeros(K,K);
   for (i in 1:t) {
      sse_2 = sse_2 + crossprod(t(Btemp[i,]));
   }

   # Draw Q, the coavariance matrix of the time-varying parameters B(t)
   Qinv = solve(sse_2 + Q_prmean);
   Qinvdraw = rWishart(1,t-1+Q_prvar,Qinv)[,,1];
   Qdraw = solve(Qinvdraw);# Qdraw is a draw from the posterior of Q
   Qchol = chol(Qdraw);
 
   # -----------------------------------------------------------------------------------------
   #STEP II: Sample Sigma from p(Sigma|y,B_t) which is i-Wishart
   # ----------------------------------------------------------------------------------------
   # Get SSE of the VAR model
   yhat = zeros(M,t);
   for (i in 1:t) {
      yhat[,i] = y[,i] - Z[((i-1)*M+1):(i*M),]%*%Btdraw[,i];
   }

   sse_2S = zeros(M,M);
   for (i in 1:t) {
      sse_2S = sse_2S + crossprod(t(yhat[,i]));
   }
   # Draw SIGMA, the VAR covariance matrix
   Sigmainv = solve(sse_2S + Sigma_prmean);
   Sigmainvdraw = rWishart(1,t+Sigma_prvar,Sigmainv)[,,1];
   Sigmadraw = solve(Sigmainvdraw);  # Sigmadraw is a draw from the posterior of SIGMA
   Sigmachol = chol(Sigmadraw);
 
   # Replicate SIGMA into a [Txp x M] matrix (not very efficient, but 
   # usefull if you want to extend to a time-varying covariance matrix SIGMA,
   # without having to interfere with the smoothing function dk.K).
   Ht = zeros(M*t,M);
   for (i in 1:t) {
      Ht[((i-1)*M+1):(i*M),] = Sigmadraw;
   }

   #----------------------------SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES -----------------
   if (irep>nburn) {
      # Save only the means of B(t), Q and SIGMA. Not memory efficient to
      # store all draws (at least for B(t) which is large). If you want to
      # store all draws, it is better to save them in a file at each iteration.
      # Use the MATLAB command 'save' (type 'help save' in the command window for more info)
      Bt_postmean[(irep-nburn),,] = Btdraw;
      Qmean[(irep-nburn),,] = Qdraw;
      Sigmamean[(irep-nburn),,] = Sigmadraw;
   } # END saving after burn-in results 
} #END main Gibbs loop (for irep = 1:nrep+nburn)
toc(); # Stop timer and print total time
#=============================GIBBS SAMPLER ENDS HERE==================================

B_t = apply(Bt_postmean,2:3,mean) # Posterior mean of B(t) (VAR regression coeff.)
B_10 = apply(Bt_postmean,2:3,quantile,0.1)
B_50 = apply(Bt_postmean,2:3,quantile,0.5)
B_90 = apply(Bt_postmean,2:3,quantile,0.9)

par(mfrow = c(ceiling(dim(B_t)[1]/3),3), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:dim(B_t)[1]) {
   plot(B_t[i,],type="l",las=1,xaxs="i",xlab="",ylab="",tck=.02,ylim=c(min(B_10[i,]),max(B_90[i,])))
   lines(B_10[i,],col=2)
   lines(B_50[i,],col=2,lty=2)
   lines(B_90[i,],col=2)
   abline(h=0,lty=2)
}

Q = apply(Qmean,2:3,mean) # Posterior mean of Q (covariance of B(t))
Q
Sigma = apply(Sigmamean,2:3,mean)  # Posterior mean of SIGMA (VAR covariance matrix)
Sigma

### END
