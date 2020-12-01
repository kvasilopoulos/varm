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
   log_lik = 0;
   for (i in 1:t) {
      R = Ht[((i-1)*p+1):(i*p),,drop=F];
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
   bdraw[t,] = t(btt) + rnorm(m)%*%chol(Vtt)
   
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
      bdraw[(t-i),] = t(bmean) + rnorm(m)%*%chol(bvar)
   }
   t(bdraw)
}

library("MASS")
library("matlab")
library("mvtnorm")
library("R.matlab")

# TVP-VAR Time varying structural VAR with stochastic volatility
# ------------------------------------------------------------------------------------
# This code implements the TVP-VAR model as in Primiceri (2005). See also
# the monograph, Section 4.2 and Section 3.3.2.
# ************************************************************************************
# The model is:
#
#     Y(t) = B0(t) + B1(t)xY(t-1) + B2(t)xY(t-2) + e(t) 
# 
#  with e(t) ~ N(0,SIGMA(t)), and  L(t)' x SIGMA(t) x L(t) = D(t)*D(t),
#             _                                          _
#            |    1         0        0       ...       0  |
#            |  L21(t)      1        0       ...       0  |
#    L(t) =  |  L31(t)     L32(t)    1       ...       0  |
#            |   ...        ...     ...      ...      ... |
#            |_ LN1(t)      ...     ...    LN(N-1)(t)  1 _|
# 
# 
# and D(t) = diag[exp(0.5 x h1(t)), .... ,exp(0.5 x hn(t))].
#
# The state equations are
#
#            B(t) = B(t-1) + u(t),            u(t) ~ N(0,Q)
#            l(t) = l(t-1) + zeta(t),      zeta(t) ~ N(0,S)
#            h(t) = h(t-1) + eta(t),        eta(t) ~ N(0,W)
#
# where B(t) = [B0(t),B1(t),B2(t)]', l(t)=[L21(t),...,LN(N-1)(t)]' and
# h(t) = [h1(t),...,hn(t)]'.
#
# ************************************************************************************
#   NOTE: 
#      There are references to equations of Primiceri, "Time Varying Structural Vector
#      Autoregressions & Monetary Policy",(2005),Review of Economic Studies 72,821-852
#      for your convenience. The definition of vectors/matrices is also based on this
#      paper.
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
numa = M*(M-1)/2; # Number of lower triangular elements of A_t (other than 0's and 1's)
# ===================================| VAR EQUATION |==============================
# Generate lagged Y matrix. This will be part of the X matrix
ylag = embed(Y,(p+1))[,-c(1:ncol(Y))]; # Y is [T x M]. ylag is [T x (Mp)]
# Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
ylag = ylag[(tau+1):(t-p),]

K = M + p*(M^2); # K is the number of elements in the state vector
# Create Z_t matrix.
Z = zeros((t-tau-p)*M,K);
for (i in 1:(t-tau-p)) {
   ztemp = eye(M);
   for (j in 1:p) {  
      xtemp = t(matrix(ylag[i,((j-1)*M+1):(j*M)]));
      xtemp = kronecker(diag(M),xtemp);
      ztemp = cbind(ztemp, xtemp);  
   }
   Z[((i-1)*M+1):(i*M),] = ztemp;
}

# Redefine FAVAR variables y
y = t(Y[(tau+p+1):t,]);
# Time series observations
t=ncol(y);# t is now 215 - p - tau = 173

#----------------------------PRELIMINARIES---------------------------------
# Set some Gibbs - related preliminaries
nrep = 500;     # Number of replications
nburn = 500;    # Number of burn-in-draws
it_print = 100; # Print in the screen every "it_print"-th iteration

#========= PRIORS:
# To set up training sample prior a-la Primiceri, use the following subroutine
# Or use uninformative values
# 1. Uninformative Prior
#A_OLS = zeros(numa,1);
#B_OLS = zeros(K,1);
#VA_OLS = diag(numa);
#VB_OLS = diag(K);
#sigma_OLS = 0*ones(M,1);

# To set up training sample prior a-la Primiceri, use the following subroutine
PRIOR = ts_prior(Y,tau,M,p);
B_OLS = PRIOR$aols
VB_OLS = PRIOR$vbar
A_OLS = PRIOR$a0
sigma_OLS = PRIOR$ssig1
VA_OLS = PRIOR$a02mo

# Set some hyperparameters here (see page 831, end of section 4.1)
k_Q = 0.01;
k_S = 0.1;
k_W = 1;

# We need the sizes of some matrices as prior hyperparameters (see page
# 831 again, lines 2-3 and line 6)
sizeW = M;   # Size of matrix W
sizeS = 1:M; # Size of matrix S

#-------- Now set prior means and variances (_prmean / _prvar)
# These are the Kalman filter initial conditions for the time-varying
# parameters B(t), A(t) and (log) SIGMA(t). These are the mean VAR
# coefficients, the lower-triangular VAR covariances and the diagonal
# log-volatilities, respectively 
# B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean = B_OLS;
B_0_prvar = 4*VB_OLS;

# A_0 ~ N(A_OLS, 4Var(A_OLS))
A_0_prmean = A_OLS;
A_0_prvar = 4*VA_OLS;

# log(sigma_0) ~ N(log(sigma_OLS),I_n)
sigma_prmean = sigma_OLS;
sigma_prvar = 4*diag(M);

# Note that for IW distribution I keep the _prmean/_prvar notation....
# Q is the covariance of B(t), S is the covariance of A(t) and W is the
# covariance of (log) SIGMA(t)
# Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean = ((k_Q)^2)*tau*VB_OLS;
Q_prvar = tau;

# W ~ IG(k2_W*(1+dimension(W))*I_n,(1+dimension(W)))
W_prmean = ((k_W)^2)*ones(M,1);
W_prvar = 2;

# S ~ IW(k2_S*(1+dimension(S)*Var(A_OLS),(1+dimension(S)))
S_prmean = cell(M-1,1);
S_prvar = zeros(M-1,1);
ind = 1;

for (ii in 2:M) {
  # S is block diagonal as in Primiceri (2005)
  S_prmean[[ii-1]] = ((k_S)^2)*(1 + sizeS[ii-1])*VA_OLS[((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind];
  S_prvar[ii-1] = 1 + sizeS[ii-1];
  ind = ind + ii;
}

#========= INITIALIZE MATRICES:
# Specify covariance matrices for measurement and state equations
consQ = 0.0001;
consS = 0.0001;
consH = 0.01;
consW = 0.0001;
Ht = kronecker(ones(t,1),consH*diag(M));   # Initialize Htdraw, a draw from the VAR covariance matrix
Htchol = kronecker(ones(t,1),sqrt(consH)*diag(M)); # Cholesky of Htdraw defined above
Qdraw = consQ*diag(K);     # Initialize Qdraw, a draw from the covariance matrix Q
Sdraw = consS*diag(numa);  # Initialize Sdraw, a draw from the covariance matrix S
Sblockdraw = cell(M-1,1); # ...and then get the blocks of this matrix (see Primiceri)
ijc = 1;
for (jj in 2:M) {
  Sblockdraw[[jj-1]] = Sdraw[((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc];
  ijc = ijc + jj;
}
Wdraw = consW*ones(M,1); # Initialize Wdraw, a draw from the covariance matrix W
Btdraw = zeros(K,t);     # Initialize Btdraw, a draw of the mean VAR coefficients, B(t)
Atdraw = zeros(numa,t);  # Initialize Atdraw, a draw of the non 0 or 1 elements of A(t)
Sigtdraw = zeros(t,M);   # Initialize Sigtdraw, a draw of the log-diagonal of SIGMA(t)
sigt = kronecker(ones(t,1),0.01*diag(M));   # Matrix of the exponent of Sigtdraws (SIGMA(t))
statedraw = 5*ones(t,M);       # initialize the draw of the indicator variable # (of 7-component mixture of Normals approximation)
Zs = kronecker(ones(t,1),diag(M));

# Storage matrices for posteriors and stuff
Bt_postmean = zeros(nrep,K,t);    # regression coefficients B(t)
At_postmean = zeros(nrep,numa,t); # lower triangular matrix A(t)
Sigt_postmean = zeros(nrep,t,M);  # diagonal std matrix SIGMA(t)
Qmean = zeros(nrep,K,K);          # covariance matrix Q of B(t)
Smean = zeros(nrep,numa,numa);    # covariance matrix S of A(t)
Wmean = zeros(nrep,M,1);          # covariance matrix W of SIGMA(t)

sigmean = zeros(nrep,t,M);    # mean of the diagonal of the VAR covariance matrix
cormean = zeros(nrep,t,numa); # mean of the off-diagonal elements of the VAR cov matrix
sig2mo = zeros(nrep,t,M);     # squares of the diagonal of the VAR covariance matrix
cor2mo = zeros(nrep,t,numa);  # squares of the off-diagonal elements of the VAR cov matrix
#----------------------------- END OF PRELIMINARIES ---------------------------

#==============================================================================================
#====================================== START SAMPLING ========================================
#==============================================================================================
print('Number of iterations');
for (irep in 1:(nrep + nburn)) {   # GIBBS iterations starts here
  # Print iterations
  if (irep%%it_print == 0) {
    print(irep);
  }
  
  # -----------------------------------------------------------------------------------------
  #   STEP I: Sample B from p(B|y,A,Sigma,V) (Drawing coefficient states, pp. 844-845)
  # -----------------------------------------------------------------------------------------

  Btdrawc = carter_kohn(y,Z,Ht,Qdraw,K,M,t,B_0_prmean,B_0_prvar);
  Btdraw = Btdrawc;
  
  #=====| Draw Q, the covariance of B(t) (from iWishart)
  # Take the SSE in the state equation of B(t)
  Btemp = t(Btdraw[,2:t]) - t(Btdraw[,1:(t-1)]);
  
  sse_2 = matrix(0,K,K);
  for (i in 1:(t-1)){
    sse_2 = sse_2 + matrix(Btemp[i,],ncol=1)%*%matrix(Btemp[i,],nrow=1);
  }
  
  # ...and subsequently draw Q, the covariance matrix of B(t)
  Qinv = solve(sse_2 + Q_prmean);
  Qinvdraw = rWishart(1,t+Q_prvar,Qinv)[,,1];
  Qdraw = solve(Qinvdraw);  # this is a draw from Q
  
  #-------------------------------------------------------------------------------------------
  #   STEP II: Draw A(t) from p(At|y,B,Sigma,V) (Drawing coefficient states, p. 845)
  #-------------------------------------------------------------------------------------------

  # Substract from the data y(t), the mean Z x B(t)
  yhat = matrix(0,M,t);
  for (i in 1:t) {
    yhat[,i] = matrix(y[,i],ncol=1) - Z[((i-1)*M+1):(i*M),]%*%Btdraw[,i];
  }
  
  # This part is more tricky, check Primiceri
  # Zc is a [M x M(M-1)/2] matrix defined in (A.2) page 845, Primiceri
  Zc = - t(yhat);
  sigma2temp = exp(Sigtdraw);
  
  Atdraw = NULL;
  ind = 1;
  for (ii in 2:M) {
    # Draw each block of A(t)
    Atblockdraw = carter_kohn(t(matrix(yhat[ii,])),Zc[,1:(ii-1),drop=F],sigma2temp[,ii,drop=F],Sblockdraw[[(ii-1)]],sizeS[(ii-1)],1,t,A_0_prmean[((ii-1)+(ii-3)*(ii-2)/2):ind,,drop=F],A_0_prvar[((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind,drop=F]);
    Atdraw = rbind(Atdraw, Atblockdraw); # Atdraw is the final matrix of draws of A(t)
    ind = ind + ii;
  }
  
  #=====| Draw S, the covariance of A(t) (from iWishart)
  # Take the SSE in the state equation of A(t)
  Attemp = t(Atdraw[,2:t, drop=F]) - t(Atdraw[,1:(t-1),drop=F]);
  sse_2 = matrix(0,numa,numa);
  for (i in 1:(t-1)) {
    sse_2 = sse_2 + crossprod(t(Attemp[i,]));
  }
  # ...and subsequently draw S, the covariance matrix of A(t) 
  ijc = 1;
  for (jj in 2:M) {
    Sinv = solve(sse_2[((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc] + S_prmean[[jj-1]]);
    Sinvblockdraw = rWishart(1,t+S_prvar[jj-1],Sinv)[,,1];
    Sblockdraw[[jj-1]] = solve(Sinvblockdraw); # this is a draw from S
    ijc = ijc + jj;
  }


  #------------------------------------------------------------------------------------------
  #   STEP III: Draw diagonal VAR covariance matrix log-SIGMA(t)
  #------------------------------------------------------------------------------------------

  capAt = matrix(0,M*t,M);
  for (i in 1:t) {
    capatemp = diag(M);
    aatemp = Atdraw[,i,drop=F];
    ic=1;
    for (j in 2:M) {
      capatemp[j,1:(j-1)] = t(aatemp[ic:(ic+j-2),1,drop=F]);
      ic = ic + j - 1;
    }
    capAt[((i-1)*M+1):(i*M),] = capatemp;
  }
  
  # yhat is the vector y(t) - Z x B(t) defined previously. Multiply yhat
  # with capAt, i.e the lower triangular matrix A(t). Then take squares
  # of the resulting quantity (saved in matrix y2)
  y2 = NULL;
  for (i in 1:t) {
    ytemps = capAt[((i-1)*M+1):(i*M),,drop=F]%*%yhat[,i,drop=F];
    y2 = cbind(y2,  ytemps^2); 
  }
  
  yss = t(log(y2 + 1e-6));
  for (j in 1:M) {
     T = length(Sigtdraw[,j]);
     # normal mixture
     pi = c(0.0073, .10556, .00002, .04395, .34001, .24566, .2575);
     mi = c(-10.12999, -3.97281, -8.56686, 2.77786, .61942, 1.79518, -1.08819) - 1.2704; # means already adjusted!! %%
     sigi = c(5.79596, 2.61369, 5.17950, .16735, .64009, .34023, 1.26261);
     sqrtsigi = sqrt(sigi);
     
     # sample S from a 7-point distrete distribution
     temprand = runif(T);
     
     q = repmat(pi,T,1)*dnorm(repmat(matrix(yss[,j]),1,7),repmat(matrix(Sigtdraw[,j]),1,7)+repmat(matrix(mi),T,1), repmat(matrix(sqrtsigi),T,1));
     q = q/repmat(rowSums(q),1,7);
     S = 7 - apply(repmat(matrix(temprand),1,7)<t(apply(q,1,cumsum)),1,sum)+1;
     
     vart = zeros(1,1,T);
     yss1 = zeros(T,1);
     for (i in 1:T) {
        imix = S[i];
        vart[1,1,i] = sigi[imix];
        yss1[i,1] = yss[i,j] - mi[imix];
     }
     h = carter_kohn(t(matrix(yss1)),ones(T,1),matrix(c(vart)),matrix(Wdraw[j,,drop=F]),1,1,T,sigma_prmean[j,drop=F],sigma_prvar[j,j,drop=F]);
     Sigtdraw[,j]  = h
     statedraw[,j] = S
  }
  sigt = exp(.5*Sigtdraw);
  
  e2 = Sigtdraw[2:nrow(Sigtdraw),,drop=F] - Sigtdraw[1:(nrow(Sigtdraw)-1),,drop=F];
  W1 = rep(W_prvar+t-p-1,nrow(Wdraw))
  W2 = W_prmean + apply(e2^2,2,sum)
  Winvdraw = Wdraw
  for (i in 1:nrow(Wdraw)) {
    Winvdraw[i,] = rgamma(1,W1[i]/2,2/W2[i]);
  }
  Wdraw = 1/Winvdraw;
  
  # Create the VAR covariance matrix H(t). It holds that: A(t) x H(t) x A(t)' = SIGMA(t) x SIGMA(t) '
  Ht = matrix(0,M*t,M);
  Htsd = matrix(0,M*t,M);
  for (i in 1:t) {
    inva = solve(capAt[((i-1)*M+1):(i*M),]);
    stem = diag(sigt[i,]);
    Hsd = inva%*%stem;
    Hdraw = Hsd%*%t(Hsd);
    Ht[((i-1)*M+1):(i*M),] = Hdraw;  # H(t)
    Htsd[((i-1)*M+1):(i*M),] = Hsd;  # Cholesky of H(t)
  }

  #----------------------------SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES -----------------
  if (irep > nburn) {               
    # Save only the means of parameters. Not memory efficient to
    # store all draws (at least for the time-varying parameters vectors,
    # which are large). If you want to store all draws, it is better to
    # save them in a file at each iteration. Use the MATLAB command 'save'
    # (type 'help save' in the command window for more info)
    Bt_postmean[(irep-nburn),,] = Btdraw;   # regression coefficients B(t)
    At_postmean[(irep-nburn),,] = Atdraw;   # lower triangular matrix A(t)
    Sigt_postmean[(irep-nburn),,] = Sigtdraw;  # diagonal std matrix SIGMA(t)
    Qmean[(irep-nburn),,] = Qdraw;     # covariance matrix Q of B(t)
    ikc = 1;
    for (kk in 2:M) {
      Sdraw[((kk-1)+(kk-3)*(kk-2)/2):ikc,((kk-1)+(kk-3)*(kk-2)/2):ikc]=Sblockdraw[[(kk-1)]];
      ikc = ikc + kk;
    }
    Smean[(irep-nburn),,] = Sdraw;    # covariance matrix S of A(t)
    Wmean[(irep-nburn),,] = Wdraw;    # covariance matrix W of SIGMA(t)
    
    # Get time-varying correlations and variances
    stemp6 = zeros(M,1);
    stemp5 = NULL;
    stemp7 = NULL;
    for (i in 1:t) {
       std = sqrt(diag(Ht[((i-1)*M+1):(i*M),]));
       stemp8 = Ht[((i-1)*M+1):(i*M),]/(std%*%t(std));
       stemp7a = NULL;
       ic = 1;
       for (j in 1:M) {
          if (j>1) {
             stemp7a = c(stemp7a,stemp8[j,1:ic]) 
             ic = ic+1;
          }
          stemp6[j,1] = sqrt(Ht[((i-1)*p+j),j]);
       }
       stemp5 = rbind(stemp5, t(stemp6)); 
       stemp7 = rbind(stemp7, t(stemp7a)); 
    }
    sigmean[(irep-nburn),,] = stemp5;
    cormean[(irep-nburn),,] = stemp7; 
    sig2mo[(irep-nburn),,] = stemp5^2;
    cor2mo[(irep-nburn),,] = stemp7^2;
  } # END saving after burn-in results 
} # END main Gibbs loop (for irep = 1:nrep+nburn)

#=============================GIBBS SAMPLER ENDS HERE==================================
Bt_postmean1 = apply(Bt_postmean,2:3,mean);  # Posterior mean of B(t) (VAR regression coeff.)
Bt_05 = apply(Bt_postmean,2:3,quantile,0.05)
Bt_95 = apply(Bt_postmean,2:3,quantile,0.95)
Sigt_postmean1 = apply(Sigt_postmean,2:3,mean);  # Posterior mean of SIGMA(t) (VAR variances)
Sigt_05 = apply(Sigt_postmean,2:3,quantile,0.05);
Sigt_95 = apply(Sigt_postmean,2:3,quantile,0.95);

At_postmean1 = apply(At_postmean,2:3,mean);  # Posterior mean of A(t) (VAR covariances)
Qmean1 = apply(Qmean,2:3,mean);   # Posterior mean of Q (covariance of B(t))
Smean1 = apply(Smean,2:3,mean);   # Posterior mean of S (covariance of A(t))
Wmean1 = apply(Wmean,2:3,mean);   # Posterior mean of W (covariance of SIGMA(t))
          
sigmean1 = apply(sigmean,2:3,mean);
cormean1 = apply(cormean,2:3,mean);
sig2mo1 = apply(sig2mo,2:3,mean);
cor2mo1 = apply(cor2mo,2:3,mean);

par(mfcol = c(7,3), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:nrow(Bt_postmean1)) {
   plot(Bt_postmean1[i,],type="l",xaxs="i",las=1,xlab="",ylab="",ylim=c(min(Bt_05[i,]),max(Bt_95[i,])),tck=0.02)
   lines(Bt_05[i,],col="steelblue4")
   lines(Bt_95[i,],col="steelblue4")
   abline(h=0,lty=2)
}

par(mfcol = c(M,1), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:M) {
   plot(Sigt_postmean1[,i],type="l",xaxs="i",las=1,xlab="",ylab="",ylim=c(min(Sigt_05[,i]),max(Sigt_95[,i])),tck=0.02)
   lines(Sigt_05[,i],col="steelblue4")
   lines(Sigt_95[,i],col="steelblue4")
   abline(h=0,lty=2)
}

### END
