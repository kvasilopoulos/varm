
library("MASS")
library("matlab")
library("mvtnorm")
library("R.matlab")
hierch_tvp_dgp = function(T=NULL,M=NULL,p=NULL) {
   if (is.null(T)) {
      T = 200;
   }
   if (is.null(M)) {
      M = 3;
   }
   if (is.null(p)) {
      p = 1;
   }
   K = p*(M^2);
   
   Q = 0.01*eye(K);
   R = 0.001*eye(K);
   Sigma = matrix(c(1.0000,-0.50,-0.25,
                    -0.50,1.25,-0.375,
                    -0.25,-0.375,1.3125),ncol=M);
   
   theta_t = zeros(T+p,K);
   beta_t = zeros(T+p,K);
   A_0 = eye(K);
   
   theta_0 = matrix(c(0.7,0,0.35,0,0.7,0,0,0.65,0.7))
   for (i in 1:(T+p)) {
      if (i==1) {
         theta_t[i,] = t(theta_0) +  rnorm(K)%*%chol(R);
      } else {
         theta_t[i,] = theta_t[(i-1),] + rnorm(K)%*%chol(R);
      }
   }
   
   for (i in 1:(T+p)) {
      beta_t[i,] = theta_t[i,]%*%A_0 + rnorm(K)%*%chol(Q);        
   }
   
   y = rbind(runif(p*M), zeros(T,M));
   for (i in (p+1):(T+p)) {
      y[i,] = y[(i-1),]%*%matrix(beta_t[i,],M,M) + rnorm(M)%*%chol(Sigma);
   }
   x=y[(p+1):T,];
   return = list(x=x,beta_t=beta_t,theta_t=theta_t)
}
carter_kohn_hom2 = function(y,Z,Ht,Qt,m,p,t,B0,V0) {
  # Carter and Kohn (1994), On Gibbs sampling for state space models.

  # Kalman Filter
  bp = B0;
  Vp = V0;
  bt = matrix(0,t,m);
  Vt = matrix(0,m^2,t);
  log_lik = 0;
  R = Ht;
  H = Z;
  for (i in 1:t) {
    cfe = y[,i] - H%*%bp;   # conditional forecast error
    f = H%*%Vp%*%t(H) + R;    # variance of the conditional forecast error
    inv_f = ginv(f);
    #log_lik = log_lik + log(det(f)) + t(cfe)%*%inv_f%*%cfe;
    btt = bp + Vp%*%t(H)%*%inv_f%*%cfe;
    Vtt = Vp - Vp%*%t(H)%*%inv_f%*%H%*%Vp;
    if (i < t) {
      bp = btt;
      Vp = Vtt + Qt;
    }
    bt[i,] = t(btt);
    Vt[,i] = matrix(Vtt,m^2,1);
  }

  # draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
  bdraw = matrix(0,t,m);
  bdraw[t,] = rmvnorm(1,btt,Vtt);
  
  # Backward recurssions
  for (i in 1:(t-1)) {
    bf = t(bdraw[(t-i+1),]);
    btt = t(bt[(t-i),]);
    Vtt = matrix(Vt[,(t-i)],m,m);
    f = Vtt + Qt;
    inv_f = ginv(f);
    cfe = bf - btt;
    bmean = btt + t(Vtt%*%inv_f%*%t(cfe));
    bvar = Vtt - Vtt%*%inv_f%*%Vtt;
    bdraw[(t-i),] = rmvnorm(1,bmean,bvar); #bmean' + randn(1,m)*chol(bvar);
  }
  bdraw = t(bdraw);
  return(list(bdraw=bdraw,log_lik=log_lik))
}

# TVP-VAR Time varying structural VAR with constant covariance matrix
# ------------------------------------------------------------------------------------
# This code implements the Hierarchical Homoskedastic TVP-VAR
# ************************************************************************************
# The model is:
#
#     Y(t) = B0(t) + beta(t) x Y(t-1) + e(t) 
# 
# with e(t) ~ N(0, Sigma).
# The state equation is
#
#     beta(t) = theta(t) x A_0 + u(t)
#     theta(t) = theta(t-1) + eta(t)
#
# with u(t) ~ N(0, Q) and eta(t) ~ N(0,R).
# ************************************************************************************
  
#----------------------------------LOAD DATA----------------------------------------
  
Yraw = as.matrix(read.table("/Users/user/Dropbox/18. Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/R/Yraw.dat"))
Y = as.matrix(Yraw)

# Number of observations and dimension of X and Y
T=dim(Y)[1];  # T is the number of observations in Y
M=dim(Y)[2];  # M is the number of series Y
  
# Number of factors & lags:
p = 1; # p is number of lags in the VAR part
# ===================================| VAR EQUATION |==============================
# Generate lagged Y matrix. This will be part of the X matrix
ylag = embed(Y,(p+1))[,-c(1:ncol(Y))]; # Y is [T x K]. ylag is [T x (nk)]
# Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T

K = p*(M^2); # K is the number of elements in the state vector
# Create Z_t = I kron X
Z = zeros((T-p)*M,K);
for (i in 1:(T-p)) {
  # No constant
  ztemp = NULL; #eye(M);
  for (j in 1:p) {
    xtemp = ylag[i,((j-1)*M+1):(j*M),drop=F];
    xtemp = kronecker(diag(M),xtemp);
    ztemp = rbind(ztemp, xtemp);  #ok<AGROW>
  }
  Z[((i-1)*M+1):(i*M),] = ztemp;
}

# Redefine FAVAR variables y
y = t(Y[(p+1):T,]);
# Time series observations
T=ncol(y);

#----------------------------PRELIMINARIES---------------------------------
# Set some Gibbs - related preliminaries
nburn = 100;   # Number of burn-in-draws
nrep = nburn+100;    # Number of replications
it_print = 100; # Print in the screen every "it_print"-th iteration
  
#========= PRIORS:
#-------- Now set prior means and variances (_prmean / _prvar)
# theta_0 ~ N(B_OLS, 4Var(B_OLS))
theta_0_prmean = zeros(K,1);
theta_0_prvar = 4*diag(K);

# Q ~ IW(I,K+1)
Q_prmean = diag(K);
Q_prvar = K+1;

# R ~ IW(I,K+1))
R_prmean = diag(K);
R_prvar = K+1;

# Sigma ~ IW(I,M+1)
Sigma_prmean = diag(M);
Sigma_prvar = M+1;
  
#========= INITIALIZE MATRICES:
# Specify covariance matrices for measurement and state equations
Q = 0.0001*diag(K);
Qinv = solve(Q);
R = Q;
Rinv = solve(R);
Sigma = 0.1*diag(M);
Sigmainv = solve(Sigma);
theta_true = zeros(T+p,K);
theta_t = theta_true
beta_t = zeros(T,K);
A_0 = diag(K);

#Storage space for posterior means        
beta_t_mean = zeros(nrep-nburn,T,K);
theta_t_mean = zeros(nrep-nburn,T,K);
Q_mean = zeros(nrep-nburn,K,K);
R_mean = zeros(nrep-nburn,K,K);
Sigma_mean = zeros(nrep-nburn,M,M);
A_0_mean = zeros(nrep-nburn,K,K);

#----------------------------- END OF PRELIMINARIES ---------------------------
  
#==============================================================================================
#====================================== START SAMPLING ========================================
#==============================================================================================
print('Number of iterations');
for (irep in 1:nrep) {    # GIBBS iterations starts here
  # Print iterations
  if (irep%%it_print == 0) {
    print(paste0(round(irep/nrep*100,2),"%"))
  }
  #----------------------------------------------------------------------
  # Step 1: Sample beta_t ~ N(beta_post_mean,beta_post_var)
  #----------------------------------------------------------------------
  for (i in 1:T) {
    beta_post_var = solve(Qinv + t(Z[((i-1)*M+1):(i*M),])%*%Sigmainv%*%Z[((i-1)*M+1):(i*M),]);
    beta_post_mean = beta_post_var%*%(Qinv%*%(A_0%*%theta_t[i,]) + t(Z[((i-1)*M+1):(i*M),])%*%Sigmainv%*%y[,i]);
    beta_t[i,] = beta_post_mean + t(chol(beta_post_var))%*%rnorm(K);    
  }
  #----------------------------------------------------------------------
  # Step 2: Sample Q ~ iW(v_q,S_q)
  #----------------------------------------------------------------------
  sse_q = 0;
  for (i in 1:T) {
    sse_q =  sse_q + crossprod(t(matrix(beta_t[i,]-theta_t[i,]%*%A_0)))
  }
  v_q = T + Q_prvar;
  S_q = solve(Q_prmean + sse_q);
  Qinv = rWishart(1,v_q,S_q)[,,1];
  Q = solve(Qinv);
    
  #----------------------------------------------------------------------
  # Step 3: Sample theta_t using Carter and kohn
  #----------------------------------------------------------------------    
  kal1 = carter_kohn_hom2(t(beta_t),A_0,Q,R,K,K,T,theta_0_prmean,theta_0_prvar);
  theta_tc=kal1$bdraw
  log_lik=kal1$log_lik
  theta_t = t(theta_tc);
  
  #----------------------------------------------------------------------
  # Step 4: Sample R ~ iW(v_r,S_r)
  #----------------------------------------------------------------------
  sse_r = 0;
  theta_temp = theta_t[2:T,] - theta_t[1:(T-1),];
  for (i in 1:(T-1)) {
    sse_r =  sse_r + crossprod(t(matrix(theta_temp[i,])));
  }
  v_r = T + R_prvar;
  S_r = solve(R_prmean + sse_r);
  Rinv = rWishart(1,v_r,S_r)[,,1];
  R = solve(Rinv);
  
  #----------------------------------------------------------------------
  # Step 5: Sample A_0 ~ N(A_post_mean,A_post_var)
  #----------------------------------------------------------------------

    #For simplicity assume that A_0 is the identity matrix. Can you write
    #your own code to draw A_0 from the Normal distribution?

  # ---------------------------------------------------------------------
  #   Step 3: Sample Sigma from M(Sigma|y,B_t) which is i-Wishart
  # ---------------------------------------------------------------------
  yhat = matrix(0,M,T);
  for (i in 1:T) {
    yhat[,i] = y[,i] - Z[((i-1)*M+1):(i*M),]%*%beta_t[i,];
  }

  sse_S = matrix(0,M,M);
  for (i in 1:T) {
    sse_S = sse_S + crossprod(t(yhat[,i]));
  }

  Sinv = solve(sse_S + Sigma_prmean);
  Sigmainv = rWishart(1, T+Sigma_prvar, Sinv)[,,1];
  Sigma = solve(Sigmainv);
  Sigmachol = chol(Sigma);

  #-----------------------SAVE AFTER-BURN-IN DRAWS-----------------------
  if (irep > nburn) {
    beta_t_mean[(irep-nburn),,] = beta_t;
    theta_t_mean[(irep-nburn),,] = theta_t;
    Q_mean[(irep-nburn),,] = Q;
    R_mean[(irep-nburn),,] = R;
    Sigma_mean[(irep-nburn),,] = Sigma;
    A_0_mean[(irep-nburn),,] = A_0;
  } # END saving after burn-in results
} # END main Gibbs loop (for irep = 1:nrep+nburn)

#=============================GIBBS SAMPLER ENDS HERE==================================
Beta_t = apply(beta_t_mean,2:3,mean)
Beta_10 = apply(beta_t_mean,2:3,quantile,.10)
Beta_50 = apply(beta_t_mean,2:3,quantile,.50)
Beta_90 = apply(beta_t_mean,2:3,quantile,.90)
Theta_t = apply(theta_t_mean,2:3,mean)
Q = apply(Q_mean,2:3,mean)
R = apply(R_mean,2:3,mean)
S = apply(Sigma_mean,2:3,mean)
A0 = apply(A_0_mean,2:3,mean)

BB = rbind(Beta_10,Beta_90)
par(mfcol = c(M,M), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:ncol(Beta_t)) {
   plot(Beta_t[,i],type="l",xlab="",ylab="",main=paste0(colnames(Yraw)[i],"-",colnames(Yraw)[j]),xaxs="i",las=1,ylim=c(min(BB),max(BB)),tck=0.02)
   lines(Beta_10[,i],lty=1,col="steelblue4")
   lines(Beta_50[,i],lty=2,col="steelblue4")
   lines(Beta_90[,i],lty=1,col="steelblue4")
   abline(h=0,lty=2)
}

### END
