library("psych") # for trace (tr)

carter_kohn = function(y,z,Ht,Qt,m,p,t,B0,V0) {
   # Carter and Kohn (1994), On Gibbs sampling for state space models.
   # Kalman Filter
   bp = B0;
   Vp = V0;
   bt = zeros(t,m);
   Vt = zeros(m^2,t);
   log_lik = 0;
   for (i in 1:t) {
      R = Ht[((i-1)*p+1):(i*p),];
      H = z[((i-1)*p+1):(i*p),];
      #F = eye(m);
      cfe = y[,i] - H%*%bp;   # conditional forecast error
      f = H%*%Vp%*%t(H) + R;    # variance of the conditional forecast error
      inv_f = solve(f);
      log_lik = log_lik + log(det(f)) + t(cfe)%*%inv_f%*%cfe;
      btt = bp + Vp%*%t(H)%*%inv_f%*%cfe;
      Vtt = Vp - Vp%*%t(H)%*%inv_f%*%H%*%Vp;
      if (i < t) {
         bp = btt;
         Vp = Vtt + Qt;
      }
      bt[i,] = t(bp);
      Vt[,i] = matrix(Vp,m^2,1);
   }

   # draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
   bdraw = zeros(t,m);
   bdraw[t,] = mvrnorm(1,btt,Vtt);
   
   # Backward recurssions
   for (i in 1:(t-1)) {
      bf = bdraw[(t-i+1),];
      btt = bt[(t-i),];
      Vtt = matrix(Vt[,(t-i)],m,m);
      f = Vtt + Qt;
      inv_f = solve(f);
      cfe = bf - btt;
      bmean = btt + Vtt%*%inv_f%*%cfe;
      bvar = Vtt - Vtt%*%inv_f%*%Vtt;
      bdraw[(t-i),] = mvrnorm(1,bmean,bvar); #t(bmean) + rnorm(m)%*%chol(bvar);
   }
   bdraw = t(bdraw);
   return = list(bdraw=bdraw,log_lik=log_lik)
}
tvpvarsim = function(T=NULL,N=NULL,L=NULL) {

   #TVPVARSIM Summary of this function goes here
   if (is.null(T)) {
      T=100   #Number of time series observations (T)
   }
   if (is.null(N)) {
      N = 4;   #Number of cross-sectional observations (N)
   }
   if (is.null(L)) {
      L = 1;             #Lag order
   }
   m=L*N*N;
   
   sigma_sd = matrix(c(1.0000,   -0.5000,  -0.2500,   -0.1250,
            -0.5000,    1.2500,   -0.3750,   -0.1875,
            -0.2500,   -0.3750,    1.3125,   -0.3437,
            -0.1250,   -0.1875,   -0.3437,    1.3281),ncol=N)

   b1 = matrix(c(0.7, 0, 0.35, 0, 0.75,0, 0.7, 0, 0, 0,
         0, 0.65, 0.7, 0, 0.75,0.4, 0, 0, 0.7, 0, 0, 0, 0, 0, 0.8),ncol=N);
   b1 = b1[1:N,1:N];
   b0 = c(b1);
   b = zeros(m,T+L);
   Qchol = diag(0.01*round(b0));
   for (i in 1:(T+L)) {  
      if (i==1) {            
         b[,i] = b0 + Qchol%*%rnorm(m);
      } else {
         b[,i] = b[,(i-1)] + Qchol%*%rnorm(m);
      }
   }

   y = rbind(runif(L*N), zeros(T,N))
   for (i in (L+1):(T+L)) {
      y[i,] = y[(i-1),]%*%matrix(b[1:(N*N),i],N,N) + rnorm(N)%*%chol(sigma_sd);
   }
   y=y[(L+1):T,];
   return = list(y=y,PHI=b)
}

# This code implements the model selection in a TVP-VAR.
#-----------------------------LOAD DATA------------------------------------
# Load Korobilis (2008) quarterly data
# load ydata.dat;
# load yearlab.dat;
# 
# Y = ydata - repmat(mean(ydata),size(ydata,1),1);

tvpvar = tvpvarsim()
#tvpvar = simvargp()
b = tvpvar$PHI
Y = tvpvar$y

#Y = as.matrix(read.table("/Users/user/Dropbox/Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/MATLAB/BVAR_Gibbs_DONE/Yraw.dat"))
#y = Y

# Number of observations and dimension of Y
t=nrow(Y);
M=ncol(Y);
p = M; # p is the dimensionality of Y
plag = 1; # plag is number of lags in the VAR part
numa = p*(p-1)/2;
# ===================================| VAR EQUATION |=========================
# Generate lagged Y matrix. This will be part of the X matrix
# Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
ylag = embed(Y,plag+1)[,-c(1:M)]; # Y is [T x m]. ylag is [T x (nk)]

m = plag*(p^2); # m is the number of elements in the state vector
# Create X_t matrix as in Primiceri equation (4). Of course I have reserved
# "X" for the data I use to extract factors, hence name this matrix "Z".
Z = zeros((t-plag)*p,m);
for (i in 1:(t-plag)) {
   ztemp = NULL; #eye(p);
   for (j in 1:plag) {       
      xtemp = ylag[i,((j-1)*p+1):(j*p)];
      xtemp = kronecker(eye(p),xtemp);
      ztemp = cbind(ztemp, xtemp);  #ok<AGROW>
   }
   Z[((i-1)*p+1):(i*p),] = ztemp;
}

# Redefine VAR variables y
y = t(Y[(plag+1):t,]);

# Time series observations
t=ncol(y);
#----------------------------PRELIMINARIES---------------------------------
# Set some Gibbs - related preliminaries
nrep = 50;  # Number of replications
nburn = 50;   # Number of burn-in-draws
nthin = 1;   # Consider every thin-th draw (thin value)
it_print = 100;  #Print in the screen every "it_print"-th iteration

#========= PRIORS:
# Use uninformative values
B_OLS = zeros(m,1);
VB_OLS = eye(m);

# Set some hyperparameters here
k_Q = 0.01;

# We need the sizes of some matrices as prior hyperparameters
sizeQ = m; # Size of matrix Q

#-------- Now set prior means and variances (_prmean / _prvar)
# B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean = B_OLS;
B_0_prvar = 4*VB_OLS;

# Note that for IW distribution I keep the _prmean/_prvar notation,
# but these are scale and shape parameters...
# Q ~ IW(k2_Q%*%size(subsample)%*%Var(B_OLS),size(subsample))
Q_prmean = ((k_Q)^2)*(1 + sizeQ)*VB_OLS;
Q_prvar = 1 + sizeQ;

# gamma_j ~ Bernoulli(1,p_j), for j = 1,...,m
p_j = 0.5*ones(m,1);

# sigma ~ iWishart(a,eta)
a = 0;
eta = 1*eye(p);
eta_inv = 0*solve(eta);

#========= INITIALIZE MATRICES:
# Specify covariance matrices for measurement and state equations%*%/
consQ = 0.0001;
consH = 0.0001;
Qdraw = consQ*eye(m);
Qchol = sqrt(consQ)*eye(m);
sigma_sd = matrix(c(1.0000, -0.5000, -0.2500, -0.1250, 
                    -0.5000, 1.2500, -0.3750, -0.1875,
                    -0.2500, -0.3750, 1.3125, -0.3437, 
                    -0.1250, -0.1875, -0.3437, 1.3281),4);
sigma_sd = matrix(c(1.0000, -0.5000, -0.2500, 
                    -0.5000, 1.2500, -0.3750,
                    -0.2500, -0.3750, 1.3125),M);
Ht = kronecker(ones(t,1),sigma_sd);
Htsd = kronecker(ones(t,1),chol(sigma_sd));
Btdraw = zeros(m,t);
gamma = ones(m,1);

# Storage matrices for posteriors and stuff
Bt_postmean = zeros(m,t);
Qmean = zeros(m,m);
gamma_draws = zeros(nrep,m);
sigma_draws = zeros(nrep,p,p);
#----------------------------- END OF PRELIMINARIES ---------------------------

#====================================== START SAMPLING ========================================
#==============================================================================================
tic(); # This is just a timer
print('Number of iterations');
for (irep in 1:(nrep+nburn)) { # GIBBS iterations starts here
   # Print iterations
   if (irep%%it_print==0) {
      print(paste0(irep/it_print*100,"%"));
      toc();
   }
   
   # -----------------------------------------------------------------------------------------
   #   STEP I: Sample B from Normal
   # -----------------------------------------------------------------------------------------
   ck = carter_kohn(y,Z%*%diag(c(gamma)),Ht,Qdraw,m,p,t,B_0_prmean,B_0_prvar);
   Btdrawc=ck$bdraw
   log_lik=ck$log_lik
   Btdraw = Btdrawc;

   Btemp = t(Btdraw[,2:t])-t(Btdraw[,1:(t-1)]);
   sse_2 = zeros(m,m);
   for (i in 1:(t-1)) {
      sse_2 = sse_2 + crossprod(t(Btemp[i,]));
   }

   Qinv = solve(sse_2 + Q_prmean);
   Qinvdraw = rWishart(1,t+Q_prvar,Qinv)[,,1];
   Qdraw = solve(Qinvdraw);
   Qchol = chol(Qdraw);

   #-----------------------------------------------------------------------------------------
   #   STEP II: Sample gamma from Bernoulli
   #-----------------------------------------------------------------------------------------
   for (j in 1:m) {
      theta = diag(gamma)*Btdraw;
      theta_j_star = theta;
      theta_j_star_star = theta;
      theta_j_star[j,] = Btdraw[j,];
      theta_j_star_star[j,] = 0;
      x_star_star=x_star=NULL;
      for (i in 1:t) {
         xtemp1 = solve(Htsd[((i-1)*p+1):(i*p),])%*%(y[,i] - Z[((i-1)*p+1):(i*p),]%*%theta_j_star[,i]);       
         xtemp2 = solve(Htsd[((i-1)*p+1):(i*p),])%*%(y[,i] - Z[((i-1)*p+1):(i*p),]%*%theta_j_star_star[,i]);
         x_star = cbind(x_star, xtemp1); #ok<AGROW>
         x_star_star = cbind(x_star_star, xtemp2); #ok<AGROW>
      }
      c_j = p_j[j]*exp(-1/2*tr(crossprod(x_star)));
      d_j = (1-p_j[j])*exp(-1/2*tr(crossprod(x_star_star)));
      p_j_tilde = c_j/(c_j+d_j);
      gamma[j,1] = rbinom(1,1,p_j_tilde);
   }

   #-----------------------------------------------------------------------------------------
   #   STEP II: Sample SIGMA from inverse-Wishart
   #-----------------------------------------------------------------------------------------        
   sse_2=zeros(p,p);
   for (i in 1:t) {
      theta = Btdraw[,1]*gamma;
      sse_2 = sse_2 + crossprod(t(y[,i]-Z[((i-1)*p+1):(i*p),]%*%theta));
   }
   # Draw SIGMA
   R_1 = solve(eta_inv + sse_2);
   R_2 = (a+t);
   rd = rWishart(1,R_2,R_1)[,,1];
   sigma=solve(rd);

   Ht = kronecker(ones(t,1),sigma);
   Htsd = kronecker(ones(t,1),chol(sigma));

   #----------------------------SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES -----------------
   if (irep > nburn) {
      Bt_postmean = Bt_postmean + Btdraw;
      gamma_draws[(irep-nburn),] = gamma;
      Qmean = Qmean + Qdraw;
      sigma_draws[(irep-nburn),,] = sigma;
   } # END saving after burn-in results 
} #END main Gibbs loop (for irep = 1:nrep+nburn)
toc(); # Stop timer and print total time

#=============================GIBBS SAMPLER ENDS HERE==================================
Bt_postmean = Bt_postmean/nrep;           # mean of time-varying parameters B(t)
Bt_postmean
Qmean = Qmean/nrep;                       # mean of the matrix Q, which is the covariance of B(t) 
Qmean
Qmeangamma_mean = colMeans(gamma_draws);  # average restriction indices
gamma_mean
SIGMA_mean = apply(sigma_draws,2:3,mean)  # mean of VAR covariance matrix
SIGMA_mean

### END
