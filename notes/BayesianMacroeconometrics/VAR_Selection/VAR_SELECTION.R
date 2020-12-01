library("MASS")
library("matlab")
library("mvtnorm")
library("R.matlab")
simvargp = function(T=NULL,N=NULL,L=NULL,PHI=NULL,PSI=NULL,const=TRUE) {
   #--------------------------------------------------------------------------
   #   PURPOSE:
   #      Get matrix of Y generated from a VAR model
   #--------------------------------------------------------------------------
   #   INPUTS:
   #     T     - Number of observations (rows of Y)
   #     N     - Number of series (columns of Y)
   #     L     - Number of lags
   #
   #   OUTPUT:
   #     y     - [T x N] matrix generated from VAR(L) model
   # -------------------------------------------------------------------------
   
   #-----------------------PRELIMINARIES--------------------

   if (is.null(T)) {
      T=100   #Number of time series observations (T)
   }
   if (is.null(N)) {
      N = 4;   #Number of cross-sectional observations (N)
   }
   if (is.null(L)) {
      L = 1;             #Lag order
   }
   if (is.null(PHI)) {
      PHI = diag(N)
   }
   if (is.null(PSI)) {
      PSI = diag(N)
      PSI[upper.tri(PSI)] = 0.5
   }
   sigma = solve(crossprod(t(PSI)));

   #----------------------GENERATE--------------------------
   # Set storage in memory for y
   # First L rows are created randomly and are used as 
   # starting (initial) values 
   y = rbind(runif(L*N), zeros(T,N));
            
   # Now generate Y from VAR (L,PHI,PSI)
   for (nn in (L+1):(T+L)) {    
      u = t(chol(sigma))%*%rnorm(N);
      y[nn,] = y[(nn-1),]%*%PHI + t(u);
   }
   return = list(y=y,PHI=PHI)
}

#--------------------------------------------------------------------------
# Variable selection in a VAR model
#--------------------------------------------------------------------------
# The VAR can be written as
# 
# Y = Z x GAMMA x BETA + e,  e ~ N(0,I x SIGMA)
#
# where GAMMA is a diagonal matrix of 0-1 restriction indices, BETA are the
# regression coefficients and SIGMA the covariance matrix. This code
# samples recursively BETA - GAMMA - SIGMA from their posteriors. For more
# information, check the manual
#--------------------------------------------------------------------------
# I use artificial data to demonstrate how restriction search can lead to
# efficient inference. I generate a 5 variable VAR with only T=50 time
# series observations, one lag and no constant. While the unrestricted OLS 
# estimates perform poorly, variable selection gives better estimates. This
# is because we restrict irrelevant variables, then there are more degrees 
# of freedom to estimate the parameters on the relevant variables. Run the 
# code and see the final results printed in the command window.
#--------------------------------------------------------------------------

# Generate artificial data
var = simvargp()
PHI = var$PHI
Y = var$y

Y = as.matrix(read.table("/Users/user/Dropbox/Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/MATLAB/BVAR_Gibbs_DONE/Yraw.dat"))
y = Y
#--------------------------DATA HANDLING-----------------------------------
Traw = nrow(y)
m = ncol(y);
for (i in 1:m) {
   y[,i] = Y[,i]-mean(Y[,i])
}
plag=1; #choose number of lags

# Generate lagged Y matrix. This will be part of the X matrix
n=(plag*m)*m; 
p=(plag*m);
x2 = embed(y,plag+1)[,-c(1:m)];
x = kronecker(eye(m),x2);

# Form y matrix accordingly
# Delete first "L" rows to match the dimensions of X matrix
y2 = y[(plag+1):Traw,]; # This is the final Y matrix used for the VAR
y = c(y2);
#Traw was the dimesnion of initial data. T is the number of actual 
#time series observations of ymat and xmat (final Y & X)
T=Traw-plag;
head(y)
# ----------------Gibbs related preliminaries
nsave = 1000;# Number of draws to save
nburn = 1000;  # Number of draws to discard
ntot = nsave + nburn;  # Number of total draws

beta_draws = zeros(nsave,n);
gamma_draws = zeros(nsave,n);
sigma_draws = zeros(nsave,m,m);

# ----------------Set priors
# beta ~ N(b0,D0)
b0 = 0*ones(n,1);
D0 = 9*eye(n);
D0_inv = solve(D0);

# gamma_j ~ Bernoulli(1,p_j), for j = 1,...,p
p_j = 0.5*ones(n,1);

# sigma ~ iWishart(a,eta)
a = m;
eta = 1*eye(m);
eta_inv = solve(eta);

# Initialize parameters
gamma = ones(n,1);
beta_OLS = solve(crossprod(x))%*%crossprod(x,y);
beta_OLS2 = solve(crossprod(x2))%*%crossprod(x2,y2);
sse = crossprod(y2-x2%*%beta_OLS2);
sigma = sse/(T-(p-1));

#-----------Gibbs iterations
print('Running variable selection for VARs')
print('Iterations')
tic();
for (irep in 1:ntot) {
   if (irep%%500==0) {
      print(paste0(irep/ntot*100,"%"))
      toc();
   }

   V = kronecker(solve(sigma),eye(T));
   x_star = ones(T*m,1)%*%t(gamma)*x;
   
   # ------Step 1:
   # update BETA ~ Normal
   beta_var = solve(D0_inv + t(x_star)%*%V%*%x_star);
   beta_mean = beta_var%*%(D0%*%b0 + t(x_star)%*%V%*%y);
   beta = beta_mean + t(chol(beta_var))%*%rnorm(n);

   beta_mat = zeros(p,m); # beta_mat is the matrix of regr. coefficients
   for (i in 1:m) {
      beta_mat[,i] = beta[((i-1)*p+1):(i*p),];
   }
   
   # ------Step 2: 
   # It is preferable to draw each gamma(j) in random order
   ind_perm = sample(1:n)
   # update GAMMA ~ Bernoulli
   for (kk in 1:n) {
      j = ind_perm[kk];
      theta = beta*gamma;
      theta_j_star = theta;
      theta_j_star_star = theta;
      theta_j_star[j] = beta[j];
      theta_j_star_star[j] = 0;
      c_j = p_j[j]*exp(-0.5*t(y-x%*%theta_j_star)%*%V%*%(y-x%*%theta_j_star));
      d_j = (1-p_j[j])*exp(-0.5*t(y - x%*%theta_j_star_star)%*%V%*%(y-x%*%theta_j_star_star));
      p_j_tilde = c_j/(c_j+d_j);
      gamma[j] = rbinom(1,1,p_j_tilde);
   }

   gamma_mat = zeros(p,m); # gamma_mat concatenates the restriction indices in matrix form
   for (i in 1:m) {
      gamma_mat[,i] = gamma[((i-1)*p+1):(i*p),];
   }
   
   # ------Step 3:  
   theta = beta_mat*gamma_mat;
   # update SIGMA ~ iWishart
   R_1 = ginv(eta_inv + crossprod(y2-x2%*%theta));
   R_2 = (a + T);
   rd = rWishart(1,R_2,R_1)[,,1];
   sigma = solve(rd);
 
   # Save draws
   if (irep > nburn) {
      beta_draws[(irep-nburn),] = beta;
      gamma_draws[(irep-nburn),] = gamma;
      sigma_draws[(irep-nburn),,] = sigma;
   }
}

# The first column is the mean of the restriction indices,
# the second gives the means of the regression coefficients B,
# the third column gives the OLS estimate of B,
# the fourth column gives the true (generated) value of B
matrix(colMeans(gamma_draws),m)
matrix(colMeans(beta_draws),m)
matrix(beta_OLS,m)
PHI

### END
