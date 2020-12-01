
# BVAR_FULL.m 
library("MASS")
library("matlab")
library("mvtnorm")
library("R.matlab")
# This code replicates the results from the 1st empirical illustration 
# in Koop and Korobilis (2009).
#
# You can chose 6 different priors. For some priors, analytical results are
# available, so Monte Carlo Integration is used. For other priors, you need
# to use the Gibbs sampler. For Gibbs sampler models I take a number of
# 'burn-in draws', so that I keep only the draws which have converged.
#
# The specification of the prior hyperparmeters are in the file
# prior_hyper.m. See there for details.
#
# The convention used here is that ALPHA is the K x M matrix of VAR coefficients,
# alpha is the KM x 1 column vector of vectorized VAR coefficients, i.e.
# alpha = vec(ALPHA), and SIGMA is the M x M VAR covariance matrix.
#--------------------------------------------------------------------------
# Bayesian estimation, prediction and impulse response analysis in VAR
# models using posterior simulation. Dependent on your choice of forecasting,
# the VAR model is:
#
# In this code we provide direct (as opposed to iterated) forecasts
# Direct h-step ahead foreacsts:
#  Y(t+h) = A0 + Y(t) x A1 + ... + Y(t-p+1) x Ap + e(t+h)
#
# so that in this case there are also p lags of Y (from 0 to p-1).
#
# In any of the two cases, the model is written as:
#
# Y(t) = X(t) x A + e(t)
#
# where e(t) ~ N(0,SIGMA), and A summarizes all parameters. Note that we
# also use the vector a which is defined as a=vec(A).
#--------------------------------------------------------------------------
# NOTES: The code sacrifices efficiency for clarity. It follows the
#        theoretical equations in the monograph and the manual.
#
# AUTHORS: Gary Koop and Dimitris Korobilis
# CONTACT: dikorombilis@yahoo.gr
#--------------------------------------------------------------------------

#------------------------------LOAD DATA-----------------------------------
# Load Quarterly US data on inflation, unemployment and interest rate, 
# 1953:Q1 - 2006:Q3
Yraw = as.matrix(read.table("/Users/user/Dropbox/Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/R/Yraw.dat"))
#Yraw = Yraw[29:189,];
# or Simulate data from a simple VAR Data Generating process
#Yraw = bvardgp();

# In any case, name the data you load 'Yraw', in order to avoid changing the
# rest of the code. Note that 'Yraw' is a matrix with T rows by M columns,
# where T is the number of time series observations (usually months or
# quarters), while M is the number of VAR dependent macro variables.

#----------------------------PRELIMINARIES---------------------------------
# Define specification of the VAR model
constant = 1;     # 1: if you desire intercepts, 0: otherwise 
p = 1;            # Number of lags on dependent variables
forecasting = 1;  # 1: Compute h-step ahead predictions, 0: no prediction
repfor = 50;   # Number of times to obtain a draw from the predictive 
               # density, for each generated draw of the parameters
h = 2;         # Number of forecast periods
impulses = 1;  # 1: compute impulse responses, 0: no impulse responses
ihor = 30;     # Horizon to compute impulse responses

# Set prior for BVAR model:
prior = 1
# prior = 1 --> Diffuse ('Jeffreys') (M-C Integration)
# prior = 2 --> Minnesota(M-C Integration)
# prior = 3 --> Normal-Wishart (M-C Integration)  
# prior = 4 --> Independent Normal-Wishart(Gibbs sampler)
# prior = 5 --> SSVS in mean-Wishart(Gibbs sampler)
# prior = 6 --> SSVS in mean-SSVS in covariance (Gibbs sampler)

# Gibbs-related preliminaries
nsave = 5000; # Final number of draws to save
nburn = 5000; # Draws to discard (burn-in)
              # For models using analytical results, there are no convergence issues (You
              # are not adviced to change the next 3 lines)

if (prior ==1 || prior == 2 || prior == 3) {
   nburn = 0*nburn;
}

ntot = nsave + nburn;  # Total number of draws
it_print = 500; # Print on the screen every "it_print"-th iteration
#--------------------------DATA HANDLING-----------------------------------
# Get initial dimensions of dependent variable
Traw = nrow(Yraw);
M = ncol(Yraw);

# The model specification is different when implementing direct forecasts,
# compared to the specification when computing iterated forecasts.
if (forecasting==1) {
   if (h<=0) { # Check for wrong (incompatible) user input
      print('You have set forecasting, but the forecast horizon choice is wrong')
   }
   # Now create VAR specification according to forecast method

   Y1 = Yraw[(h+1):nrow(Yraw),];
   Y2 = Yraw[2:(nrow(Yraw)-h),];
   Traw = Traw - h - 1;

} else {
   Y1 = Yraw;
   Y2 = Yraw;
}

# Generate lagged Y matrix. This will be part of the X matrix
Ylag = embed(Y2,p+1)[,-c(1:M)]; # Y is [T x M]. ylag is [T x (Mp)]
# Now define matrix X which has all the R.H.S. variables (constant, lags of
# the dependent variable and exogenous regressors/dummies).
# Note that in this example I do not include exogenous variables (other macro
# variables, dummies, or trends). You can load a file with exogenous
# variables, call them, say W, and then extend variable X1 in line 133, as:
# X1 = cbind(ones(Traw-p,1), Ylag[(p+1):Traw,], W[(p+1):Traw,]); and line 135 as:
# X1 = [Ylag(p+1:Traw,:)  W(p+1:Traw,:)];
if (constant) {
   X1 = cbind(ones(Traw-p,1), Ylag);
} else {
   X1 = Ylag;  #ok<UNRCH>
}
# Get size of final matrix X
Traw3 = nrow(X1)
K = ncol(X1);

# Create the block diagonal matrix Z
Z1 = kronecker(eye(M),X1);

# Form Y matrix accordingly
# Delete first "LAGS" rows to match the dimensions of X matrix
Y1 = Y1[(p+1):Traw,]; # This is the final Y matrix used for the VAR

# Traw was the dimesnion of the initial data. T is the number of actual 
# time series observations of Y and X
T = Traw - p;

#========= FORECASTING SET-UP:
# Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
if (forecasting==1) {
   Y_pred = zeros(nsave*repfor,M); # Matrix to save prediction draws
   PL = zeros(nsave,1);# Matrix to save Predictive Likelihood
   
   # Direct forecasts, we only need to keep the last observation for evaluation
   Y = Y1[1:(nrow(Y1)-1),]; 
   X = X1[1:(nrow(X1)-1),];
   Z = kronecker(eye(M),X);
   T = T - 1;

} else { # if no prediction is present, keep all observations
   Y = Y1;
   X = X1;
   Z = Z1;
}

#========= IMPULSE RESPONSES SET-UP:
# Create matrices to store forecasts
if (impulses == 1) {
   imp = zeros(nsave,M,M,ihor);
   bigj = zeros(M,M*p);
   bigj[1:M,1:M] = eye(M);
}

#-----------------------------PRELIMINARIES--------------------------------
# First get ML estimators
A_OLS = solve(crossprod(X))%*%crossprod(X,Y); # This is the matrix of regression coefficients
a_OLS = c(A_OLS);# This is the vector of parameters, i.e. it holds
# that a_OLS = vec(A_OLS)
SSE = crossprod(Y-X%*%A_OLS);# Sum of squared errors
SIGMA_OLS = SSE/(T-K+1);

# Initialize Bayesian posterior parameters using OLS values
alpha = a_OLS;  # This is the single draw from the posterior of alpha
ALPHA = A_OLS;  # This is the single draw from the posterior of ALPHA
SSE_Gibbs = SSE;# This is the SSE based on each draw of ALPHA
SIGMA = SIGMA_OLS; # This is the single draw from the posterior of SIGMA
IXY =  kronecker(eye(M),crossprod(X,Y));

# Storage space for posterior draws
alpha_draws = zeros(nsave,K*M);# save draws of alpha
ALPHA_draws = zeros(nsave,K,M);# save draws of alpha
SIGMA_draws = zeros(nsave,M,M);# save draws of ALPHA


#-----------------Prior hyperparameters for bvar model
# load file which sets hyperparameters for chosen prior
# Define hyperparameters under different priors
# I indicate the exact line in which the prior hyperparameters can be found

if (prior == 1) { # Diffuse
   # I guess there is nothing to specify in this case!
   # Posteriors depend on OLS quantities
} else if (prior == 2) { # Minnesota-Whishart
   # Prior mean on VAR regression coefficients
   A_prior = rbind(zeros(1,M), 0.9*eye(M), matrix(0,(p-1)*M,M));  #<---- prior mean of ALPHA (parameter matrix) 
   a_prior = c(A_prior);#<---- prior mean of alpha (parameter vector)
   
   # Minnesota Variance on VAR regression coefficients
   # First define the hyperparameters 'a_bar_i'
   a_bar_1 = 0.5;
   a_bar_2 = 0.5;
   a_bar_3 = 10^2;
   
   # Now get residual variances of univariate p_MIN-lag autoregressions. Here
   # we just run the AR(p) model on each equation, ignoring the constant
   # and exogenous variables (if they have been specified for the original VAR model)
   p_MIN = 6;
   sigma_sq = zeros(M,1); # vector to store residual variances
   for (i in 1:M) {
      # Create lags of dependent variables   
      Ylag_i = rbind(matrix(rep(0,p^2),nrow=p),matrix(embed(Yraw[,i],p+1)[,-1],ncol=p));
      Ylag_i = Ylag_i[(p_MIN+1):Traw,];
      
      X_i = cbind(ones(Traw-p_MIN,1), Ylag_i);
      Y_i = Yraw[(p_MIN+1):Traw,i];
      # OLS estimates of i-th equation
      alpha_i = solve(crossprod(X_i))%*%crossprod(X_i,Y_i);
      sigma_sq[i,1] = (1/(Traw-p_MIN))*crossprod(Y_i-X_i%*%alpha_i);
   }
   # Now define prior hyperparameters.
   # Create an array of dimensions K x M, which will contain the K diagonal
   # elements of the covariance matrix, in each of the M equations.
   V_i = zeros(K,M);
   
   # index in each equation which are the own lags
   ind = zeros(M,p);
   for (i in 1:M) {
      ind[i,] = seq(i+constant,K,by=M);
   }
   for (i in 1:M) {  # for each i-th equation
      for (j in 1:K) {  # for each j-th RHS variable
         if (constant) {
            if (j==1) { # if there is constant, use this code
               V_i[j,i] = a_bar_3*sigma_sq[i,1]; # variance on constant                
            } else if (length(which(j==ind[i,]))>0) {
               V_i[j,i] = a_bar_1/(ceiling((j-1)/M)^2); # variance on own lags           
               # Note: the "ceil((j-1)/M)" command finds the associated lag  number for each parameter
            } else {
               for (kj in 1:M) {
                  if (length(which(j==ind[kj,]))>0) {
                     ll = kj;                   
                  }
               } # variance on other lags  
               V_i[j,i] = (a_bar_2*sigma_sq[i,1])/((ceiling((j-1)/M)^2)*sigma_sq[ll,1]);           
            }
         } else { # if no constant is defined, then use this code
            if (length(which(j==ind[i,]))>0) {
               V_i[j,i] = a_bar_1/(ceiling(j/M)^2); # variance on own lags
            } else {
               for (kj in 1:M) {
                  if (length(which(j==ind[kj,]))>0) {
                     ll = kj;
                  }
               } # variance on other lags  
               V_i[j,i] = (a_bar_2*sigma_sq[i,1])/((ceiling(j/M)^2)*sigma_sq[ll,1]);            
            }
         }
      }
   }
   
   # Now V is a diagonal matrix with diagonal elements the V_i
   V_prior = diag(c(V_i));  # this is the prior variance of the vector alpha
   
   #NOTE: No prior for SIGMA. SIGMA is simply a diagonal matrix with each
   #diagonal element equal to sigma_sq(i). See Kadiyala and Karlsson (1997)
   SIGMA = diag(c(sigma_sq));
   
} else if (prior == 3) { # Normal-Whishart
   # Hyperparameters on a ~ N(a_prior, SIGMA x V_prior)
   A_prior = 0*ones(K,M);   #<---- prior mean of ALPHA (parameter matrix)
   a_prior = c(A_prior);    #<---- prior mean of alpha (parameter vector)
   V_prior = 10*eye(K);     #<---- prior variance of alpha
   
   # Hyperparameters on solve(SIGMA) ~ W(v_prior,solve(S_prior))
   v_prior = M+1;          #<---- prior Degrees of Freedom (DoF) of SIGMA
   S_prior = eye(M);       #<---- prior scale of SIGMA
   inv_S_prior = solve(S_prior);    
   
} else if (prior == 4) {  # Independent Normal-Wishart
   n = K*M; # Total number of parameters (size of vector alpha)
   a_prior = 0*ones(n,1);   #<---- prior mean of alpha (parameter vector)
   V_prior = 10*eye(n);     #<---- prior variance of alpha
   
   # Hyperparameters on solve(SIGMA) ~ W(v_prior,solve(S_prior))
   v_prior = M+1;           #<---- prior Degrees of Freedom (DoF) of SIGMA
   S_prior = eye(M);        #<---- prior scale of SIGMA
   inv_S_prior = solve(S_prior);
   
} else if (prior == 5 || prior == 6) { # SSVS on alpha, Wishart or SSVS on SIGMA    
   n = K*M; # Total number of parameters (size of vector alpha)
   # mean of alpha
   a_prior = zeros(n,1);
   # This is the std of the OLS estimate ofalpha. You can use this to 
   # scale tau_0 and tau_1 (see below) if you want.
   sigma_alpha = sqrt(diag(kronecker(SIGMA,solve(crossprod(X)))));
   # otherwise, set ' sigma_alpha = ones(n,1); '
   # SSVS variances for alpha
   tau_0 = 0.1*sigma_alpha;   # Set tau_[0i], tau_[1i]
   tau_1 = 10*sigma_alpha;
   # Priors on SIGMA
   if (prior == 6) { #SSVS on SIGMA
      # SSVS variances for non-diagonal elements of SIGMA     
      kappa_0 = 0.1; # Set kappa_[0ij], kappa_[1ij]
      kappa_1 = 6;
      # Hyperparameters for diagonal elements of SIGMA (Gamma density)
      a_i = .01;
      b_i = .01;
   } else if (prior == 5) { #Wishart on SIGMA
      # Hyperparameters on solve(SIGMA) ~ W(v_prior,solve(S_prior))
      v_prior = M+1;             #<---- prior Degrees of Freedom (DoF) of SIGMA
      S_prior = eye(M);          #<---- prior scale of SIGMA  
      inv_S_prior = solve(S_prior);
   }
   
   # Hyperparameters for Gamma ~ BERNOULLI(m,p_i), see eq. (14)
   p_i = .5;
   # Hyperparameters for Omega_[j] ~ BERNOULLI(j,q_ij), see eq. (17)
   q_ij = .5;
   
   # Initialize Gamma and Omega vectors
   gammas = ones(n,1);       # vector of Gamma  
   omega=cell(1,M-1);
   for (kk_1 in 1:(M-1)) {
      omega[[kk_1]] = ones(kk_1,1);	# Omega_j
   }
   
   # Set space in memory for some vectors that we are using in SSVS
   gamma_draws = zeros(nsave,n); # vector of gamma draws
   omega_draws = zeros(nsave,.5*M*(M-1)); # vector of omega draws    
} else {
   print('Wrong choice of prior, prior = 1-6 only')
}
#-------------------- Prior specification ends here

#========================== Start Sampling ================================
#==========================================================================
tic();
print('Number of iterations');
for (irep in 1:ntot) { #Start the Gibbs "loop"
   if (irep%%it_print==0) { # print iterations
      print(paste0(irep/ntot*100,"%"));
      toc();
   }
   
   #--------- Draw ALPHA and SIGMA with Diffuse Prior
   if (prior == 1) {
      # Posterior of alpha|SIGMA,Data ~ Normal
      V_post = kronecker(SIGMA,solve(crossprod(X)));
      alpha = a_OLS + t(chol(V_post))%*%rnorm(K*M);# Draw alpha
      ALPHA = matrix(alpha,K,M); # Create draw of ALPHA 
   
      # Posterior of SIGMA|Data ~ iW(SSE_Gibbs,T-K) 
      SIGMA = solve(rWishart(1,T-K,solve(SSE_Gibbs))[,,1]);# Draw SIGMA

      #--------- Draw ALPHA and SIGMA with Minnesota Prior
   } else if (prior == 2) {
      #Draw ALPHA
      for (i in 1:M) {
         inv_V_prior = solve(V_prior[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)])
         V_post = solve( inv_V_prior + 1/SIGMA[i,i]*crossprod(X));
         a_post = V_post%*%(inv_V_prior%*%a_prior[((i-1)*K+1):(i*K)] + 1/SIGMA[i,i]*t(X)%*%Y[,i]);
         alpha[((i-1)*K+1):(i*K)] = a_post + chol(V_post)%*%rnorm(K); # Draw alpha
      }
      ALPHA = matrix(alpha,K,M); # Create draw in terms of ALPHA

      # SIGMA in this case is a known matrix, whose form is decided in
      # the prior (see prior_hyper.m)

      #--------- Draw ALPHA and SIGMA with Normal-Wishart Prior
   } else if (prior == 3) {
      # ******Get all the required quantities for the posteriors 
      V_post = solve(solve(V_prior) + crossprod(X));
      A_post = V_post%*%(solve(V_prior)%*%A_prior + crossprod(X)%*%A_OLS);
      a_post = c(A_post);
 
      S_post = SSE + S_prior + t(A_OLS)%*%crossprod(X)%*%A_OLS + t(A_prior)%*%solve(V_prior)%*%A_prior - t(A_post)%*%(solve(V_prior) + crossprod(X))%*%A_post;
      v_post = T + v_prior;
  
      # This is the covariance for the posterior density of alpha
      COV = kronecker(SIGMA,V_post);
  
      # Posterior of alpha|SIGMA,Data ~ Normal
      alpha = a_post + t(chol(COV))%*%rnorm(K*M);  # Draw alpha
      ALPHA = reshape(alpha,K,M); # Draw of ALPHA
  
      # Posterior of SIGMA|ALPHA,Data ~ iW(solve(S_post),v_post)
      SIGMA = solve(rWishart(1,v_post,solve(S_post))[,,1]);# Draw SIGMA
  
      #--------- Draw ALPHA and SIGMA with Independent Normal-Wishart Prior
   } else if (prior == 4) {
      VARIANCE = kronecker(solve(SIGMA),eye(T));
      V_post = solve(V_prior + t(Z)%*%VARIANCE%*%Z);
      a_post = V_post%*%(V_prior%*%a_prior + t(Z)%*%VARIANCE%*%c(Y));
      alpha = a_post + chol(V_post)%*%rnorm(n); # Draw of alpha
      ALPHA = matrix(alpha,K,M); # Draw of ALPHA
  
      # Posterior of SIGMA|ALPHA,Data ~ iW(solve(S_post),v_post)
      v_post = T + v_prior;
      S_post = S_prior + crossprod(Y - X%*%ALPHA);
      SIGMA = solve(rWishart(1,v_post,solve(S_post))[,,1]);# Draw SIGMA
  
      #--------- Draw ALPHA and SIGMA using SSVS prior 
   } else if (prior == 5 || prior == 6) {
      # Draw SIGMA
      if (prior == 5) { # Wishart
         # Posterior of SIGMA|ALPHA,Data ~ iW(solve(S_post),v_post)
         v_post = T + v_prior;
         S_post = inv_S_prior + crossprod(Y - X%*%ALPHA);
         SIGMA = solve(rWishart(1,v_post,solve(S_post))[,,1]);# Draw SIGMA
      } else if (prior == 6) { # SSVS
         # Draw psi|alpha,gamma,omega,DATA from the GAMMA dist.
         # Get S_[j] - upper-left [j x j] submatrices of SSE
         # The following loop creates a cell array with elements S_1,
         # S_2,...,S_j with respective dimensions 1x1, 2x2,...,jxj
         S=cell(1,M);
         for (kk_2 in 1:M) {
            S[[kk_2]] = SSE_Gibbs[1:kk_2,1:kk_2];
         }
         # Set also SSE =(s_[i,j]) & get vectors s_[j]=(s_[1,j] , ... , s_[j-1,j])
         s=cell(1,M-1);
         for (kk_3 in 2:M) {
            s[[kk_3-1]] = SSE_Gibbs[1:(kk_3 - 1),kk_3];
         }
         # Parameters for Heta|omega ~ N_[j-1](0,D_[j]*R_[j]*D_[j]), see eq. (15)
         # Create and update h_[j] matrix
         # If omega_[ij] = 0 => h_[ij] = kappa0, else...
         hh=cell(1,M-1);
         for (kk_4 in 1:(M-1)) {
            omeg = omega[[kk_4]];
            het = hh[[kk_4]];
            for (kkk in 1:nrow(omeg)) {  
               if (omeg[kkk,1] == 0) {
                  het[kkk] = kappa_0;
               } else {
                  het[kkk] = kappa_1;
               }
            }
            hh[[kk_4]] = het;
         }
         # D_j = diag(hh_[1j],...,hh_[j-1,j])
         D_j=cell(1,M-1);
         for (kk_5 in 1:(M-1)) {
            D_j[[kk_5]] = diag(length(hh[[kk_5]]))*hh[[kk_5]];
         }
         # Now create covariance matrix D_[j]*R_[j]*D_[j], see eq. (15)
         DD_j=cell(1,M-1);
         for (kk_6 in 1:(M-1)) {
            DD = D_j[[kk_6]];
            DD_j[[kk_6]] = DD%*%DD;
         }
         # Create B_[i] matrix
         B=cell(1,M);
         for (rr in 1:M) {  
            if (rr == 1) {
               B[[rr]] = b_i + 0.5*(SSE[rr,rr]);
            } else if (rr > 1) {
               s_i = s[[rr-1]];
               S_i = S[[rr-1]];
               DiDi = DD_j[[rr-1]];
               B[[rr]] = b_i + 0.5*(SSE_Gibbs[rr,rr] - t(s_i)%*%solve(S_i + solve(DiDi))%*%s_i);
            }
         }
         # Now get B_i from cell array B, and generate (psi_[ii])^2
         B_i = unlist(B);
         psi_ii_sq = zeros(M,1);
         for (kk_7 in 1:M) {	  
            psi_ii_sq[kk_7,1] = rgamma(1,(a_i + 0.5*T),B_i[kk_7]);
         }

         # Draw eta|psi,phi,gamma,omega,DATA from the [j-1]-variate NORMAL dist.
         eta = cell(1,M-1);
         for (kk_8 in 1:(M-1)) { 
            s_i = s[[kk_8]];
            S_i = S[[kk_8]];
            DiDi = DD_j[[kk_8]];
            miu_j = - sqrt(psi_ii_sq[kk_8+1])*(solve(S_i + solve(DiDi))%*%s_i);
            Delta_j = solve(S_i + solve(DiDi));
            eta[[kk_8]] = miu_j + t(chol(Delta_j))%*%rnorm(kk_8);
         }

         # Draw omega|eta,psi,phi,gamma,omega,DATA from BERNOULLI dist.
         omega_vec = NULL; #temporary vector to store draws of omega
         for (kk_9 in 1:(M-1)) { 
            omeg_g = omega[[kk_9]];
            eta_g = eta[[kk_9]];
            for (nn in 1:nrow(omeg_g)) {  # u_[ij1], u_[ij2], see eqs. (32 - 33)  
               u_ij1 = (1/kappa_0)*exp(-0.5*((eta_g[nn])^2)/((kappa_0)^2))*q_ij;
               u_ij2 = (1/kappa_1)*exp(-0.5*((eta_g[nn])^2)/((kappa_1)^2))*(1-q_ij);
               ost = u_ij1/(u_ij1 + u_ij2);
               omeg_g[nn,1] = rbinom(1,1,1-ost);
               omega_vec = rbind(omega_vec, omeg_g[nn,1]); #ok<AGROW>
            }
            omega[[kk_9]] = omeg_g; #ok<AGROW>
         }

         # Create PSI matrix from individual elements of "psi_ii_sq" and "eta"
         PSI_ALL = zeros(M,M);
         for (nn_1 in 1:M) {  # first diagonal elements
            PSI_ALL[nn_1,nn_1] = sqrt(psi_ii_sq[nn_1,1]);
         }
         for (nn_2 in 1:(M-1)) { # Now non-diagonal elements
            eta_gg = eta[[nn_2]];
            for (nnn in 1:nrow(eta_gg)) {
               PSI_ALL[nnn,(nn_2+1)] = eta_gg[nnn];
            }
         }
         # Create SIGMA
         SIGMA = solve(crossprod(t(PSI_ALL)));  
      } # END DRAWING SIGMA 

      # Draw alpha  
      # Hyperparameters for alpha|gamma ~ N_[m](0,D*D)
      h_i = zeros(n,1);# h_i is tau_0 if gamma=0 and tau_1 if gamma=1
      for (nn_3 in 1:n) {
         if (gammas[nn_3,1] == 0) {
            h_i[nn_3,1] = tau_0[nn_3];
         } else if (gammas[nn_3,1] == 1) {  
            h_i[nn_3,1] = tau_1[nn_3];
         }
      }
      D = c(h_i)*eye(n); # Create D. Here D=diag(h_i) will also do
      DD = D%*%D;# Prior covariance matrix for Phi_m
      isig = solve(SIGMA);
      psi_xx = kronecker(solve(SIGMA),crossprod(X));
      V_post = solve(psi_xx + solve(DD));

     visig=c(isig);
     a_post = V_post%*%(IXY%*%visig + solve(DD)%*%a_prior);
     alpha = a_post + t(chol(V_post))%*%rnorm(n); # Draw alpha
     alpha = a_post + t(chol(V_post))%*%rnorm(n); # Draw alpha
     ALPHA = matrix(alpha,K,M); # Draw of ALPHA
  
     # Draw gamma|phi,psi,eta,omega,DATA from BERNOULLI dist. 
     for (nn_6 in 1:n) {
        u_i1 = (1/tau_0[nn_6])*exp(-0.5*(alpha[nn_6]/(tau_0[nn_6]))^2)*p_i;  
        u_i2 = (1/tau_1[nn_6])*exp(-0.5*(alpha[nn_6]/(tau_1[nn_6]))^2)*(1-p_i);
        gst = u_i1/(u_i1 + u_i2);
        gammas[nn_6,1] = rbinom(1,1,1-gst); #ok<AGROW>
     }
      
     # Save new Sum of Squared Errors (SSE) based on draw of ALPHA  
     SSE_Gibbs = crossprod(Y-X%*%ALPHA);
  }
  # =============Estimation ends here
  
  # ****************************|Predictions, Responses, etc|***************************
  if (irep > nburn) {
      #=========FORECASTING:
      if (forecasting==1) {
         Y_temp = zeros(repfor,M);
         # compute 'repfor' predictions for each draw of ALPHA and SIGMA
         for (ii in 1:repfor) {
            xxx = ifelse(2<(M*(p-1)+1),X[T,2:(M*(p-1)+1)],NA)
            X_fore = c(1, Y[T,], xxx);
            X_fore = na.omit(X_fore)
            # Forecast of T+1 conditional on data at time T
            Y_temp[ii,] = X_fore%*%ALPHA + rnorm(M)%*%chol(SIGMA);
         }
         # Matrix of predictions
         Y_pred[(((irep-nburn)-1)*repfor+1):((irep-nburn)*repfor),] = Y_temp;
         # Predictive likelihood
         PL[(irep-nburn),] = dmvnorm(Y1[(T+1),],X[T,]%*%ALPHA,SIGMA);
         if (PL[(irep-nburn),] == 0) {
            PL[(irep-nburn),] = 1;
         }
      } # end forecasting
      #=========Forecasting ends here

      #=========IMPULSE RESPONSES:
      if (impulses==1) {
         Bv = zeros(M,M,p);
         for (i_1 in 1:p) {
            Bv[,,i_1] = ALPHA[(1+((i_1-1)*M + 1)):(i_1*M+1),];               
         }
         
         # st dev matrix for structural VAR
         shock = t(chol(SIGMA));
         d = diag(diag(shock));
         shock = solve(d)%*%shock;
         
         neq=dim(Bv)[1];
         nvar=dim(Bv)[2];
         nlag=dim(Bv)[3];
         
         responses=zeros(nvar,neq,ihor);
         responses[,,1]=t(shock); # need lower triangular, last innovation untransformed
         for (it in 2:ihor) {
            for (ilag in 1:min(nlag,it-1)) {
               responses[,,it]=responses[,,it]+Bv[,,ilag]%*%responses[,,(it-ilag)];
            }
         }

         # Restrict to policy shocks
         imp[(irep-nburn),,,] = responses
      }

      #----- Save draws of the parameters
      alpha_draws[(irep-nburn),] = alpha;
      ALPHA_draws[(irep-nburn),,] = ALPHA;
      SIGMA_draws[(irep-nburn),,] = SIGMA;
      if (prior == 5 || prior == 6) {
         gamma_draws[(irep-nburn),] = gammas; #ok<AGROW>
         if (prior == 6) {
            omega_draws[(irep-nburn),] = omega_vec; #ok<AGROW>
         }
      }
   } # end saving results
} #end the main Gibbs for loop
#====================== End Sampling Posteriors ===========================

#Posterior mean of parameters:
ALPHA_mean = apply(ALPHA_draws,2:3,mean); #posterior mean of ALPHA
ALPHA_mean
SIGMA_mean = apply(SIGMA_draws,2:3,mean); #posterior mean of SIGMA
SIGMA_mean

#Posterior standard deviations of parameters:
ALPHA_std = apply(ALPHA_draws,2:3,sd); #posterior std of ALPHA
SIGMA_std = apply(SIGMA_draws,2:3,std); #posterior std of SIGMA
#or you can use 'ALPHA_COV = cov(alpha_draws,1);' to get the full
#covariance matrix of the posterior of alpha (of dimensions [KM x KM] )

if (prior == 5 || prior == 6) {
   # Find average of restriction indices Gamma
   gammas = apply(gamma_draws,2,mean);
   gammas_mat = matrix(gammas,K,M);
   if (prior == 6) {
      # Find average of restriction indices Omega
      omega = apply(omega_draws,2,mean);
      omega_mat = zeros(M,M);
      for (nn_5 in 1:(M-1)) {
         ggg = omega[((nn_5-1)*(nn_5)/2 + 1):(nn_5*(nn_5+1)/2)];
         omega_mat[1:length(ggg),(nn_5+1)] = ggg;
      }
   }
}

# mean prediction and log predictive likelihood
if (forecasting == 1) {
   Y_pred_mean = apply(Y_pred,2,mean); # mean prediction
   Y_pred_std = apply(Y_pred,2,std);# std prediction
   log_PL = apply(log(PL),2,mean);
   # This are the true values of Y at T+h:
   true_value = Y1[(T+1),];
   #(subsequently you can easily also get MSFE and MAFE)
  
   #======PLOT posterior predictive
   
   par(mfcol = c(M,1), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
   for (i in 1:M) {
      plot(density(Y_pred[,i]),main=colnames(Yraw)[i],xlab="",ylab="",las=1,xaxs="i",tck=.02);
   }
}

# You can also get other quantities, like impulse responses
if (impulses==1) {
   # Set quantiles from the posterior density of the impulse responses
   qus = c(0.1,0.5,0.9)
   imp_resp = apply(imp,2:4,quantile,qus)
   par(mfcol = c(M,M), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
   for (j in 1:M) {
      for (i in 1:M) {
         plot(imp_resp[2,i,j,],type="l",xlab="",ylab="",main=paste0(colnames(Yraw)[i],"-",colnames(Yraw)[j]),las=1,xaxs="i",tck=.02,ylim=c(min(imp_resp),max(imp_resp)))
         lines(imp_resp[1,i,j,],col=2)
         lines(imp_resp[3,i,j,],col=2)
         abline(h=0,lty=2)
      }
   }
}
 
# Clear screen and print elapsed time
toc()

### END
