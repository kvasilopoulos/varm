
### BVAR_GIBBS
library("MASS")
library("matlab")
library("mvtnorm")
library("R.matlab")

bvardgp = function(T=NULL,N=NULL,L=NULL,PHI=NULL,PSI=NULL,const=TRUE) {
   #--------------------------------------------------------------------------
   #   PURPOSE:
   # Get matrix of Y generated from a VAR model
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
      T = 400;           #Number of time series observations (T)
   }
   if (is.null(N)) {
      N = 6;             #Number of cross-sectional observations (N)
   }
   if (is.null(L)) {
      L = 1;             #Lag order
   }
   if (is.null(PHI)) {
      PHI = rbind(rep(1,N),0.5*diag(N))
   }
   if (is.null(PSI)) {
      PSI = matrix(c(1.70, 0.02, 0.16, 0.08, 0.19, -0.03,
                     0.02, 0.66, 0.16, 0.03, 0.26, 0.07,
                     0.16, 0.16, 0.94, 0.25, 0.06, 0.24,
                     0.08, 0.03, 0.25, 1.45, 0.06, 0.38,
                     0.19, 0.26, 0.06, 0.06, 0.64, 0.16,
                     -0.03, 0.07, 0.24, 0.38, 0.16, 1.33),N,N);
   }
   
   #---------------------------------------
   # Ask user if a constant is desired
   if (const) {
      const = 1;
      print(paste0('VAR with ',L ,'-lag(s) and an intercept generated'))
   } else {
      const = 0;
      print(paste0('VAR with ',L ,'-lag(s) and NO intercept generated'))
   }
   
   #----------------------GENERATE--------------------------
   # Set storage in memory for y
   # First L rows are created randomly and are used as 
   # starting (initial) values 
   y = rbind(runif(L*N), zeros(T,N));
   
   # Now generate Y from VAR (L,PHI,PSI)
   u = t(chol(PSI))%*%matrix(rnorm(N*T),N,T);
   ylag = embed(y,L+1)[,-c(1:N)];
   y = cbind(const, ylag)%*%PHI + t(u);
   y
}
Yraw = bvardgp()
Yraw = as.matrix(read.table("/Users/user/Dropbox/18. Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/R/Yraw.dat"))

constant = 1;        # 1: if you desire intercepts, 0: otherwise 
p = 1;               # Number of lags on dependent variables

forecasting = 1;     # 1: Compute h-step ahead predictions, 0: no prediction
forecast_method = 1; # 0: Direct forecasts, 1: Iterated forecasts
h = 2;               # Number of forecast periods
repfor = 10;         # Number of times to obtain a draw from the predictive 
                     # density, for each generated draw of the parameters                     
impulses = 1;        # 1: compute impulse responses, 0: no impulse responses
ihor = 30;           # Horizon to compute impulse responses

# Set prior for BVAR model:
prior = 2;  # prior = 1 --> Indepependent Normal-Whishart Prior
            # prior = 2 --> Indepependent Minnesota-Whishart Prior

# Gibbs-related preliminaries
nsave = 5000;        # Final number of draws to save
nburn = 5000;        # Draws to discard (burn-in)
ntot = nsave+nburn;  # Total number of draws
it_print = 100;      # Print on the screen every "it_print"-th iteration

#--------------------------DATA HANDLING-----------------------------------
# Get initial dimensions of dependent variable
Traw = dim(Yraw)[1]
M = dim(Yraw)[2]

# The model specification is different when implementing direct forecasts,
# compared to the specification when computing iterated forecasts.
if (forecasting==1) {
  if (h<=0) {
    print('You have set forecasting, but the forecast horizon choice is wrong')
  }
  # Now create VAR specification according to forecast method
  if (forecast_method==0) {       # Direct forecasts
    Y1 = Yraw[(h+1):nrow(Yraw),];
    Y2 = Yraw[2:(nrow(Yraw)-h),];
    Traw = Traw-h-1;
  } else if (forecast_method==1) {   # Iterated forecasts
    Y1 = Yraw;
    Y2 = Yraw;
  } else {
    print('Wrong choice of forecast_method')
  }
} else {
  Y1 = Yraw;
  Y2 = Yraw;
}

# Generate lagged Y matrix. This will be part of the X matrix
Ylag = embed(Y2,p+1)[,-c(1:M)]; # Y is [T x M]. ylag is [T x (Mp)]

# Now define matrix X which has all the R.H.S. variables (constant, lags of
# the dependent variable and exogenous regressors/dummies)
if (constant) {
  X1 = cbind(1, Ylag)
} else {
  X1 = Ylag;
}

# Get size of final matrix X
Traw3 = dim(X1)[1]
K = dim(X1)[2]

# Create the block diagonal matrix Z
Z1 = kronecker(diag(M),X1);

# Form Y matrix accordingly
# Delete first "LAGS" rows to match the dimensions of X matrix
Y1 = Y1[(p+1):Traw,]; # This is the final Y matrix used for the VAR

# Traw was the dimesnion of the initial data. T is the number of actual 
# time series observations of Y and X (we lose the p-lags)
T = Traw - p;

#========= FORECASTING SET-UP:
# Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
if (forecasting==1) {
  Y_pred = matrix(0,nsave*repfor,M) # Matrix to save prediction draws
  PL = matrix(0,nsave,1)            # Matrix to save Predictive Likelihood
  if (forecast_method==0) { # Direct forecasts, we only need to keep the 
    Y = Y1[1:(nrow(Y1)-1),]                             # last observation
    X = X1[1:(nrow(X1)-1),]
    Z = kronecker(diag(M),X)
    T = T - 1;
  } else if (forecast_method==1) {              # Iterated forecasts, we keep the last h observations
    Y = Y1[1:(nrow(Y1)-h),]
    X = X1[1:(nrow(X1)-h),]
    Z = kronecker(diag(M),X)
    T = T - h
  }
} else {
  Y = Y1
  X = X1
  Z = Z1
}

#========= IMPULSE RESPONSES SET-UP:
# Create matrices to store forecasts
if (impulses == 1) {
  # Impulse response horizon
   imp = zeros(nsave,M,M,ihor);
   bigj = zeros(M,M*p);
   bigj[1:M,1:M] = eye(M);
}

#-----------------------------PRELIMINARIES--------------------------------
# First get ML estimators
A_OLS = solve(crossprod(X))%*%crossprod(X,Y); # This is the matrix of regression coefficients
a_OLS = c(A_OLS);         # This is the vector of parameters, i.e. it holds that a_OLS = vec(A_OLS)
SSE = crossprod(Y-X%*%A_OLS);   # Sum of squared errors
SIGMA_OLS = SSE/(T-K-1)

# Initialize Bayesian posterior parameters using OLS values
alpha = a_OLS;     # This is the single draw from the posterior of alpha
ALPHA = A_OLS;     # This is the single draw from the posterior of ALPHA
SSE_Gibbs = SSE;   # This is the single draw from the posterior of SSE
SIGMA = SIGMA_OLS; # This is the single draw from the posterior of SIGMA

# Storage space for posterior draws
alpha_draws = zeros(nsave,K*M);
ALPHA_draws = zeros(nsave,K,M);
SIGMA_draws = zeros(nsave,M,M);

#-----------------Prior hyperparameters for bvar model
n = K*M; # Total number of parameters (size of vector alpha)
# Define hyperparameters
if (prior == 1) { # Normal-Wishart
  a_prior = 0*ones(n,1); #<---- prior mean of alpha (parameter vector)
  V_prior = 10*eye(n);    #<---- prior variance of alpha

  # Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
  v_prior = M;             #<---- prior Degrees of Freedom (DoF) of SIGMA
  S_prior = diag(M);       #<---- prior scale of SIGMA
  inv_S_prior = solve(S_prior);
} else if (prior == 2) { # Minnesota-Whishart
  # Prior mean on VAR regression coefficients
  A_prior = rbind(matrix(0,1,M), 0.9*diag(M), matrix(0,(p-1)*M,M));  #<---- prior mean of ALPHA (parameter matrix) 
  a_prior = as.matrix(c(A_prior));               #<---- prior mean of alpha (parameter vector)

  # Minnesota Variance on VAR regression coefficients
  # First define the hyperparameters 'a_bar_i'
  a_bar_1 = 0.5;
  a_bar_2 = 0.5;
  a_bar_3 = 10^2;

  # Now get residual variances of univariate p-lag autoregressions. Here
  # we just run the AR(p) model on each equation, ignoring the constant
  # and exogenous variables (if they have been specified for the original VAR model)
  sigma_sq = zeros(M,1); # vector to store residual variances
  for (i in 1:M) {
    # Create lags of dependent variable in i-th equation
    Ylag_i = matrix(embed(Yraw[,i],(p+1))[,-1],ncol=p);
    Ylag_i = Ylag_i[1:(Traw-p),];
    # Dependent variable in i-th equation
    Y_i = Yraw[(p+1):Traw,i];
    # OLS estimates of i-th equation
    alpha_i = solve(t(Ylag_i)%*%Ylag_i)%*%(t(Ylag_i)%*%Y_i);
    sigma_sq[i,1] = (1/(Traw-p+1))*t(Y_i - Ylag_i%*%alpha_i)%*%(Y_i - Ylag_i%*%alpha_i);
  }
                               
  # Now define prior hyperparameters.
  # Create an array of dimensions K x M, which will contain the K diagonal
  # elements of the covariance matrix, in each of the M equations.
  V_i = matrix(0,K,M);
                               
  # index in each equation which are the own lags
  ind = matrix(0,M,p);
  for (i in 1:M) {
    ind[i,] = na.omit(seq(i+constant,K,3)[1:2])
  }
  for (i in 1:M) {  # for each i-th equation
    for (j in 1:K) {   # for each j-th RHS variable
      if (constant==1) { # if there is constant in the model use this code:
        if (j==1) {
          V_i[j,i] = a_bar_3*sigma_sq[i,1]; # variance on constant                
        } else if (length(which(j==ind[i,]))>0) {
          V_i[j,i] = a_bar_1/(p^2); # variance on own lags           
        } else {
          for (kj in 1:M) {
            if (length(which(j==ind[kj,]))>0) {
              ll = kj;                   
            }
          } # variance on other lags   
         V_i[j,i] = (a_bar_2*sigma_sq[i,1])/((p^2)*sigma_sq[ll,1])
        }
      } else { # if there is no constant in the model use this:
        if (length(which(j==ind[i,]))>0) {
          V_i[j,i] = a_bar_1/(p^2); # variance on own lags
        } else {
          for (kj in 1:M) {
            if (length(which(j==ind[kj,]))>0) {
              ll = kj
            }
          } # variance on other lags  
          V_i[j,i] = (a_bar_2*sigma_sq[i,1])/((p^2)*sigma_sq[ll,1]);
        }
      }
    }
  }
                               
  # Now V is a diagonal matrix with diagonal elements the V_i
  V_prior = diag(c(V_i))  # this is the prior variance of the vector alpha
  
  # Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
  v_prior = M;
  S_prior = diag(M);
  inv_S_prior = solve(S_prior);
}

#========================== Start Sampling ================================
#==========================================================================
tic();
for (irep in 1:ntot) {  #Start the Gibbs "loop"
  if (irep%%100 == 0) {
    print(paste0(irep/ntot*100,"%"))
    toc();
  }
  VARIANCE = kronecker(solve(SIGMA),diag(T));
  inv_V_prior = solve(V_prior)
  V_post = solve(inv_V_prior + t(Z)%*%VARIANCE%*%Z);
  a_post = V_post%*%(inv_V_prior%*%a_prior + t(Z)%*%VARIANCE%*%c(Y));
  alpha = a_post + t(chol(V_post))%*%rnorm(n); # Draw of alpha
  ALPHA = matrix(alpha,K,M); # Draw of ALPHA
  
  # Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
  v_post = T + v_prior;
  S_post = S_prior + crossprod(Y-X%*%ALPHA);
  SIGMA = solve(rWishart(1, v_post, solve(S_post))[,,1]); # Draw SIGMA
  
  # Store results  
  if (irep > nburn) {
    #=========FORECASTING:
    if (forecasting==1) {
      if (forecast_method == 0) {   # Direct forecasts
        Y_temp = matrix(0,repfor,M);
        # compute 'repfor' predictions for each draw of ALPHA and SIGMA
        for (ii in 1:repfor) {
          X_fore = matrix(c(1, Y[T,], X[T,2:(M*(p-1)+1)]),nrow=1);
          # Forecast of T+1 conditional on data at time T
          Y_temp[ii,] = X_fore%*%ALPHA + matrix(rnorm(M,0,1),nrow=1)%*%chol(SIGMA);
        }
        # Matrix of predictions
        Y_pred[(((irep-nburn)-1)*repfor+1):((irep-nburn)*repfor),] = Y_temp;
        # Predictive likelihood
        PL[(irep-nburn),] = dmvnorm(Y1[T+1,], X[T,]%*%ALPHA, SIGMA)
        if (PL[(irep-nburn),] == 0) {
          PL[(irep-nburn),] = 1;
        }
      } else if (forecast_method == 1) {   # Iterated forecasts
        # The usual way is to write the VAR(p) model in companion
        # form, i.e. as VAR(1) model in order to estimate the
        # h-step ahead forecasts directly (this is similar to the 
        # code we use below to obtain impulse responses). Here we 
        # just iterate over h periods, obtaining the forecasts at 
        # T+1, T+2, ..., T+h iteratively.
        Y_temp2 = matrix(NA, nrow=repfor, ncol=M)
        for (ii in 1:repfor) {
          # Forecast of T+1 conditional on data at time T
          xxx = ifelse((M*(p-1)+1)>2,X[T,2:(M*(p-1)+1)],NA)
          X_fore = c(1, Y[T,], xxx);
          X_fore = matrix(c(na.omit(X_fore)),nrow=1)
          Y_hat = X_fore%*%ALPHA + rnorm(M,0,1)%*%chol(SIGMA);
          Y_temp = Y_hat;
          X_new_temp = X_fore;
          if ((h-1)!=0) {
             for (i in 1:(h-1)) {  # Predict T+2, T+3 until T+h                   
               if (i <= p) {
                 # Create matrix of dependent variables for                       
                 # predictions. Dependent on the horizon, use the previous                       
                 # forecasts to create the new right-hand side variables
                 # which is used to evaluate the next forecast.                       
                 X_new_temp = c(1, Y_hat)#, X_fore[2:(M*(p-i)+1)]);
                 
                 # This gives the forecast T+i for i=1,..,p      
                 Y_temp = X_new_temp%*%ALPHA + rnorm(M)%*%chol(SIGMA);                       
                 Y_hat = c(Y_hat, Y_temp);
               } else {
                 X_new_temp = c(1, Y_hat[,1:(M*p)]);
                 Y_temp = X_new_temp%*%ALPHA + rnorm(M)%*%chol(SIGMA);              
                 Y_hat = c(Y_hat, Y_temp);
               }
             } #  the last value of 'Y_temp' is the prediction T+h
             Y_temp2[ii,] = Y_temp;
          }
        }

        # Matrix of predictions               
        Y_pred[(((irep-nburn)-1)*repfor+1):((irep-nburn)*repfor),] = Y_temp2;
        # Predictive likelihood
        PL[(irep-nburn),] = dmvnorm(Y1[(T+h),],X_new_temp%*%ALPHA,SIGMA);
        if (PL[(irep-nburn),] == 0) {
           PL[(irep-nburn),] = 1;
        }
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
  } # end saving results
} #end the main Gibbs for loop
#====================== End Sampling Posteriors ===========================

#Posterior mean of parameters:
ALPHA_mean = apply(ALPHA_draws,2:3,mean); #posterior mean of ALPHA
ALPHA_mean
SIGMA_mean = apply(SIGMA_draws,2:3,mean); #posterior mean of SIGMA
SIGMA_mean

# mean prediction and log predictive likelihood
if (forecasting == 1) {
  Y_pred_mean = apply(Y_pred,2,mean);
  log_PL = apply(log(PL),2,mean);
  #This are the true value of Y at T+h:
  if (forecast_method == 0) {
    true_value = Y1[(T+1),];
  } else if (forecast_method == 1) {
    true_value = Y1[(T+h),];
  } #(subsequently you can easily also get MSFE and MAFE)
}

# You can also get other quantities, like impulse responses
if (impulses==1) {
   # Set quantiles from the posterior density of the impulse responses
   qus = c(0.1,0.5,0.9)
   imp_resp = apply(imp,2:4,quantile,qus)
   par(mfcol = c(M,M), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
´   for (j in 1:M) {
      for (i in 1:M) {
         plot(imp_resp[2,i,j,],type="l",xlab="",ylab="",main=paste0(colnames(Yraw)[i],"-",colnames(Yraw)[j]),las=1,xaxs="i",tck=.02,ylim=c(min(imp_resp),max(imp_resp)))
         lines(imp_resp[1,i,j,],col=2)
         lines(imp_resp[3,i,j,],col=2)
         abline(h=0,lty=2)
      }
   }
}

### END
