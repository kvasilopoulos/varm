
### BVAR_Analytical
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
   ylag = embed(y, L+1)[,-c(1:N)];
   y = cbind(const, ylag)%*%PHI + t(u);
   y
}
Yraw = bvardgp()

Yraw = as.matrix(read.table("/Users/user/Dropbox/Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/R/Yraw.dat"))
# In any case, name the data you load 'Yraw', in order to avoid changing the
# rest of the code. Note that 'Yraw' is a matrix with T rows by M columns,
# where T is the number of time series observations (usually months or
# quarters), while M is the number of VAR dependent macro variables.
#----------------------------PRELIMINARIES---------------------------------
# Define specification of the VAR model
constant = 1;        # 1: if you desire intercepts, 0: otherwise
p = 1;               # Number of lags on dependent variables
forecasting = 1;     # 1: Compute h-step ahead predictions, 0: no prediction
forecast_method = 0; # 0: Direct forecasts
                     # 1: Iterated forecasts
h = 4;               # Number of forecast periods
impulses = 1;        # 1: compute impulse responses, 0: no impulse responses
ihor = 30;           # Horizon to compute impulse responses

# Set prior for BVAR model:
prior = 2;  # prior = 1 --> Noninformative Prior
            # prior = 2 --> Minnesota Prior
            # prior = 3 --> Natural conjugate Prior

#--------------------------DATA HANDLING-----------------------------------
# Get initial dimensions of dependent variable
Traw=dim(Yraw)[1]
M=dim(Yraw)[2]

if (forecasting==1) {
  if (h<=0) {
    print('You have set forecasting, but the forecast horizon choice is wrong')
  }

  # Now create VAR specification according to forecast method
  if (forecast_method==0) {   # Direct forecasts
    Y1 = Yraw[(h+1):nrow(Yraw),]
    Y2 = Yraw[2:(nrow(Yraw)-h),]
    Traw = Traw-h-1
  } else if (forecast_method==1) {   # Iterated forecasts
    Y1 = Yraw
    Y2 = Yraw
  } else {
    print('Wrong choice of forecast_method')
  }
} else {
  Y1 = Yraw
  Y2 = Yraw
}

Ylag = embed(Y2,(p+1))[,-c(1:ncol(Y2))]
if (constant) {
  X1 = cbind(1,Ylag)
} else {
  X1 = Ylag
}

Traw3 = dim(X1)[1]
K = dim(X1)[2]
Z1 = kronecker(diag(M),X1)
Y1 = Y1[(p+1):Traw,]
T = Traw-p

#========= FORECASTING SET-UP:
# Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
if (forecasting==1) {
  if (forecast_method==0) { # Direct forecasts, we only need to keep the
    Y = Y1[1:(nrow(Y1)-1),] # last observation
    X = X1[1:(nrow(X1)-1),]
    Z = kronecker(diag(M),X)
    T = T - 1
  } else {  # Iterated forecasts, we keep the last h observations
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


### PRIORS
A_OLS = solve(crossprod(X))%*%crossprod(X,Y)
a_OLS = as.matrix(c(A_OLS))
SSE = crossprod(Y-X%*%A_OLS)
SIGMA_OLS = SSE/(T-K)

if (prior==1){
  # Posteriors depend on OLS quantities
} else if (prior==2) {
  A_prior = zeros(K,M)
  a_prior = c(A_prior)

  # Hyperparameters on the Minnesota variance of alpha
  a_bar_1 = 0.5
  a_bar_2 = 0.5
  a_bar_3 = 10^2

  # Now get residual variances of univariate p-lag autoregressions. Here
  # we just run the AR(p) model on each equation, ignoring the constant
  # and exogenous variables (if they have been specified for the original
  # VAR model)
  sigma_sq = matrix(0,M,1) # vector to store residual variances
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
  V_i = matrix(0,K,M)

  # index in each equation which are the own lags
  ind = matrix(0,M,p)

  for (i in 1:M) {
     ind[i,] = na.omit(seq(i+constant,K,3)[1:2])
  }

  for (i in 1:M) {  # for each i-th equation
    for (j in 1:K) {  # for each j-th RHS variable
      if (constant==1) {
        if (j==1) {
          V_i[j,i] = a_bar_3*sigma_sq[i,1] # variance on constant
        } else if (length(which(j==ind[i,]))>0) {
          V_i[j,i] = a_bar_1/(ceiling((j-1)/M)^2) # variance on own lags
        } else {
          for (kj in 1:M) {
            if (length(which(j==ind[kj,]))>0) {
              ll = kj
            }
          }
          V_i[j,i] = (a_bar_2*sigma_sq[i,1])/((ceiling((j-1)/M)^2)*sigma_sq[ll,1])
        }
      } else {
        if (length(which(j==ind[i,]))>0) {
          V_i[j,i] = a_bar_1/(ceiling((j-1)/M)^2) # variance on own lags
        } else {
          for (kj in 1:M) {
            if (length(which(j==ind[kj,]))>0) {
              ll = kj
            }
          }
          V_i[j,i] = (a_bar_2*sigma_sq(i,1))/((ceiling((j-1)/M)^2)*sigma_sq[ll,1])
        }
      }
    }
  }

  # Now V is a diagonal matrix with diagonal elements the V_i
  V_prior = diag(c(V_i))  # this is the prior variance of the vector a

  # SIGMA is equal to the OLS quantity
  SIGMA = SIGMA_OLS

} else if (prior == 3) { # Normal-Wishart (nat conj)
  # Hyperparameters on a ~ N(a_prior, SIGMA x V_prior)
  A_prior = zeros(K,M)
  a_prior = c(A_prior)
  V_prior = 10*diag(K)

  # Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
  v_prior = M
  S_prior = diag(M)
  inv_S_prior = solve(S_prior)
}

#### POSTERIORS
# Posterior hyperparameters of ALPHA and SIGMA with Diffuse Prior
if (prior == 1) {
  # Posterior of alpha|Data ~ Multi-T(kron(SSE,inv(X'X)),alpha_OLS,T-K)
  V_post = solve(t(X)%*%X)
  a_post = a_OLS
  A_post = matrix(a_post,K,M)

  # posterior of SIGMA|Data ~ inv-Wishart(SSE,T-K)
  S_post = SSE
  v_post = T-K

  # Now get the mean and variance of the Multi-t marginal posterior of alpha
  alpha_mean = a_post
  alpha_var = (1/(v_post - M - 1))*kronecker(S_post,V_post)

  #--------- Posterior hyperparameters of ALPHA and SIGMA with Minnesota Prior
} else if (prior == 2) {
  # ******Get all the required quantities for the posteriors
  V_post = solve(solve(V_prior) + kronecker(solve(SIGMA),t(X)%*%X) )
  a_post = V_post %*% (solve(V_prior)%*%a_prior + kronecker(solve(SIGMA),t(X)%*%X)%*%a_OLS )
  A_post = matrix(a_post,K,M)

  # In this case, the mean is a_post and the variance is V_post
  alpha_mean = a_post
  #--------- Posterior hyperparameters of ALPHA and SIGMA with Normal-Wishart Prior
} else if (prior == 3) {
  # ******Get all the required quantities for the posteriors
  # For alpha
  V_post = solve(solve(V_prior) + t(X)%*%X );
  A_post = V_post %*% (solve(V_prior)%*%A_prior + t(X)%*%X%*%A_OLS )
  a_post = c(A_post)

  # For SIGMA
  S_post = SSE + S_prior + t(A_OLS)%*%t(X)%*%X%*%A_OLS + t(A_prior)%*%solve(V_prior)%*%A_prior - t(A_post)%*%( solve(V_prior) + t(X)%*%X )%*%A_post
  v_post = T + v_prior

  # Now get the mean and variance of the Multi-t marginal posterior of alpha
  alpha_mean = a_post
  alpha_var = (1/(v_post-M-1))*kronecker(S_post,V_post)
}

#======================= PREDICTIVE INFERENCE =============================
#==========================================================================
xxx = ifelse(2<(M*(p-1)+1),X[T,2:(M*(p-1)+1)],NA)
X_tplus1 = c(1, Y[T,], xxx)
X_tplus1 = na.omit(X_tplus1)

# As a point forecast use predictive mean
Pred_mean = X_tplus1%*%A_post

nsave = 20000
#========= FORECASTING:
PRED = matrix(NA,ncol=M,nrow=nsave)
#========= IMPULSE RESPONSES SET-UP:
# Create matrices to store forecasts
if (impulses == 1) {
   imp = zeros(nsave,M,M,ihor);
   bigj = zeros(M,M*p);
   bigj[1:M,1:M] = eye(M);
}
for (irep in 1:nsave) {
   ALPHA = matrix(mvrnorm(1,c(alpha_mean),Sigma=V_post),ncol=M)
   if (forecasting==1) {
      PRED[irep,] = X_tplus1%*%ALPHA
   }

   if (impulses==1) {
      Bv = zeros(M,M,p);
      for (i_1 in 1:p) {
         Bv[,,i_1] = ALPHA[(1+((i_1-1)*M + 1)):(i_1*M+1),];
      }

      # st dev matrix for structural VAR
      shock = t(chol(S_post));
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
      imp[irep,,,] = responses
   }
}

if (forecasting==1) {
   par(mfcol = c(M,1), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
   for (i in 1:M) {
      plot(density(PRED[,i]),xlab="",ylab="",las=1,xaxs="i",tck=.02,col="steelblue4",yaxt="n",main=colnames(Yraw)[i])
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
### END
