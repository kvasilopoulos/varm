
### SSVS_VAR_CONST
library("matlab")
library("R.matlab")

### FUNCTIONS:
simvardgp = function (T=NULL,N=NULL,L=NULL,PHI=NULL,PSI=NULL,const=TRUE) {
   #--------------------------------------------------------------------------
   #PURPOSE:
   #Get matrix of Y generated from a VAR model
   #--------------------------------------------------------------------------
   #INPUTS:
   #  T  - Number of observations (rows of Y)
   #  N  - Number of series (columns of Y)
   #  L  - Number of lags
   #
   #OUTPUT:
   #  y  - [T x N] matrix generated from VAR(L) model
   # -------------------------------------------------------------------------
   
   #-----------------------PRELIMINARIES--------------------
   if (is.null(T)) {
      T = 1000; #Number of time series observations (T)
   }
   if (is.null(N)) {
      N = 6; #Number of cross-sectional observations (N)
   }
   if (is.null(L)) {
      L = 1; #Lag order
   }
   if (is.null(PSI)) {
      PSI = eye(N)
      PSI[1,-1] = 0.5
   }
   if (is.null(PHI)) {
      PHI = rbind(rep(1,6),0.9*eye(N))
   }
   
   
   #----------------------GENERATE--------------------------
   # Set storage in memory for y
   # First L rows are created randomly and are used as 
   # starting (initial) values 
   y = rbind(matrix(runif(L*N),L,N), zeros(T,N));
   
   # Now generate Y from VAR (L,PHI,PSI)
   sigma = solve(crossprod(t(PSI)))
   u = chol(sigma)%*%matrix(rnorm(N*T),N);
   ylag = embed(y,L+1)[,-c(1:ncol(y))];
   y = cbind(const, ylag)%*%PHI + t(u);
   y
}

# SSVS_VAR_CONST.m
#
#This code implements the Stochastic Search Variable Selection
#algorithm of George and McCulloch (1993), JASA, for VAR models. The
#difference with SSVS_VAR.m is that the constant in this case is left
#unrestricted (i.e. no mixture prior, just the Normal prior with
#large variance).
#
#--------------------------------------------------------------------------
# Based on George, Sun and Ni (2008) "Bayesian Stochastic Search for VAR 
# Model Restrictions", Journal of Econometrics.
#
# Note that this code is based on the WORKING paper version of the paper
# and hence reference to equations in the comments may be different from the journal article.
#
# In the comments you will see references to specific equations in the 
# working paper mentioned above. You can also check the .pdf file I provide
# (Program SSVS) for implementation details. I use the notation of the
# monograph (Koop and Korobilis, 2009).
#--------------------------------------------------------------------------

#-------------------------LOAD DATA----------------------------------------
# Simulate a VAR model with the same parameters as in the first simulation 
# of the George, Sun and Ni (2008) paper
#y = simvardgp();
Yraw = as.matrix(read.table("/Users/user/Dropbox/Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/R/Yraw.dat"))
y = Yraw
p=1; #choose number of lags

#--------------------------DATA HANDLING-----------------------------------
Traw = nrow(y)
M = ncol(y);
h = 0; # Define if there are exogenous variables (h/=0) or not (h=0).Leave value to zero for this demonstrating code.

# Generate lagged Y matrix. This will be part of the X matrix
ylag = embed(y,p+1)[,-c(1:M)]; # Y is [T x M]. ylag is [T x (pM)]

# Insert "ylag" defined above to form X matrix
if (h) {# If we have exogenous variables (Z), then X is [T x (1 + h + pM)]
   X = cbind(ones(Traw-p,1), ylag, z[(p+1):Traw,]);
} else {# If no Z is defined, then X is [T x (1 + pM)]
   X = cbind(ones(Traw-p,1), ylag);
}

# Form y matrix accordingly
# Delete first "p" rows to match the dimensions of X matrix
Y = y[(p+1):Traw,]; # This is the final Y matrix used for the VAR

#Traw was the dimesnion of initial data. T is the number of actual 
#time series observations of Y and X (final Y & X)
T=Traw-p;

n = (1 + h + p*M)*M;  # n is the number of alpha parameter elements (is the
                      # number of rows of the "stacked" column vector of parameters)  

# In this program we assume that the constants are unrestricted
# and that the lag coefficients are restricted
m = (n - M); # choose number of restrictions, n minus the # of intercepts
non = n - m; # unrestricted alphas (= M intercepts)
#----------------------------PRELIMINARIES---------------------------------
# Set some Gibbs - related preliminaries
nsave = 100;  # Number of draws to save
nburn = 100;  # Number of burn-in-draws
ntot = nsave + nburn; # Total number of draws
thin = 3;       # Consider every thin-th draw (thin value)
it_print = 500; # Print every "it_print"-th iteration

# Set storage space in computer memory for parameter matrices
alpha_draws = zeros(n,nsave); # store alphas draws
psi_ii_sq = zeros(M,1); # vector of psi^2 drawn at each iteration
psi_ii_sq_draws = zeros(M,nsave); # store psi^2 draws
gammas_draws = zeros(m,nsave); # store gamma draws
omega_draws = zeros(0.5*M*(M-1),nsave); # store omega draws
psi_mat_draws = zeros(M,M,nsave); # store draws of PSI matrix
alpha_mat_draws = zeros(n/M,M,nsave); # store draws of ALPHA matrix

# Set storage space in memory for cell arrays. These can be put
# individually before each "for" loopS=cell(1,M);
S=cell(1,M);
s=cell(1,M-1);
omega=cell(1,M-1);
h=cell(1,M-1);
D_j=cell(1,M-1);
R_j=cell(1,M-1);
DRD_j=cell(1,M-1);
B=cell(1,M);
eta=cell(1,M-1);

#---------------Prior hyperparameters for ssvs algorithm-------------------
# First get ML estimators, see eq. (8)
ALPHA_OLS = solve(crossprod(X))%*%crossprod(X,Y);
SSE = t(Y - X%*%ALPHA_OLS)%*%(Y - X%*%ALPHA_OLS);
SSE = crossprod(Y - X%*%ALPHA_OLS)
# Stack columns of ALPHA_M
alpha_OLS_vec=matrix(ALPHA_OLS,n,1);  # vector of "MLE" alphas (=vec(ALPHA_OLS))

# Variances for the "Phi mixture", see eq.(13)
tau_0 = 0.1; #*ones(m,1);  % Set tau_[0i], tau_[1i]
tau_1 = 9;  #*ones(m,1);

# Variances for the "Eta mixture", see eq.(16)  
kappa_0 = 0.1;#*ones(m,M-1); % Set kappa_[0ij], kappa_[1ij]
kappa_1 = 6; #*ones(m,M-1);

# Hyperparameters for Phi_non ~ N_[n-m](f_non,M_non), see eq.(11)
f_non = 1*ones(non,1);
c = 0.1;
M_non = c*eye(non);

# Hyperparameters for Phi_m ~ N_[m](f_m,DRD), see eq.(11)
f_m = 0*ones(m,1); #Mean of restricted alphas 

f_tot = zeros((n/M),M); #Prior mean of all alphas(restr. & unrestr.)
for (ged in 1:M) {
   nnn_9 = 1 + (m/M)*(ged-1);
   f_tot[,ged] = c(f_non[ged,], f_m[nnn_9:((m/M)*ged)]);
}
f_tot = matrix(f_tot,nrow(f_tot)*ncol(f_tot),1);

# Matrices of restrictions on Phi and Psi parameters, R & R_[j]
R = eye(m);# Set here (square) R matrix of m restrictions
# Create R_[j] = R_1, R_2, R_3,... matrices of restrictions
for (kk in 1:(M-1)) {	# Set matrix R_j of restrictions on psi
   R_j[[kk]] = eye(kk);
}

# Initialize Gamma and Omega vectors
gammas = ones(m,1); # vector of Gamma

for (kk_1 in 1:(M-1)) {
   omega[[kk_1]] = ones(kk_1,1);	# Omega_j
}

# Hyperparameters for Gamma ~ BERNOULLI(m,p_i), see eq. (14)
p_i = .5;

# Hyperparameters for Omega_[j] ~ BERNOULLI(j,q_ij), see eq. (17)
q_ij = .7;

# Hyperparameters for (Psi)^2 ~ GAMMA(a_i , b_i), see eq. (18)
a_i = .01;
b_i = .01;
#***************End of Preliminaries & PriorSpecification******************


#========================== Start Sampling ================================
tic(); # Start the timer
#************************** Start the Gibbs "loop" ************************
print('Number of iterations');
for (irep in 1:ntot) {
   if (irep%%it_print==0) {
      print(paste(irep/ntot*100,"%"));
      toc();
   }
   # STEP 1.: ----------------------------Draw "psi"----------------------
   # Draw psi|alpha,gamma,omega,DATA from the GAMMA dist.
   #----------------------------------------------------------------------

   #   Get S_[j] - upper-left [j x j] submatrices of SSE
   # The following loop creates a cell array with elements S_1,
   # S_2,...,S_j with respective dimensions 1x1, 2x2,...,jxj
   for (kk_2 in 1:M) {
      S[[kk_2]] = SSE[1:kk_2,1:kk_2];
   }

   # Set also SSE =(s_[i,j]) & get vectors s_[j]=(s_[1,j] , ... , s_[j-1,j])
   for (kk_3 in 2:M) {
      s[[kk_3-1]] = SSE[1:(kk_3-1),kk_3];
   }

   # Parameters for Heta|omega ~ N_[j-1](0,D_[j]*R_[j]*D_[j]), see eq. (15)
   # Create and update h_[j] matrix
   # If omega_[ij] = 0 => h_[ij] = kappa0, else...
   for (kk_4 in 1:(M-1)) {
      omeg = omega[[kk_4]];
      het = h[[kk_4]];
      for (kkk in 1:nrow(omeg)) {
         if (omeg[kkk,] == 0) {
            het[kkk] = kappa_0;
         } else {
            het[kkk] = kappa_1;
         }
      }
      h[[kk_4]] = het;
   }

   # D_j = diag(h_[1j],...,h[j-1,j])
   for (kk_5 in 1:(M-1)) {
      D_j[[kk_5]] = diag(length(h[[kk_5]]))*h[[kk_5]]
   }

   # Now create covariance matrix D_[j]*R_[j]*D_[j], see eq. (15)
   for (kk_6 in 1:(M-1)) {
      DD = D_j[[kk_6]];
      RR = R_j[[kk_6]];
      DRD_j[[kk_6]] = DD%*%RR%*%DD;
   }

   # Create B_[i] matrix
   for (rr in 1:M) {
      if (rr == 1) {
         B[[rr]] = b_i + 0.5*(SSE[rr,rr]);
      } else if (rr > 1) {
         s_i = s[[rr-1]];
         S_i = S[[rr-1]];
         DiRiDi = DRD_j[[rr-1]];
         B[[rr]] = b_i + 0.5*(SSE[rr,rr] - t(s_i)%*%solve(S_i + solve(DiRiDi))%*%s_i);
      }
   }

   # Now get B_i from cell array B, and generate (psi_[ii])^2
   B_i = unlist(B);
   for (kk_7 in 1:M) {
      # If you have the Statistics toolbox, use "gamrnd" instead
      psi_ii_sq[kk_7,1] = rgamma(1,(a_i + 0.5*T),B_i[kk_7]);
   }

   # STEP 2.: ----------------------------Draw "eta"----------------------
   # Draw eta|psi,alpha,gamma,omega,DATA from the [j-1]-variate NORMAL dist.
   #----------------------------------------------------------------------
   for (kk_8 in 1:(M-1)) {
      s_i = s[[kk_8]];
      S_i = S[[kk_8]];
      DiRiDi = DRD_j[[kk_8]];
      Delta_j = solve(S_i + solve(DiRiDi));
      miu_j = - sqrt(psi_ii_sq[kk_8+1])*Delta_j%*%s_i;
      eta[[kk_8]] = miu_j + t(chol(Delta_j))%*%rnorm(kk_8);
   }

   # STEP 3.: --------------------------Draw "omega"----------------------
   # Draw omega|eta,psi,alpha,gamma,omega,DATA from BERNOULLI dist.
   #---------------------------------------------------------------------- 
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
      omega[[kk_9]] = omeg_g;
   }

   # STEP 4.: --------------------------Draw "alpha"------------------------
   # Draw alpha|gamma,Sigma,omega,DATA from NORMAL dist.
   # ---------------------------------------------------------------------- 

   # Create PSI matrix from individual elements of "psi_ii_sq" and "eta"
   PSI_ALL = zeros(M,M);
   for (nn_1 in 1:M) {
      PSI_ALL[nn_1,nn_1] = sqrt(psi_ii_sq[nn_1,1]);
   }

   for (nn_2 in 1:(M-1)) {
      eta_gg = eta[[nn_2]];
      for (nnn in 1:nrow(eta_gg)) {
         PSI_ALL[nnn,(nn_2+1)] = eta_gg[nnn];
      }
   }

   # Hyperparameters for Phi_m|gamma ~ N_[m](0,D*R*D), see eq.(12)
   h_i = zeros(m,1);# h_i is tau_0 if gamma=0 and tau_1 if gamma=1
   for (nn_3 in 1:m) {
      if (gammas[nn_3,1] == 0) {
         h_i[nn_3,1] = tau_0;
      } else if (gammas[nn_3,1] == 1) {
         h_i[nn_3,1] = tau_1;
      }
   }
   D = diag(c(h_i)); # Create D. Here D=diag(h_i) will also do
   DRD = D%*%R%*%D;# Prior covariance matrix for Phi_m
   psi_xx = kronecker(crossprod(t(PSI_ALL)),crossprod(X));
   temp1 = zeros(nrow(M_non),nrow(DRD));
   temp2 = zeros((n/M),M);
   M_vec = diag(M_non);
   DRD_vec = diag(DRD);
   for (sed in 1:M) {
      nnn_9 = 1 + (m/M)*(sed-1);
      temp2[,sed] = c(M_vec[sed], DRD_vec[nnn_9:((m/M)*sed)]);
   }
   temp2 = matrix(temp2,nrow(temp2)*ncol(temp2),1);
   temp2 = diag(c(temp2));
   Delta_alpha = solve(psi_xx + solve(temp2));
 
   miu_alpha = Delta_alpha%*%(psi_xx%*%alpha_OLS_vec + solve(temp2)%*%f_tot);
   alphas = miu_alpha + t(chol(Delta_alpha))%*%rnorm(n,1);
 
   alpha_mat = matrix(alphas,n/M,M);
   alpha_temp = alpha_mat[2:(n/M),];
   
   # STEP 5.: --------------------------Draw "gamma"----------------------
   # Draw gamma|alpha,psi,eta,omega,DATA from BERNOULLI dist.
   #---------------------------------------------------------------------- 
   for (nn_6 in 1:m) {
      u_i1 = (1/tau_0)*exp(-0.5*(c(alpha_temp)[nn_6]/(tau_0))^2)*p_i;
      u_i2 = (1/tau_1)*exp(-0.5*(c(alpha_temp)[nn_6]/(tau_1))^2)*(1-p_i);
      gst = u_i1/(u_i1 + u_i2);
      gammas[nn_6,1] = rbinom(1,1,1-gst);
   }
 
   # Save new Sum of Squared Errors (SSE)
   SSE = crossprod(Y - X%*%alpha_mat);
 
   # Store matrices
   if (irep>nburn) {
      alpha_draws[,irep-nburn] = alphas;
      psi_ii_sq_draws[,irep-nburn] = psi_ii_sq;
      psi_mat_draws[,,irep-nburn] = PSI_ALL;
      alpha_mat_draws[,,irep-nburn] = alpha_mat;
      gammas_draws[,irep-nburn] = gammas;
      omega_draws[,irep-nburn] = omega_vec;
   }
}

# ========================== Finished Samplin g============================
# =========================================================================
print('Finished sampling successfully')
 
# Do thining in case of high correlation
thin_val = seq(1,floor(nsave/thin),thin)
 
alpha_draws = alpha_draws[,thin_val];
psi_ii_sq_draws = psi_ii_sq_draws[,thin_val];
psi_mat_draws = psi_mat_draws[,,thin_val];
alpha_mat_draws = alpha_mat_draws[,,thin_val];
gammas_draws = gammas_draws[,thin_val];
omega_draws = omega_draws[,thin_val];

# Find average of restriction indices Gamma (note that the constant is not included!)
gammas = t(t(apply(gammas_draws,1,mean)));
gammas_mat = zeros(p*M,M);
for (nn_5 in 1:M) {
   nnn_1 = 1 + p*M*(nn_5-1);
   gammas_mat[,nn_5] = gammas[nnn_1:(p*M*nn_5)];
}
 
# Find average of restriction indices Omega
omega = t(t(apply(omega_draws,1,mean)));
omega_mat = zeros(M,M);
for (nn_5 in 1:(M-1)) {
   ggg = omega[((nn_5-1)*(nn_5)/2 + 1):(nn_5*(nn_5+1)/2),];
   omega_mat[1:length (ggg),(nn_5+1)] = ggg;
}
toc()

### END
