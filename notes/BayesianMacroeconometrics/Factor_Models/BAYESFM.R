#BFM - Bayesian Factor Model

# Generate data
X = bfmdgp();
X = as.matrix(read.table("/Users/user/Dropbox/18. Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/R/Yraw.dat"))

T=dim(X)[1]
N=dim(X)[2]

#Demean xraw
for (i in 1:N) {
   X[,i] = X[,i]-mean(X[,i])
}

# Number of factors
K=2;
#----------------------------PRELIMINARIES---------------------------------
# Set some Gibbs - related preliminaries
nods = 500; # Number of draws to keep
bid = 200;   # Number of burn-in-draws
ndraws = nods + bid; # Total number of draws

# store draws in:
Ldraw=array(0,c((nods-bid)/thin,N,K));
Fdraw=array(0,c((nods-bid)/thin,T,K));  
Rdraw=array(0,c((nods-bid)/thin,N,1));
#********************************************************
# STANDARDIZE for PC only
X_st = X
for (i in 1:N) {
   X_st[,i] = X[,i]/sd(X[,i])
}

# First step - extract PC from X
EX=extract(X_st,K);
F0=EX$fac
Lf=EX$lam

# Transform factors and loadings for LE normalization
Q=qr(t(Lf));
ql=qr.Q(Q)
rl=qr.R(Q)
Lf=rl;  # do not transpose yet, is upper triangular
F0=F0%*%ql;

# Need identity in the first K columns of Lf, call them A for now
A=Lf[,1:K];
Lf=rbind(diag(K),t(solve(A)%*%Lf[,(K+1):N]));
F0=F0%*%A;

# Obtain R:
e=X_st-F0%*%t(Lf);
R=t(e)%*%e/T;
R=diag(diag(R));
L=Lf;
F=matrix(rnorm(T*K,0,1),T,K);
#------------------PRIORS:--------------
L_var_prior=diag(K);
           
# Hyperparameters for R(n,n) ~ IG(a/2,b/2)
mu_sigeps = 0.5;
a = 3;
b = (1/(2*mu_sigeps));
v0 = 5;
s0 = 0.1;

#==========================================================================
#========================== Start Sampling ================================
#==========================================================================

#************************** Start the Gibbs "loop" ************************
print('Number of iterations');
for (rep in 1:nods) {
  if (rep%%100 == 0) {
    print(paste0(rep/nods*100,"%"));
  }

  # STEP 1. =========|DRAW FACTORS
  for (t in 1:T) {
    F_var = solve(diag(K) + t(L)%*%solve(R)%*%L);
    F_mean = F_var%*%t(L)%*%solve(R)%*%matrix(X[t,],ncol=1);
    Fd = F_mean + t(chol(F_var))%*%rnorm(K);
    F[t,1:K]=t(Fd);
  }
  
  # STEP 2. =========|DRAW LOADINGS + SIGMA
  for (n in 1:N) {
    ed=X[,n]-F%*%L[n,];
    #R[n,n] = invgamrnd((T/2) + a, solve(solve(b) + .5*crossprod(ed)));
    
    # draw R(n,n)
    R_2 = (s0 + crossprod(ed)/2);
    R_1 = (v0 + T/2);
    rd = rgamma(1,R_1,R_2);
    Rd = solve(rd);
    R[n,n]=Rd;
                                            
    # draw L(n,1:K):
    if (n <= K) {
      Lg=zeros(K,1);
      L_var_post=solve(solve(diag(n))+1/(R[n,n])*crossprod(F[,1:n]));
      Lpostmean = L_var_post%*%(1/(R[n,n])*t(F[,1:n])%*%X[,n]);
      if (n<1) {
        Lg[1:(n-1),] = Lpostmean[1:(n-1),] + t(chol(L_var_post[1:(n-1),1:(n-1)]))%*%rnorm(n-1);
      }
      Lg[n,] = rnorm(1,Lpostmean[n],L_var_post[n,n]);
      L[n,1:K]=t(Lg);
    } else if (n>K) {
      L_var_post=solve(solve(L_var_prior)+1/R[n,n]*crossprod(F));
      Lpostmean=(1/R[n,n])*L_var_post%*%t(F)%*%X[,n];
      Ld=Lpostmean + t(chol(L_var_post))%*%rnorm(K);
      L[n,1:K]=t(Ld);
    }
  }

  if (rep>bid) {
    Ldraw[(rep-bid),,]=L;
    Fdraw[(rep-bid),,]=F;
    Rdraw[(rep-bid),,]=diag(R);
  }
}

# Average over Gibbs draws
Ldraw2=apply(Ldraw,2:3,mean);
Ldraw2
Fdraw2=apply(Fdraw,2:3,mean);
Fdraw2
Rdraw2=apply(Rdraw,2:3,mean);
Rdraw2

par(mfcol = c(K,1), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:K) {
   plot(Fdraw2[,i],type="l",xaxs="i",las=1,xlab="",ylab="",tck=.02)
}

### END
