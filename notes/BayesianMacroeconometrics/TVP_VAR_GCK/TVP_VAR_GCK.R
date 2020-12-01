
gck=function(yg,gg,hh,capg,f,capf,sigv,kold,t,ex0,vx0,nvalk,kprior,kvals,p,kstate) {
# GCK Function to implement Gerlach, Carter and Kohn (2000), 'Efficient 
# Bayesian Inference for Dynamic Mixture Models', JASA. 
#--------------------------------------------------------------------------
# The dynamic mixture model is of the form
# 
# Y[t] = g[t] + h[t] x X[t] + G[t] x u[t] 
# X[t] = f[t] + F[t] x X[t-1] + GAMMA[t] x K[t] x v[t]
#
# where u[t] and v[t] are the errors (normal or conditionally normal).
# Conmditional on the parameters g[t], h[t], G[t], f[t], F[t] amd GAMMA[t],
# this algorithm provides estimates of the indicators K[t] as in Gerlach,
# Carter and Kohn (GCK), Recursion 1, pages 821-822.
#--------------------------------------------------------------------------
# Inputs: yg = p by t data matrix
# gg = p by t matrix of det terms in measurement equation (usually set to zeros).
# hh = tp by m matrix (note that this is transpose of GKCs definition)
# capg = t*p by p matrix containing gt which is akin to the st dev of measurement equation
# f = m by t matrix of det terms in the state equation (usually set to zero in my work)
# capf = tm by m matrix from state equation (set to identies for random walk evolution of states)
# sigv is the standard deviation (or sigv*sigv') is the state equation error variance when 
# regime change occurs (i.e. Kt=1). sigv will be an m by m matrix 
# kold is the previous draw of k, which is t by 1 (i.e. this code only allows for one k)
# t = nunber of observations
# kstate = dimensionality of state vector
# ex0, vx0 = mean and variance for initial values for state equation (m by 1 and m by m)
# nvalk = number of values k can take on -- usually 2 for 0/1
# kprior = prior probabilities for each value for k = nvalk by 1 (this is okay for 
# Bernoulli case, but in general may make this nvalk by t) 
# kvals are the values k can take on -- usually 0/1

# GCK's Step 1 on page 821
lpy2n=0;
mu = zeros(t*kstate,1);
omega = zeros(t*kstate,kstate);
for (i in (t-1):1) {
   gatplus1 = sigv*kold[(i+1)];
   ftplus1 = capf[(kstate*i+1):(kstate*(i+1)),];
   cgtplus1 = capg[(i*p+1):((i+1)*p),];
   htplus1 = t(hh[(i*p+1):((i+1)*p),]);
   
   rtplus1 = (t(htplus1)%*%gatplus1)%*%t(t(htplus1)%*%gatplus1) + crossprod(t(cgtplus1));
   rtinv = solve(rtplus1);
   btplus1 = crossprod(t(gatplus1))%*%htplus1%*%rtinv;
   atplus1 = (eye(kstate) - btplus1%*%t(htplus1))%*%ftplus1;
   if (kold[(i+1)] == 0) {
      ctplus1 = zeros(kstate,kstate);
   } else {
      cct = gatplus1%*%(eye(kstate) - t(gatplus1)%*%htplus1%*%rtinv%*%t(htplus1)%*%gatplus1)%*%t(gatplus1);
      ctplus1 = t(chol(cct));
   }
   otplus1 = omega[(kstate*i+1):(kstate*(i+1)),];
   dtplus1 = t(ctplus1)%*%otplus1%*%ctplus1 + eye(kstate);
   omega[(kstate*(i-1)+1):(kstate*i),] = t(atplus1)%*%(otplus1 - otplus1%*%ctplus1%*%solve(dtplus1)%*%t(ctplus1)%*%otplus1)%*%atplus1 + t(ftplus1)%*%htplus1%*%rtinv%*%t(htplus1)%*%ftplus1;
   satplus1 = (eye(kstate) - btplus1%*%t(htplus1))%*%f[,(i+1)] - btplus1%*%gg[,(i+1)];
   mutplus1 = mu[(kstate*i+1):(kstate*(i+1)),];
   mu[(kstate*(i-1)+1):(kstate*i),] = t(atplus1)%*%(eye(kstate)-otplus1%*%ctplus1%*%solve(dtplus1)%*%t(ctplus1)) %*% (mutplus1-otplus1%*%(satplus1+btplus1%*%yg[,(i+1)]))+t(ftplus1)%*%htplus1%*%rtinv%*%(yg[,(i+1)]-gg[,(i+1)]-t(htplus1)%*%f[,(i+1)]);
}

# GCKs Step 2 on pages 821-822
kdraw = kold;
ht = t(hh[1:p,]);
ft = capf[1:kstate,];
gat = zeros(kstate,kstate);
# Note: this specification implies no shift in first period -- sensible
rt = t(ht)%*%ft%*%vx0%*%t(ft)%*%ht + t(ht)%*%gat%*%t(gat)%*%ht+ capg[1:p,]%*%t(capg[1:p,]);
rtinv = solve(rt);
jt = (ft%*%vx0%*%t(ft)%*%ht + crossprod(t(gat))%*%ht)%*%rtinv;
mtm1 = (eye(kstate) - jt%*%t(ht))%*%(f[,1] + ft%*%ex0) + jt%*%(yg[,1] - gg[,1]);
vtm1 = ft%*%vx0%*%t(ft)+ crossprod(t(gat)) - jt%*%rt%*%t(jt);
lprob = zeros(nvalk,1);
for (i in 2:t) {
   ht = t(hh[((i-1)*p+1):(i*p),]);
   ft = capf[(kstate*(i-1)+1):(kstate*i),];
   for (j in 1:nvalk) {
      gat = kvals[j,1]*sigv;
      rt = t(ht)%*%ft%*%vtm1%*%t(ft)%*%ht + t(ht)%*%crossprod(t(gat))%*%ht + crossprod(t(capg[((i-1)*p+1):(i*p),]));
      rtinv = solve(rt);
      jt = (ft%*%vtm1%*%t(ft)%*%ht + crossprod(t(gat))%*%ht)%*%rtinv;
      mt = (eye(kstate) - jt%*%t(ht))%*%(f[,i] + ft%*%mtm1) + jt%*%(yg[,i] - gg[,i]);
      vt = ft%*%vtm1%*%t(ft) + crossprod(t(gat)) - jt%*%rt%*%t(jt);
      lpyt = -.5*log(det(rt)) - .5*t(yg[,i] - gg[,i] - t(ht)%*%(f[,i] + ft%*%mtm1))%*%rtinv%*%(yg[,i] - gg[,i] - t(ht)%*%(f[,i] + ft%*%mtm1));
      if (det(vt)<=0) {
         tt = zeros(kstate,kstate);
      } else {
         tt = t(chol(vt));
      }
      ot = omega[(kstate*(i-1)+1):(kstate*i),];
      mut = mu[(kstate*(i-1)+1):(kstate*i),];
      tempv = eye(kstate) + t(tt)%*%ot%*%tt;
      lpyt1n = -.5*log(det(tempv)) -.5*(t(mt)%*%ot%*%mt - 2*t(mut)%*%mt - t(mut - ot%*%mt)%*%tt%*%solve(tempv)%*%t(tt)%*%(mut - ot%*%mt));
      lprob[j,1] = log(kprior[j,1]) + lpyt1n + lpyt;
      if (i == 2) {
         lpy2n = lpyt1n + lpyt;
      }
   }
   pprob = exp(lprob)/sum(exp(lprob));
   tempv = runif(1);
   tempu = 0;
   for (j in 1:nvalk) {
      tempu = tempu + pprob[j];
      if (is.nan(tempu)) {
         next
      } else {
         if (tempu > tempv) {
            kdraw[i] = kvals[j,1];
            break
         }
      }
   }

   gat = kdraw[i]*sigv;
   rt = t(ht)%*%ft%*%vtm1%*%t(ft)%*%ht + t(ht)%*%crossprod(t(gat))%*%ht + crossprod(t(capg[((i-1)*p+1):(i*p),]));
   rtinv = solve(rt);
   jt = (ft%*%vtm1%*%t(ft)%*%ht + crossprod(t(gat))%*%ht)%*%rtinv;
   mtm1 = (eye(kstate) - jt%*%t(ht))%*%(f[,i] + ft%*%mtm1) + jt%*%(yg[,i] - gg[,i]);
   vtm1 = ft%*%vtm1%*%t(ft) + crossprod(t(gat)) - jt%*%rt%*%t(jt);
}
kdraw
}
gck1=function(yg,gg,hh,f, capf, sigv,kold,t,ex0,vx0,nvalk,kprior,kvals,p,kstate,sdraw) {
   # GCK.m modified to handle stochastic volatility case (i.e. mean and variance depend on sdraw)
   # Subroutine implements Gerlach, Carter and Kohn (JASA, 2000) and uses their notation
   # I have set it up so that their K enters only their Gamma matrix (i.e. state equation error covariance)
   # I have also set it up so that state equation has constant error covariance over time -- except for changepoints
   # The parts for Gamma matrix should be only part which are application specific
   # I have tried to follow their notation as far as possible and references to Lemmas and pages refer to this paper
   # One difference is that I have assumed error in state and measurement equation independent and, thus, used the slightly simpler
   # formulation of Giordani and Kohn in their appendices
   # Inputs: yg = p by t data matrix
   # gg = p by t matrix of det terms in measurement equation (usually set to zeros).
   # hh = tp by m matrix (note that this is transpose of GKCs definition)
   # f = m by t matrix of det terms in the state equation (usually set to zero in my work)
   # capf = tm by m matrix from state equation (set to identies for random walk evolution of states)
   # sigv is the standard deviation (or sigv*sigv') is the state equation error variance when regime change occurs (i.e. Kt=1)
   # sigv will be an m by m matrix 
   # kold is the previous draw of k, which is t by 1 (i.e. this code only allows for one k)
   # t = nunber of observations
   # kstate = dimensionality of state vector
   # ex0, vx0 = mean and variance for initial values for state equation (m by 1 and m by m)
   # nvalk = number of values k can take on -- usually 2 for 0/1
   # kprior = prior probabilities for each value for k = nvalk by 1 (this is okay for Bernoulli case, but in general may make this nvalk by t) 
   # kvals are the values k can take on -- usually 0/1
   mi = zeros(7,1);
   vi = zeros(7,1);
   
   # Set the values for the mixing distribution from KSC page 371
   vi[1,1]=5.79596;   vi[2,1] = 2.61369;  vi[3,1] = 5.17950;  vi[4,1] = 0.16735; vi[5,1] = 0.64009; vi[6,1] = 0.34023; vi[7,1] = 1.26261; 
   mi[1,1]=-10.12999; mi[2,1] = -3.97281; mi[3,1] = -8.56686; mi[4,1] = 2.77786; mi[5,1] = 0.61942; mi[6,1] = 1.79518; mi[7,1] = -1.08819; 
   
   capg = zeros(t*p,p);
   for (i in 1:t) {
      for (j in 1:p) {
         imix = sdraw[i,j];
         capg[((i-1)*p+j),j] = sqrt(vi[imix,1]);
         yg[j,i] = yg[j,i] - mi[imix,1] + 1.2704;  
      }
   }
   
   # GCK's Step 1 on page 821@
   mu = zeros(t*kstate,1);
   omega = zeros(t*kstate,kstate);
   gatplus1 = zeros(kstate,kstate);
   for (i in (t-1):1) {
      gatplus1 = sigv*kold[(i+1)];
      ftplus1 = capf[(kstate*i+1):(kstate*(i+1)),];
      cgtplus1 = capg[(i*p+1):((i+1)*p),];
      htplus1 = t(hh[(i*p+1):((i+1)*p),]);
      rtplus1 = (t(htplus1)%*%gatplus1)%*%(t(htplus1)%*%t(gatplus1)) + crossprod(t(cgtplus1));
      rtinv = solve(rtplus1);
      btplus1 = crossprod(t(gatplus1))%*%htplus1%*%rtinv;
      atplus1 = (eye(kstate) - btplus1%*%t(htplus1))%*%ftplus1;
      if (kold[(i+1)] == 0) {
         ctplus1 = zeros(kstate,kstate);
      } else {
         cct = gatplus1%*%(eye(kstate) - t(gatplus1)%*%htplus1%*%rtinv%*%t(htplus1)%*%gatplus1)%*%t(gatplus1);
         ctplus1 = t(chol(cct));
      }
      otplus1 = omega[(kstate*i+1):(kstate*(i+1)),];
      dtplus1 = t(ctplus1)%*%otplus1%*%ctplus1 + eye(kstate);
      omega[(kstate*(i-1)+1):(kstate*i),] = t(atplus1)%*%(otplus1 - otplus1%*%ctplus1%*%solve(dtplus1)%*%t(ctplus1)%*%otplus1)%*%atplus1 + t(ftplus1)%*%htplus1%*%rtinv%*%t(htplus1)%*%ftplus1;
      satplus1 = (eye(kstate) - btplus1%*%t(htplus1))%*%f[,(i+1)] - btplus1%*%gg[,(i+1)];
      mutplus1 = mu[(kstate*i+1):(kstate*(i+1)),];
      mu[(kstate*(i-1)+1):(kstate*i),] = t(atplus1)%*%(eye(kstate) - otplus1%*%ctplus1%*%solve(dtplus1)%*%t(ctplus1))%*%(mutplus1 - otplus1%*%(satplus1 + btplus1%*%yg[,(i+1)])) + t(ftplus1)%*%htplus1%*%rtinv%*%(yg[,(i+1)] - gg[,(i+1)] - t(htplus1)%*%f[,(i+1)]);  
   }
   
   # GCKs Step 2 on pages 821-822
   kdraw = kold;
   ht = t(hh[1:p,]);
   ft = capf[1:kstate,];
   gat = zeros(kstate,kstate);
   # Note: this specification implies no shift in first period -- sensible
   rt = t(ht)%*%ft%*%vx0%*%t(ft)%*%ht + t(ht)%*%gat%*%t(gat)%*%ht + crossprod(t(capg[1:p,]));
   rtinv = solve(rt);
   jt = (ft%*%vx0%*%t(ft)%*%ht + gat%*%t(gat)%*%ht)%*%rtinv;
   mtm1 = (eye(kstate) - jt%*%t(ht))%*%(f[,1] + ft%*%ex0) + jt%*%(yg[,1] - gg[,1]);
   vtm1 = ft%*%vx0%*%t(ft) + gat%*%t(gat) - jt%*%rt%*%t(jt);
   lprob = zeros(nvalk,1);
   for (i in 2:t) {
      ht = t(hh[((i-1)*p+1):(i*p),]);
      ft = capf[(kstate*(i-1)+1):(kstate*i),];
      for (j in 1:nvalk) {
         gat = kvals[j,1]*sigv;
         rt = t(ht)%*%ft%*%vtm1%*%t(ft)%*%ht + t(ht)%*%gat%*%t(gat)%*%ht + crossprod(t(capg[((i-1)*p+1):(i*p),]));
         rtinv = solve(rt);
         jt = (ft%*%vtm1%*%t(ft)%*%ht + gat%*%t(gat)%*%ht)%*%rtinv;
         mt = (eye(kstate) - jt%*%t(ht))%*%(f[,i] + ft%*%mtm1) + jt%*%(yg[,i] - gg[,i]);
         vt = ft%*%vtm1%*%t(ft) + gat%*%t(gat) - jt%*%rt%*%t(jt);
         lpyt = -.5*log(det(rt)) - .5*t(yg[,i] - gg[,i] - t(ht)%*%(f[,i] + ft%*%mtm1)) %*%rtinv%*%(yg[,i] - gg[,i] - t(ht)%*%(f[,i] + ft%*%mtm1));
         if (det(vt) <= 0) {
            tt = zeros(kstate,kstate);
         } else {
            tt = t(chol(vt));
         }
         ot = omega[(kstate*(i-1)+1):(kstate*i),];
         mut = mu[(kstate*(i-1)+1):(kstate*i),];
         tempv = eye(kstate) + t(tt)%*%ot%*%tt;
         lpyt1n = -.5*log(det(tempv)) -.5*(t(mt)%*%ot%*%mt - 2*t(mut)%*%mt - t(mut - ot%*%mt)%*%tt%*%solve(tempv)%*%t(tt)%*%(mut - ot%*%mt));
         lprob[j,1] = log(kprior[j,1]) + lpyt1n + lpyt; 
      }
      pprob = exp(lprob)/sum(exp(lprob));
      tempv = runif(1);
      tempu = 0;
      for (j in 1:nvalk) {
         tempu = tempu + pprob[j];
         if (is.nan(tempu)) {
            next
         } else {
            if (tempu > tempv) {
               kdraw[i] = kvals[j,1];
               break
            }
         }
      }
      gat = kdraw[i]*sigv;
      rt = t(ht)%*%ft%*%vtm1%*%t(ft)%*%ht + t(ht)%*%gat%*%t(gat)%*%ht + crossprod(t(capg[((i-1)*p+1):(i*p),]));
      rtinv = solve(rt);
      jt = (ft%*%vtm1%*%t(ft)%*%ht + gat%*%t(gat)%*%ht)%*%rtinv;
      mtm1 = (eye(kstate) - jt%*%t(ht))%*%(f[,i] + ft%*%mtm1) + jt%*%(yg[,i] - gg[,i]);
      vtm1 = ft%*%vtm1%*%t(ft) + gat%*%t(gat) - jt%*%rt%*%t(jt);
   }
   kdraw
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
   bdraw[t,] = t(btt) + rnorm(m)%*%chol(Vtt,pivot=T)
   
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
      bdraw[(t-i),] = t(bmean) + rnorm(m)%*%chol(bvar,pivot=T)
   }
   t(bdraw)
}

# TVP-VAR Time varying structural VAR with stochastic volatility
# --------------------------------------------------------------------------
# This code implements the model in Primiceri (2005).
# ***********************************************************
# The model is:
#__ __
# |Y(t)|= B(t) x |Y(t-1)|+u(t) 
#-- --
# 
# 
#with u(t)~N(0,H(t)), andA(t)' x H(t) x A(t) = Sigma(t)*Sigma(t),
# __
#| 100 ... 0|
#|a21(t)10 ... 0|
# A(t) =|a31(t)a32(t) 1 ... 0|
#|............... |
#|_ aN1(t)......aNN(t)1 _|
# 
# 
# and Sigma(t) = diag{s1(t), .... ,sn(t)}.
# **************************************************************
#NOTE: 
#There are references to equations of Primiceri, "Time Varying Structural Vector
#Autoregressions & Monetary Policy",(2005),Review of Economic Studies 72,821-852
#for your convenience. The definition of vectors/matrices is also based on this paper.
# --------------------------------------------------------------------------

#-----------------------------LOAD DATA------------------------------------
# Load Korobilis (2008) quarterly data
library("MASS")
library("matlab")
library("mvtnorm")
library("R.matlab")
Yraw = as.matrix(read.table("/Users/user/Dropbox/18. Korobilis Code/Bayesian Macroeconometrics/1 R Code for Bayesian VARs/R/Yraw.dat"))
ydata=Yraw

# Demean and standardize data
t2 = nrow(ydata);
# stdffr = std(ydata(:,3));
#for (i in 1:ncol(ydata)) {
#ydata[,i] = scale(ydata[,i],T,T)
#}

Y=ydata;
# Number of observations and dimension of X and Y
t=nrow(Y);
M=ncol(Y);

# Number of factors & lags:
tau = 40; # tau is the size of the training sample
p = M; # p is the dimensionality of Y
plag = 2; # plag is number of lags in the VAR part
numa = p*(p-1)/2; # numa is the number of elements of At
# ===================================| VAR EQUATION |==============================
# Generate lagged Y matrix. This will be part of the X matrix
# Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
ylag = embed(Y,plag+1)[,-c(1:M)]; # Y is [T x M]. ylag is [T x (Mp)]
ylag = ylag[(tau+1):nrow(ylag),];

m = p + plag*(p^2); # K is the number of elements in the state vector
# Create Z_t matrix.
Z = zeros((t-tau-plag)*p,m);
for (i in 1:(t-tau-plag)) {
   ztemp = eye(p);
   for (j in 1:plag) {
      xtemp = t(matrix(ylag[i,((j-1)*p+1):(j*p)]));
      xtemp = kronecker(diag(p),xtemp);
      ztemp = cbind(ztemp, xtemp);#ok<AGROW>
   }
   Z[((i-1)*p+1):(i*p),] = ztemp;
}

# Redefine FAVAR variables y
y = t(Y[(tau+plag+1):t,]);
#yearlab = yearlab[(tau+p+1):t];
# Time series observations
t=ncol(y);# t is now 215 - p - tau = 173

#----------------------------PRELIMINARIES---------------------------------
# Set some Gibbs - related preliminaries
nburn = 500;# Number of burn-in-draws
nrep = 500;# Number of replications
it_print = 100;#Print in the screen every "it_print"-th iteration

#========= PRIORS:
#========= PRIORS ON TRANSISION PROBS (Beta) 
# 1 - Least informative Beta prior
# a_prob = sqrt(t)*0.5;
# b_prob = sqrt(t)*0.5;

# 2 - Reference Beta prior
a_prob = 1;
b_prob = 1;

# 3 - Informative prior (few breaks)
# a_prob = 0.1;
# b_prob = 10;
ap_0 = a_prob*ones(2,1);
bp_0 = b_prob*ones(2,1);

# Implied prior "sample size" for state equations
t_0 = (ap_0/(ap_0 + bp_0))*t;

#=========PRIORS ON TIME-VARYING PARAMETERS AND THEIR COVARIANCES
# To set up training sample prior a la primiceri, use the following subroutine
#PRIOR = ts_prior(Y,tau,M,p)
#B_OLS=PRIOR$aols
#VB_OLS=PRIOR$vbar
#A_OLS=PRIOR$a0
#VA_OLS=PRIOR$a02mo
#sigma_OLS=PRIOR$ssig1

B_OLS = 0*ones(m,1);
A_OLS = 0*ones(numa,1);
VA_OLS = eye(numa);
VB_OLS = eye(m);
sigma_OLS = ones(p,1);

# Set some hyperparameters here (see page 831, } of section 4.1)
k_Q = 0.01;
k_S = 0.1;
k_W = 0.01;

# We need the sizes of some matrices as prior hyperparameters (see page
# 831 again, lines 2-3 and line 6)
sizeW = p; # Size of matrix W
sizeS = 1:p; # Size of matrix S

#-------- Now set prior means and variances (_prmean / _prvar)
# B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean = B_OLS;
B_0_prvar = 4*VB_OLS;

# A_0 ~ N(A_OLS, 4Var(A_OLS))
A_0_prmean = A_OLS;
A_0_prvar = 4*VA_OLS;

# log(sigma_0) ~ N(log(sigma_OLS),I_n)
sigma_prmean = sigma_OLS;
sigma_prvar = 4*eye(p);

# Note that for IW distribution I keep the _prmean/_prvar notation, but these are scale and shape parameters...
# Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean = ((k_Q)^2)*tau*(t/t_0[1,1])*VB_OLS;
Q_prvar = tau;

# W ~ IW(k2_W*(1+dimension(W))*I_n,(1+dimension(W)))
W_prmean = ((k_W)^2)*(t/t_0[1,1])*(1 + sizeW)*eye(p);
W_prvar = 1 + sizeW;

# S ~ IW(k2_S*(1+dimension(S)*Var(A_OLS),(1+dimension(S)))
S_prmean = cell(p-1,1);
S_prvar = zeros(p-1,1);

ind = 1;
for (ii in 2:p) {
   # S is block diagonal as in Primiceri (2005)
   S_prmean[[(ii-1)]] = ((k_S)^2)*(1 + sizeS[ii-1])*VA_OLS[(((ii-1)+(ii-3)*(ii-2)/2)):ind,((ii-1)+(ii-3)*(ii-2)/2):ind];
   S_prvar[(ii-1)] = 1 + sizeS[(ii-1)];
   ind = ind + ii;
}

# Parameters of the 7 component mixture approximation to a log(chi^2) density:
q_s = c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750);
m_s = c(-11.40039, -5.24321, -9.83726, 1.50746, -0.65098, 0.52478, -2.35859);
u2_s = c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261);

#========= INITIALIZE MATRICES:
# Specify covariance matrices for measurement and state equations
numa = p*(p-1)/2;
consQ = 0.0001;
consS = 0.0001;
consH = 0.01;
consW = 0.0001;
Qdraw = consQ*eye(m);
Qchol = sqrt(consQ)*eye(m);
Ht = kronecker(ones(t,1),consH*eye(p));
Htsd = kronecker(ones(t,1),sqrt(consH)*eye(p));
Sdraw = consS*eye(numa);
Sblockdraw = cell(p-1,1);
ijc = 1;

for (jj in 2:p) {
   Sblockdraw[[(jj-1)]] = Sdraw[((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc];
   ijc = ijc + jj;
}

Wdraw = consW*eye(p);
Bdraw = zeros(m,t);
Atdraw = zeros(numa,t);
Sigtdraw = zeros(p,t);
sigt = kronecker(ones(t,1),0.01*eye(p));
statedraw = 5*ones(t,p);
Zs = kronecker(ones(t,1),eye(p));
prw = zeros(numel(q_s),1);

kdraw = 1*ones(t,2);
pdraw = .5*ones(1,2);
kold = kdraw;
kmean = zeros(nrep,t,2);
kmax = zeros(t,2);
kvals = ones(2,1);
kvals[1,1] = 0;
kprior = .5*ones(2,1);

# Storage matrices for posteriors and stuff
B_postmean = zeros(nrep,m,t);
At_postmean = zeros(nrep,numa,t);
Sigt_postmean = zeros(nrep,p,t);
Qmean = zeros(nrep,m,m);
Smean = zeros(nrep,numa,numa);
Wmean = zeros(nrep,p,p);

sigmean = zeros(nrep,t,p);
cormean = zeros(nrep,t,numa);
kpmean = zeros(nrep,1,2);
sig2mo = zeros(nrep,t,p);
cor2mo = zeros(nrep,t,numa);
kp2mo = zeros(nrep,1,2);

# Model selection
llikmax = -9999e200;
llikmean = 0;
gelfday = 0;
lpymean = 0;

#----------------------------- END OF PRELIMINARIES ---------------------------

#====================================== START SAMPLING ========================================
#==============================================================================================
tic(); # This is just a timer
print('Number of iterations');
for (irep in 1:(nrep + nburn)) { # GIBBS iterations starts here
   # Print iterations
   if (irep%%it_print==0) {
      print(paste0(round(irep/(nrep+nburn)*100,2),"%"));
      toc();
   }
   
   # -----------------------------------------------------------------------------------------
   # STEP I: Sample B from p(B|y,A,Sigma,V) (Drawing coefficient states, pp. 844-845)
   # -----------------------------------------------------------------------------------------
   
   #----------------------------------------------------------------------
   # I.1: draw K1 index and related probabilities
   #----------------------------------------------------------------------
   ap = ap_0[1,1] + sum(kdraw[,1]);
   bp = bp_0[1,1] + t - sum(kdraw[,1]);

   pdrawa = rbeta(1,ap,bp);
   pdraw[1,1] = pdrawa;
   kprior[2,1] = pdrawa;
   kprior[1,1] = 1 - kprior[2,1];

   kdrawa = gck(y,zeros(p,t),Z,Htsd,zeros(m,t),kronecker(ones(t,1),eye(m)),t(Qchol),kold[,1],t,zeros(m,1),zeros(m,m),2,kprior,kvals,p,m);
   kold[,1] = kdraw[,1] = kdrawa;

   # I.2: draw Bt and keep only stationary draws
   #----------------------------------------------------------------------
   Bdrawc = carter_kohn(y,Z,Ht,Qdraw,m,p,t,B_0_prmean,B_0_prvar,kdraw[,1]);
   Bdraw = Bdrawc;
   
   Btemp = t(Bdraw[,2:t]) - t(Bdraw[,1:(t-1)]);
   sse_2 = zeros(m,m);
   for (i in 1:(t-1)) {
      sse_2 = sse_2 + crossprod(t(Btemp[i,]));
   }

   Qinv = solve(sse_2 + Q_prmean);
   Qinvdraw = rWishart(1,t+Q_prvar,Qinv)[,,1];
   Qdraw = solve(Qinvdraw);
   Qchol = chol(Qdraw);

   #-------------------------------------------------------------------------------------------
   # STEP II: Draw At from p(At|y,B,Sigma,V) (Drawing coefficient states, p. 845)
   #-------------------------------------------------------------------------------------------
   # Substract from the data y(t), the mean Z x B(t)
   yhat = zeros(p,t);
   for (i in 1:t) {
      yhat[,i] = y[,i] - Z[((i-1)*p+1):(i*p),]%*%Bdraw[,i];
   }

   # This part is more tricky, check Primiceri
   # Zc is a [p x p(p-1)/2] matrix defined in (A.2) page 845, Primiceri
   Zc = -t(yhat);
   sigma2temp = zeros(t,p);
   for (i in 1:t) {
      sigma2temp[i,] = diag(sigt[((i-1)*p+1):(i*p),]^2);
   }

   Atdraw = NULL;
   ind = 1;
   for (ii in 2:p) {
      # Draw each block of A(t)
      Atblockdraw = carter_kohn(t(matrix(yhat[ii,])),Zc[,1:(ii-1),drop=F],sigma2temp[,ii,drop=F],Sblockdraw[[(ii-1)]],sizeS[(ii-1)],1,t,A_0_prmean[((ii-1)+(ii-3)*(ii-2)/2):ind,,drop=F],A_0_prvar[((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind,drop=F],ones(t,1));
      Atdraw = rbind(Atdraw, Atblockdraw); #ok<AGROW> # Atdraw is the final matrix of draws of A(t)
      ind = ind + ii;
   }
   
   #=====| Draw S, the covariance of A(t) (from iWishart)
   Attemp = t(Atdraw[,2:t]) - t(Atdraw[,1:(t-1)])
   sse_2 = zeros(numa,numa)
   for (i in 1:(t-1)) {
      sse_2 = sse_2 + crossprod(t(Attemp[i,]));
   }

   ijc = 1;
   for (jj in 2:p) {
      Sinv = solve(sse_2[((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc] + S_prmean[[(jj-1)]]);
      Sinvblockdraw = rWishart(1,t+S_prvar[(jj-1)],Sinv)[,,1];
      Sblockdraw[[(jj-1)]] = solve(Sinvblockdraw);
      ijc = ijc + jj;
   }
   
   #------------------------------------------------------------------------------------------
   # STEP III: Draw diagonal VAR covariance matrix "H_t" elements
   #------------------------------------------------------------------------------------------
   capAt = zeros(p*t,p);
   for (i in 1:t) {
      capatemp = eye(p);
      aatemp = Atdraw[,i];
      ic=1;
      for (j in 2:p) {
         capatemp[j,1:(j-1)] = t(aatemp[ic:(ic+j-2)]);
         ic = ic + j - 1;
      }
      capAt[((i-1)*p+1):(i*p),] = capatemp;
   }
   
   # III.1: draw "H_t"s' elements from Normal distribution
   #----------------------------------------------------------------------
   y2 = NULL;
   for (i in 1:t) {
      ytemps = capAt[((i-1)*p+1):(i*p),]%*%yhat[,i];
      y2 = cbind(y2,ytemps^2); #ok<AGROW>
   }
   
   yss = zeros(t,p);
   for (i in 1:p) {
      yss[,i] = log(y2[i,] + 0.001);
   }

   # III.1: draw K2 index and related probabilities
   #----------------------------------------------------------------------
   ap = ap_0[2,1] + sum(kdraw[,2]);
   bp = bp_0[2,1] + t - sum(kdraw[,2]);
   
   pdrawa = rbeta(1,ap,bp);
   pdraw[1,2] = pdrawa;
   kprior[2,1] = pdrawa;
   kprior[1,1] = 1 - kprior[2,1];
   Wchol = t(chol(Wdraw));
   
   kdrawa = gck1(t(yss),zeros(p,t),Zs,zeros(p,t),kronecker(ones(t,1),eye(p)),Wchol,kold[,2],t,zeros(p,1),zeros(p,p),2,kprior,kvals,p,p,statedraw);
   kdraw[,2] = kdrawa;
   kold[,2] = kdraw[,2];


   # III.2: draw statedraw (chi square approximation mixture component)
   #----------------------------------------------------------------------
   for (jj in 1:p) {
      for (i in 1:t) {
         for (j in 1:length(m_s)) {
            temp1= (1/sqrt(2*pi*u2_s[j]))*exp(-.5*(((yss[i,jj] - Sigtdraw[jj,i] - m_s[j] + 1.2704)^2)/u2_s[j]));
            prw[j,1] = q_s[j]*temp1;
         }
         prw = prw/sum(prw);
         cprw = cumsum(prw);
         trand = runif(1);
         if (trand < cprw[1]) {
           imix=1
         } else if (trand < cprw[2]) {
            imix=2;
         } else if (trand < cprw[3]) {
            imix=3;
         } else if (trand < cprw[4]) {
            imix=4;
         } else if (trand < cprw[5]) {
            imix=5;
         } else if (trand < cprw[6]) {
            imix=6;
         } else { 
            imix=7; 
         }
         statedraw[i,jj]=imix;
      }
   }

   # III.3: draw "H_t"s' elements from Normal distribution
   #---------------------------------------------------------------------- 
   #First draw volatilities conditional on sdraw
   vart = zeros(t*p,p);
   yss1 = zeros(t,p);
   for (i in 1:t) {
      for (j in 1:p) {
         imix = statedraw[i,j];
         vart[((i-1)*p+j),j] = u2_s[imix];
         yss1[i,j] = yss[i,j] - m_s[imix] + 1.2704;
      }
   }
   Sigtdraw = carter_kohn(t(yss1),Zs,vart,Wdraw,p,p,t,sigma_prmean,sigma_prvar,kdraw[,2]); 

   sigtemp = eye(p);
   sigt = zeros(p*t,p);
   for (i in 1:t) {
      for (j in 1:p) {
         sigtemp[j,j] = exp(.5*Sigtdraw[j,i]);
      }
      sigt[((i-1)*p+1):(i*p),] = sigtemp;
   }
   
   Sigttemp = t(Sigtdraw[,2:t]) - t(Sigtdraw[,1:(t-1)]);
   sse_2 = zeros(p,p);
   for (i in 1:(t-1)) {
      sse_2 = sse_2 + crossprod(t(Sigttemp[i,]));
   }
   Winv = solve(sse_2 + W_prmean);
   Winvdraw = rWishart(1,t+W_prvar,Winv)[,,1];
   Wdraw = solve(Winvdraw);

   Ht = zeros(p*t,p);
   Htsd = zeros(p*t,p);
   for (i in 1:t) {
      inva = solve(capAt[((i-1)*p+1):(i*p),]);
      stem = sigt[((i-1)*p+1):(i*p),];
      Hsd = inva%*%stem;
      Hdraw = crossprod(t(Hsd));
      Ht[((i-1)*p+1):(i*p),] = Hdraw;
      Htsd[((i-1)*p+1):(i*p),] = Hsd;
   }

   #----------------------------SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES -----------------
   if (irep > nburn) {
      B_postmean[(irep-nburn),,] = Bdraw;
      At_postmean[(irep-nburn),,] = Atdraw;
      Sigt_postmean[(irep-nburn),,] = Sigtdraw;
      Qmean[(irep-nburn),,] = Qdraw;
      ikc = 1;
      for (kk in 2:p) { 
         Sdraw[((kk-1)+(kk-3)*(kk-2)/2):ikc,((kk-1)+(kk-3)*(kk-2)/2):ikc]=Sblockdraw[[(kk-1)]]
         ikc = ikc + kk;
      }
      Smean[(irep-nburn),,] = Sdraw;
      Wmean[(irep-nburn),,] = Wdraw;
      kmean[(irep-nburn),,] = kdraw;
      kpmean[(irep-nburn),,] = pdraw;
      kp2mo[(irep-nburn),,] = pdraw^2;
      
      stemp6 = zeros(p,1);
      stemp5 = NULL;
      stemp7 = NULL;
      for (i in 1:t) {
         std = sqrt(diag(Ht[((i-1)*p+1):(i*p),]));
         stemp8 = Ht[((i-1)*p+1):(i*p),]/(std%*%t(std));
         stemp7a = NULL;
         ic = 1;
         for (j in 1:p) {
            if (j>1) {
               stemp7a = c(stemp7a,stemp8[j,1:ic]) #ok<AGROW>
               ic = ic+1;
            }
            stemp6[j,1] = sqrt(Ht[((i-1)*p+j),j]);
         }
         stemp5 = rbind(stemp5, t(stemp6)); #ok<AGROW>
         stemp7 = rbind(stemp7, t(stemp7a)); #ok<AGROW>
      }
      sigmean[(irep-nburn),,] = stemp5;
      cormean[(irep-nburn),,] = stemp7; 
      sig2mo[(irep-nburn),,] = stemp5^2;
      cor2mo[(irep-nburn),,] = stemp7^2;
   } # } saving after burn-in results 
} #} main Gibbs loop (for irep = 1:nrep+nburn)

toc(); # Stop timer and print total time
#=============================GIBBS SAMPLER ENDS HERE==================================

B_postmean1 = apply(B_postmean,2:3,mean);
B_05 = apply(B_postmean,2:3,quantile,0.05);
B_50 = apply(B_postmean,2:3,quantile,0.50);
B_95 = apply(B_postmean,2:3,quantile,0.95);
Sigt_postmean1 = apply(Sigt_postmean,2:3,mean);
Sigt_05 = apply(Sigt_postmean,2:3,quantile,0.05);
Sigt_50 = apply(Sigt_postmean,2:3,quantile,0.50);
Sigt_95 = apply(Sigt_postmean,2:3,quantile,0.95);

At_postmean1 = apply(At_postmean,2:3,mean);
Qmean1 = apply(Qmean,2:3,mean);
Smean1 = apply(Smean,2:3,mean);
Wmean1 = apply(Wmean,2:3,mean);
kmean1 = apply(kmean,2:3,mean);
kpmean1 = apply(kpmean,2:3,mean);
kp2mo1 = apply(kp2mo,2:3,mean);
sigmean1 = apply(sigmean,2:3,mean);
cormean1 = apply(cormean,2:3,mean);
sig2mo1 = apply(sig2mo,2:3,mean);
cor2mo1 = apply(cor2mo,2:3,mean)


par(mfcol = c(7,3), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:dim(B_postmean1)[1]) {
   plot(B_postmean1[i,],type="l",xaxs="i",las=1,xlab="",ylab="",tck=0.02,ylim=c(min(B_05[i,]),max(B_95[i,])))
   lines(B_50[i,],col="steelblue4",lty=2)
   lines(B_05[i,],col="steelblue4")
   lines(B_95[i,],col="steelblue4")
}

par(mfcol = c(M,1), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:M) {
   plot(Sigt_postmean1[i,],type="l",xaxs="i",las=1,ylim=c(min(Sigt_05[i,]),max(Sigt_95[i,])),xlab="",ylab="",tck=0.02)
   lines(Sigt_50[i,],col="steelblue4",lty=2)
   lines(Sigt_05[i,],col="steelblue4")
   lines(Sigt_95[i,],col="steelblue4")
}

### END
