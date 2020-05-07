##############################################################################################################
# B-CRM (Bridging Continuous Reassessment Method)   based on 	
#                Suyu Liu, Sept 1, 2014
#

MCMCprobit <- function (x, y, nburn=1000, niter=10000)
{ 
  truncNorm = function(lower,upper,mean,sd){
    problow = pnorm((lower-mean)/sd)
    probup = pnorm((upper-mean)/sd)
    u = runif(length(problow))
    p = problow + u*(probup-problow)
    p[p<exp(-37)]=exp(-37)
    p[p>1-exp(-37)]=1-exp(-37)
    z = qnorm(p,mean,sd)
    z
  }
  
  # center the covariates
  x.bar = mean(x)
  x.c = x-x.bar;
  Sxx=sum(x.c^2);
  n=length(y);
  
  # define normal prior for beta0 and uniform prior for beta1
  prior.mu0 = -1.645;
  prior.var0 = 4;
  beta1.L = 0;
  beta1.R = 4;
  
  y1.ind = (y==1);
  y0.ind = (y==0);
  
  z = numeric(n);
  Beta0 = numeric(nburn+niter);
  Beta1 = numeric(nburn+niter);
  beta0 = prior.mu0;
  beta1 = 1;
  
  for(i in 1:(nburn+niter))
  {
    # draw latent variable Z
    z[y1.ind] = truncNorm(0, 10^6, beta0 + beta1*x.c[y1.ind], 1);
    z[y0.ind] = truncNorm(-10^6, 0, beta0 + beta1*x.c[y0.ind], 1);
    
    # draw beta0 and beta1
    beta0.hat = mean(z);
    var.beta0 = 1/n;
    beta0 = rnorm(1, (prior.mu0*var.beta0+beta0.hat*prior.var0)/(var.beta0+prior.var0), sqrt(var.beta0*prior.var0/(var.beta0+prior.var0)));
    
    if(Sxx!=0)  # there are more than one unique values of X
    {
      Szx=sum((z-mean(z))*x);
      beta1.hat = Szx/Sxx;
      var.beta1 = 1/Sxx;
      beta1 = truncNorm(beta1.L, beta1.R, beta1.hat, sqrt(var.beta1));
    }
    else { beta1 = runif(1, beta1.L, beta1.R); } #draw from prior as slope is not identifiable
    
    Beta0[i]=beta0-beta1*x.bar;
    Beta1[i]=beta1;
  }
  return(c(mean(Beta0[(nburn+1):(nburn+niter)]), mean(Beta1[(nburn+1):(nburn+niter)])))
}

## generate skeletons based on the landmark data y (i.e., the number of toxicity at each dose) 
## and n (i.e., the number of patients at each dose). 
get.skel <- function(y, n)
{
  ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
  pava <- function (x, wt = rep(1, length(x))) 
  {
    n <- length(x)
    if (n <= 1) 
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol))) 
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }
  
  ndose = length(y);
  x = seq(0, 1, 1/(ndose-1)); 
  y.bin = NULL; ## change data into binary format to use MCMClogit
  x.bin = NULL;
  for(i in 1:ndose) 
  {	
    y.bin = c(y.bin, rep(1,y[i]), rep(0, n[i]-y[i]));
    x.bin = c(x.bin, rep(x[i], n[i]));
  }
  
  ###################### estimate dose-toxicity curve ###################
  ## parametric regression, Bayeisan probit model
  probitfit = MCMCprobit(x.bin, y.bin);
  beta0 = probitfit[1];
  beta1 = probitfit[2];
  pi.para = pnorm(beta0 + beta1*x);
  
  ## isotonic regression
  phat = (y+0.005)/(n+0.01); 
  phat.var = (y+0.005)*(n-y+0.005)/((n+0.01)^2*(n+0.01+1));
  pi.isoto = pava(phat, wt=1/phat.var);  ## perform the isotonic transformation using PAVA
  
  ## weighted estimate
  r = (pi.para^y*(1-pi.para)^(n-y))/(pi.isoto^y*(1-pi.isoto)^(n-y));
  w = r/(1+r);
  w[n==0]=1;
  pi.w = w*pi.para + (1-w)*pi.isoto;
  
  pi.w = pava(pi.w);
  return(pi.w);
}

sumres=table(factor(D0$x0, levels = 1:6), D0$y0)
y.land=as.vector(sumres[,2])
n.land= as.vector(apply(sumres, 1, sum))
pi.w = get.skel(y.land, n.land);
pi.w = pi.w + c(0,0,0,0.01,0,0.01)
p1 = pi.w;
ndose=length(p1)
p2 = c(pi.w[2:ndose], (pi.w[ndose]+1)/2);  
p3 = c(pi.w[1]/2, pi.w[1:(ndose-1)]); 
skeletons = rbind(p1, p2, p3);


# p1 = c(0.05, 0.07, 0.2, 0.4, 0.5, 0.55)
# ndose=length(p1)
# p2 = c(p1[2:ndose], (p1[ndose]+1)/2);  
# p3 = c(p1[1]/2, p1[1:(ndose-1)]); 
# skeletons = rbind(p1, p2, p3);


## CRM using Bayesian model averaging 
bma.crm <- function (skeletons, target, y, n, BMA=1)
{ 
  # if a single skeleton is inputed as a vector, convert it to a matrix
  if(is.vector(skeletons)) skeletons=t(as.matrix(skeletons));
  
  nskel = nrow(skeletons);
  mprior = rep(1/nskel, nskel); 
  ndose = ncol(skeletons)# prior for each model formed the skeleton
  ptox.hat = numeric(ndose);
  # posterior = likelihood x prior
  posterior <- function(alpha, p, y, n)
  {
    sigma2 = 2;
    lik=1;
    for(j in 1:length(p))
    {
      pj = p[j]^(exp(alpha));
      lik = lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
  }
  
  # the posterior mean of ptox
  posttoxf <- function(alpha, p, y, n, j) { p[j]^(exp(alpha))*posterior(alpha, p, y, n); }
  
  ### run a trial 	
  TOXSTOP = 2; ## toxicity stopptoxng boundary
   #number of toxicity at each dose level y
  #number of treated patients at each dose level n
  #stop=0; #indicate if trial stops early
  
    marginal = rep(0, nskel);
    for(k in 1:nskel)
    {
      marginal[k] = integrate(posterior,lower=-Inf,upper=Inf, skeletons[k,], y, n)$value;
    }
    
    postprob = (marginal*mprior)/sum(marginal*mprior);
    
    if(BMA==1)  ### model averaging
    {
      pj.overtox = rep(0, nskel);
      for(k in 1:nskel)
      {
        pj.overtox[k] = integrate(posterior,lower=-Inf,upper=log(log(target)/log(skeletons[k,1])), skeletons[k,], y, n)$value/marginal[k];  
      }
      p.overtox = sum(postprob*pj.overtox);
      if(p.overtox>TOXSTOP) { stop=1; break;}
      
      # calculate posterior mean of toxicity probability at each dose leavel
      ptoxj.hat = rep(0, nskel);
      for(j in 1:ndose){
        for(k in 1:nskel){ 
          ptoxj.hat[k] = integrate(posttoxf,lower=-Inf,upper=Inf, skeletons[k,], y, n, j)$value/marginal[k]; 
        }
        ptox.hat[j] = sum(postprob*ptoxj.hat);
      }
    }	

    
    diff = abs(ptox.hat-target);
    dose.best = min(which(diff==min(diff)));
    
  return(list(newDose=dose.best, pstim = ptox.hat, parameters = 0,
    gamma=0, alphaPP = 0))
}

# target = 0.3;
# pi.w = get.skel(y.land, n.land);
# p1 = pi.w;
# 
# 
# 
# ## An example
# #n.land <- c(3, 3, 3, 6, 3, 0);
# #y.land <- c(0, 0, 0, 1, 2, 0);
# #mtd.land <- 4
# 
# 
# n.land= c(as.vector(apply(table(D0$x0,D0$y0), 1, sum)),0)
# y.land= c(as.vector(table(D0$x0,D0$y0)[,2]),0)
# mtd.land <- 3



