data {
  int<lower=0> N;         // number of previous patients
  int y[N];               // binary response
  int x[N];               // dose levels
  int J;                  // number of doses
  real doses[J];          // pseudo doses

  int<lower=0> N0;        // number of previous patients NO
  int y0[N0];             // binary response
  int x0[N0];             // dose levels
  int J0;                 // number of doses
  real doses0[J0];        // pseudo doses
  
  real<lower=0, upper=1> alpha0;
  real a;                 // crm intercept
}
parameters {
  real beta; 
}
transformed parameters{
  real pstim[J];
  for (j in 1:J) {
    pstim[j] = 1 / (1 + exp(-(a + exp(beta)*doses[j])));
  }
}
model {
  real p[N];          // probabilities
  vector[N] z;        // logistic transformation
  real p0[N0];          // probabilities
  vector[N0] z0;        // logistic transformation
  
  
  for (n in 1:N){
  z[n] = a + doses[x[n]]*exp(beta);
  p[n] = 1 / (1 + exp(-z[n]));
  }
  y ~ bernoulli(p);

  // prior

  for (i in 1:N0){
  z0[i] = a + doses0[x0[i]]*exp(beta);
  p0[i] = 1 / (1 + exp(-z0[i]));
  target += bernoulli_lpmf(y0[i] | p0[i])*alpha0;
  }
  beta ~ normal(0, sqrt(1.34));       // non-info prior
}
