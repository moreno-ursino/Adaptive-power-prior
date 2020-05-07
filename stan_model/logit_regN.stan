data {
  int<lower=0> N;         // number of previous patients
  int y[N];               // binary response
  int x[N];               // dose levels
  int J;                  // number of doses
  real doses[J];          // pseudo doses

  real<lower=0, upper=1> NN0;
  real a;                 // crm intercept
}
parameters {
  real<lower=-5, upper=5> beta; 
}
model {
  real p[N];          // probabilities
  vector[N] z;        // logistic transformation
     // logistic transformation

  for (n in 1:N){
  z[n] = a + doses[x[n]]*exp(beta);
  p[n] = 1 / (1 + exp(-z[n]));
  target += bernoulli_lpmf(y[n] | p[n])*NN0;
  }
}
