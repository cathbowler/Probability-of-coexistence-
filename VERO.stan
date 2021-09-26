data{
  int<lower = 1> N;
  int Fecundity[N];
  matrix[N,6] ModelMatrix;
  int<lower = 1>P;
  int Plot[N];
}
parameters{
  real epsilon[P];
  real<lower = 0> sigma;
  real<lower = 0> lambda;
  real<lower = 0> alpha_intra;
  real<lower = 0> alpha_hygl;
  real<lower = 0> alpha_poca;
  real<lower = 0> alpha_trcy;
  real<lower = 0> alpha_peai;
  real<lower = 0> alpha_other;
}
model{
  // create a vector of predictions
  vector[N] F_hat;
  vector[N] F_hat2;

  // set priors
  sigma ~ gamma(0.001, 0.001);
  epsilon ~ gamma(sigma, sigma);
  lambda ~ gamma(0.001, 0.001);
  alpha_intra ~ normal(0, 1);
  alpha_hygl ~ normal(0,1);
  alpha_poca ~ normal(0, 1);
  alpha_trcy ~ normal(0, 1);
  alpha_peai ~ normal(0, 1);
  alpha_other ~ normal(0, 1);

  // implement the biological model
  for(i in 1:N){
   F_hat[i] = lambda / (1 + alpha_intra * ModelMatrix[i,1] + alpha_hygl * ModelMatrix[i,2] + alpha_poca * ModelMatrix[i,3] + alpha_trcy * ModelMatrix[i,4] + alpha_peai * ModelMatrix[i,5] + alpha_other * ModelMatrix[i,6]);
   F_hat2[i] = F_hat[i]*epsilon[Plot[i]];
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat2);
}
