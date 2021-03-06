data{
  int<lower = 1> N;
  int Fecundity[N];
  matrix[N,7] ModelMatrix;
  int<lower = 1>P;
  int Plot[N];
}
parameters{
  real epsilon[P];
  real sigma_0;
  real<lower = 0> lambda;
  real<lower = 0> alpha_intra;
  real<lower = 0> alpha_arca;
  real<lower = 0> alpha_plde;
  real<lower = 0> alpha_trcy;
  real<lower = 0> alpha_vero;
  real<lower = 0> alpha_peai;
  real<lower = 0> alpha_other;
}

transformed parameters{
  real<lower = 0> sigma;
  sigma = exp(sigma_0);
}

model{
  // create a vector of predictions
  vector[N] F_hat;
  vector[N] F_hat2;

  // set priors
  sigma_0 ~ normal(0, 1000);
  epsilon ~ gamma(sigma, sigma);
  lambda ~ gamma(0.001, 0.001);
  alpha_intra ~ normal(0, 1);
  alpha_arca ~ normal(0, 1);
  alpha_plde ~ normal(0, 1);
  alpha_trcy ~ normal(0, 1);
  alpha_vero ~ normal(0, 1);
  alpha_peai ~ normal(0, 1);
  alpha_other ~ normal(0, 1);

  // implement the biological model
  for(i in 1:N){
   F_hat[i] = lambda / (1 + alpha_intra * ModelMatrix[i,1] + alpha_arca * ModelMatrix[i,2] + alpha_plde * ModelMatrix[i,3] + alpha_trcy * ModelMatrix[i,4] + alpha_vero * ModelMatrix[i,5] + alpha_peai*ModelMatrix[i,6] + alpha_other * ModelMatrix[i,7]);
   F_hat2[i] = F_hat[i]*epsilon[Plot[i]];
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat2);
}
