data{
  int<lower = 1> N;
  int<lower = 0> present[N];
  int<lower = 0> spent[N];
  int<lower = 0> S;
  int species[N];
}
parameters{
  real<lower = 0,upper = 1> p[S];
}
model{
  for(i in 1:N){
    spent[i] ~ binomial(present[i], p[species[i]]);
  }
}
