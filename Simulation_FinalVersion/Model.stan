data {
  int N;
  int K;
  real Z_cov[N];
  int Z_resp[N];
  real W1[N, K];
  real W2[N, K];
  real X1[N];
  real X2[N];
  real alpha;
}

parameters {
  real phi1[2];
  real phi2[3];
  real beta0;
  real betaX1;
  real betaX2;
  real betaZ;
  real X1_new[N];
  real X2_star_new[N];
}

model {
  for (n in 1:N)
    X1[n] ~ normal(Z_cov[n]*phi1[1] + phi1[2], 1);
    
  for (n in 1:N)  
    X2[n] ~ normal(Z_cov[n]*phi1[1] + phi1[2], 1);
    
  for (n in 1:N)
    X1_new[n] ~ normal(Z_cov[n]*phi1[1] + phi1[2], 1);
  
  for (n in 1:N)
    X2_star_new[n] ~ normal(Z_cov[n]*phi1[1] + phi1[2], 1);
    
  for (n in 1:N)
    Z_resp[n] ~ bernoulli(Phi(beta0 + betaX1*X1[n] + betaX2*X2[n] + betaZ*Z_cov[n]));
  
  for (n in 1:N)
    for (k in 1: K)
      W1[n,k] ~ normal(X1_new[n],1);

 for (n in 1:N)
    for (k in 1: K)
      W2[n,k] ~ normal(X2_star_new[n],1);
}

