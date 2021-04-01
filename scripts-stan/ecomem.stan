functions {
  real spherical(real x, real phi) {
    real g;
    
    if (x <= 0) {
      reject("x negative spherical function undefined; found x = ", x);
    }
    else if (x > phi) g=0;
    else {
        g = (1 - 3.0/2.0 * x/phi + 1.0/2.0 * (x/phi)^3); 
    }
    
    return g;
  }
}

data { 
  int<lower=0> N_years; 
  int<lower=0> N_sites; 
  int<lower=0> lag;
  
  vector[N_years] Y[N_sites];
  vector[N_years] X[N_sites];
  vector[N_years] d[N_sites];
} 

parameters { 
  real<lower=1e-6> tau; 
  real<lower=1e-6> sigma;
  //vector<lower=1e-6>[N_sites] sigma;
  real<lower=1e-6> sigma_alpha;
  vector[N_sites] u;
  
  real beta;
  
  real gamma_0;
  real gamma_1;
} 

transformed parameters { 
  real phi;
  //matrix[N_years, N_sites] z;
  //matrix[N_years, N_sites] alpha;
  vector[N_years] z[N_sites];
  vector[N_years] alpha[N_sites];
  
  phi = 1/tau;
  
  for (site in 1:N_sites) {
    z[site] = rep_vector(0, N_years);
    for (year in (lag+1):N_years) {
      for (l in 1:lag){
          z[site][year] = z[site][year] + d[site][year-l] * spherical(year-l, phi); 
      }
    }
    
  }
  
  for (site in 1:N_sites) {
    alpha[site][(lag+1):N_years] = gamma_0 + z[site][(lag+1):N_years] * gamma_1;
  }
  
} 

model { 
  // priors
  tau ~ uniform(1.0/15.0, 1); 
  sigma_alpha ~ cauchy(0, 5);
  sigma ~ cauchy(0, 5);
  
  beta ~ normal(0, 10);
  
  gamma_0 ~ normal(0, 1e5);
  gamma_1 ~ normal(0, 1e5);
  
  u ~ normal(0, sigma_alpha);

  for (site in 1:N_sites){
    Y[site][(lag+1):N_years]  ~ normal(alpha[site][(lag+1):N_years] + X[site][(lag+1):N_years] * beta + u[site], sigma); 
    //Y[site][(lag+1):N_years]  ~ normal(alpha[site][(lag+1):N_years] + X[site][(lag+1):N_years] * beta, sigma[site]); 
  }
} 
