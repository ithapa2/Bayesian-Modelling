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
  int<lower=0> T;
  int<lower=0> T_clim;
  int<lower=0> clim_start;
  int<lower=0> N_sites; 
  int<lower=0> lag;
  int<lower=0> N_covars;
  int<lower=0> n_basis;
  int<lower=0> n_knots;
  
  vector[T] Y[N_sites];
//  vector[N_years] X[N_sites];
  matrix[T_clim, N_covars] X[N_sites];
  vector[T_clim] d[N_sites];
  vector[T] fire;
  
  matrix[n_knots, n_basis] B;
  matrix[n_basis, n_basis] S_inv;
} 
transformed data {
    real sigma_eta;
      
    sigma_eta = 0.1;
}

parameters { 
  //real<lower=1e-6> tau; 
  real<lower=1.0/15.0, upper=1> tau_fire; 
  real<lower=1e-6> sigma;
  //vector<lower=1e-6>[N_sites] sigma;
  real<lower=1e-6> sigma_alpha;
  vector[T] u[N_sites];
  
  vector[N_covars] beta;

  real gamma0;
  real gamma1;
  real gamma2;
  //vector[N_sites] gamma_1;
  
  vector[n_basis] eta;
  //real<lower=1e-6> sigma_eta;
} 

transformed parameters { 
  //real phi;
  real phi_fire;
  //matrix[N_years, N_sites] z;
  //matrix[N_years, N_sites] alpha;
  vector[T_clim] z[N_sites];
  vector[T] f;
  vector[T] alpha[N_sites];
  
  // phi = 1/tau;
  // 
  // for (site in 1:N_sites) {
  //   z[site] = rep_vector(0, N_years);
  //   for (year in (lag+1):N_years) {
  //     for (l in 1:lag){
  //         z[site][year] = z[site][year] + d[site][year-l] * spherical(year-l, phi); 
  //     }
  //   }
  // }
  
  vector[lag+1] w;
  real w_sum;
  
  //for (l in 1:(lag+1)){
  w = B * eta;
  //}
  
  w = exp(w);
  w_sum = sum(w);
  w = w/w_sum;
  
  
  for (site in 1:N_sites) {
    z[site] = rep_vector(0, T_clim);
    for (year in (lag+1):T_clim) {
      for (l in 0:lag){
          z[site][year] = z[site][year] + d[site][year-l] * w[l+1]; 
      }
    }
  }
  
  
  phi_fire = 1/tau_fire;
  
    f = rep_vector(0, T);
    for (year in 1:T) {
      for (prior_year in 1:(year-1)){
        if (fire[year] == 1){
          f[year] = f[year] + spherical(year-prior_year, phi_fire);
        }
      }
    }
  
  for (site in 1:N_sites) {
    for (t in 1:(clim_start+lag-1)){
      //alpha[site][1:(clim_start+lag-1)] = gamma0 + f[1:(clim_start+lag-1)] * gamma2;
      //alpha[site][(clim_start+lag):T] = gamma0 + z[site][(1+lag):T_clim] * gamma1 + f[(clim_start+lag):T] * gamma2;
      alpha[site][t] = gamma0 + f[t] * gamma2;
    }
    for (t in (clim_start+lag):T){
      alpha[site][t] = gamma0 + z[site][t-clim_start+1] * gamma1 + f[t] * gamma2;
    }
  }

} 

model { 
  
  // priors
  //tau ~ uniform(1.0/15.0, 1); 
  tau_fire ~ uniform(1.0/15.0, 1); 
  sigma_alpha ~ cauchy(1e-6, 5);
  sigma ~ cauchy(1e-6, 5);
  
  beta ~ normal(0, 10);
  
  gamma0 ~ normal(0, 1e5);
  gamma1 ~ normal(0, 1e5);
  gamma2 ~ normal(0, 1e5);
  
  for (site in 1:N_sites){
    u[site] ~ normal(0, sigma_alpha);
  }
  
  
  eta ~ multi_normal(rep_vector(0, n_basis), sigma_eta^2 * S_inv);

//  for (site in 1:N_sites){
//    Y[site][(lag+1):N_years]  ~ normal(alpha[site][(lag+1):N_years] + X[site][(lag+1):N_years] * beta + u[site][(lag+1):N_years], sigma); 
//    //Y[site][(lag+1):N_years]  ~ normal(alpha[site][(lag+1):N_years] + X[site][(lag+1):N_years] * beta + u[site][(lag+1):N_years], sigma); 
//  }
  
  for (site in 1:N_sites) {
    for (t in 1:(clim_start+lag-1)){
       Y[site][t]  ~ normal(alpha[site][t] + u[site][t], sigma); 
    }
    for (t in (clim_start+lag):T){
      Y[site][t]  ~ normal(alpha[site][t] + X[site][t-clim_start+1] * beta + u[site][t], sigma); 
      //alpha[site][t] = gamma0 + z[site][t-clim_start+1] * gamma1 + f[t] * gamma2;
    }
  }
  
} 

