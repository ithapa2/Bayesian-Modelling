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
  //int<lower=0> N_covars;
  //int<lower=0> n_basis;
  //int<lower=0> n_knots;
  
  //int<lower=0, upper=1> is_fire;
  
  vector[N_years] Y[N_sites];
  //matrix[N_years, N_covars] X[N_sites];
  //vector[N_years] d[N_sites];
  vector[N_years] fire;
  
  //matrix[n_knots, n_basis] B;
  //matrix[n_basis, n_basis] S_inv;
} 
transformed data {
    real sigma_eta;
      
    sigma_eta = 0.1;
}

parameters { 
  //real<lower=1e-6> tau; 
  real<lower=1.0/15.0, upper=1> tau_fire; 
  //real<lower=1e-6> sigma;
  vector<lower=1e-6>[N_sites] sigma;
  real<lower=1e-6> sigma_alpha;
  //vector<lower=1e-6>[N_years] u[N_sites];
  vector[N_years] u[N_sites];
  
  //vector[N_covars] beta;

  real gamma0;
  //real gamma1;
  real gamma2;
  //vector[N_sites] gamma_1;
  
  //vector[n_basis] eta;
  //real<lower=1e-6> sigma_eta;
} 

transformed parameters { 
  //real phi;
  real phi_fire;
  //matrix[N_years, N_sites] z;
  //matrix[N_years, N_sites] alpha;
  //vector[N_years] z[N_sites];
  vector[N_years] f;
  vector[N_years] alpha[N_sites];
  

    phi_fire = 1/tau_fire;
  
    f = rep_vector(0, N_years);
    for (year in 1:N_years) {
      for (prior_year in 1:(year-1)){
        if (fire[year] == 1){
          f[year] = f[year] + spherical(year-prior_year, phi_fire);
        }
      }
      // print(fire[year]);
    }

  for (site in 1:N_sites) {
    //alpha[site][1:lag] = gamma0 + f[1:lag] * gamma2;
    //alpha[site][(lag+1):N_years] = gamma_0 + z[site][(lag+1):N_years] * gamma_1[N_sites];
    alpha[site][1:N_years] = gamma0 + f[1:N_years] * gamma2;
    //alpha[site][(lag+1):N_years] =  alpha[site][(lag+1):N_years] + f[(lag+1):N_years] * gamma2;
  }
  

} 

model { 
  
  // priors
  //tau ~ uniform(1.0/15.0, 1); 
  tau_fire ~ uniform(1.0/15.0, 1); 
  sigma_alpha ~ cauchy(1e-6, 5);
  sigma ~ cauchy(1e-6, 5);
  
  //beta ~ normal(0, 10);
  
  gamma0 ~ normal(0, 1e5);
  //gamma1 ~ normal(0, 1e5);
  gamma2 ~ normal(0, 1e5);

  // gamma0 ~ normal(0, 1e5);
  // gamma1 ~ normal(0, 1e5);
  // gamma2 ~ normal(0, 1e5);

  for (site in 1:N_sites){
    u[site] ~ normal(0, sigma_alpha);
  }

  
  // for (site in 1:N_sites){
  //   u[site] ~ cauchy(0, sigma_alpha);
  // }
  
  // eta ~ multi_normal(rep_vector(0, n_basis), sigma_eta^2 * S_inv);

  for (site in 1:N_sites){
    // working
    //Y[site]  ~ normal(alpha[site] + X[site] * beta + u[site], sigma[site]);
    //print(alpha[site][1:N_years]);
    Y[site][1:N_years]  ~ normal(alpha[site][1:N_years] + u[site][1:N_years], sigma[site]); 
    //Y[site][(lag+1):N_years]  ~ normal(alpha[site][(lag+1):N_years] + X[site][(lag+1):N_years] * beta + u[site][(lag+1):N_years], sigma); 
    //Y[site][(lag+1):N_years]  ~ normal(alpha[site][(lag+1):N_years] + X[site][(lag+1):N_years] * beta, sigma); 
    //Y[site][(lag+1):N_years]  ~ normal(alpha[site][(lag+1):N_years] + X[site][(lag+1):N_years] * beta, u[site][(lag+1):N_years]); 
  }
} 



