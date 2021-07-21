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
  int<lower=0> N_covars;
  int<lower=0> n_basis;
  int<lower=0> n_knots;
  
  vector[N_years] Y[N_sites];
  matrix[N_years, N_covars] X[N_sites];
  vector[N_years] cmem[N_sites];
  vector[N_years] dmem;

  matrix[n_knots, n_basis] B;
  matrix[n_basis, n_basis] S_inv;

  int X_nmiss;
  int d_nmiss;

  int X_index[X_nmiss];
  int d_index[d_nmiss];
  
} 
transformed data {
    real sigma_eta;
    vector[N_years] logY[N_sites];

    sigma_eta = 0.1;

    for (site in 1:N_sites){
      logY[site] = log(Y[site]);
    }

   
}

parameters {  
  real<lower=1.0/15.0, upper=1> tau_dmem; 
  vector<lower=1e-6>[N_sites] sigma;
  real<lower=1e-6> sigma_alpha;
  vector[N_years-lag] u[N_sites];
  
  vector[N_covars] beta;

  real gamma0;
  real gamma1;
  vector[N_sites] gamma2;
  
  vector[n_basis] eta;

  // MVN covariates; for imputation
  cholesky_factor_cov[N_covars+2] Sigma;
  vector[N_covars+2] mu[N_sites];

  // Missings
  //matrix[X_nmiss, N_covars] X_imp[N_sites];
  vector<lower=0>[X_nmiss] X1_imp[N_sites];
  vector[X_nmiss] X2_imp[N_sites];
  vector[d_nmiss] d_imp[N_sites];
} 

transformed parameters { 
  //real phi;
  real phi_dmem;
  //matrix[N_years, N_sites] z;
  //matrix[N_years, N_sites] alpha;
  //vector[N_years-lag] z[N_sites];
  vector[N_years] z[N_sites];
  vector[N_years] f;
  vector[N_years-lag] alpha[N_sites];

  matrix[N_years, N_covars+2] dat[N_sites];
  

  vector[lag+1] exp_B_eta;
  real sum_exp_B_eta;

  //vector[lag+1] w_incr;
  vector[lag+1] w;
  //real w_sum;

  //vector[Z_nmiss] logX2_imp[N_sites];

  //logX2_imp = log(X2_imp);
  
  for (site in 1:N_sites){
    dat[site][, 1] = logY[site];
    //dat[site][, 1] = Y[site];
    dat[site][, 2:(N_covars+1)] = X[site];
    dat[site][, N_covars+2] = cmem[site];

    //dat[site][X_index, 2:(N_covars+1)] = X_imp[site];
    dat[site][X_index, 2] = X1_imp[site];
    //dat[site][X_index, 3] = logX2_imp[site];
    dat[site][X_index, 3] = X2_imp[site];
    dat[site][d_index, N_covars+2] = d_imp[site];
  }
  


  
  exp_B_eta = exp(B * eta);
  sum_exp_B_eta = sum(exp_B_eta);

  w = exp_B_eta / sum_exp_B_eta;
  //for (t in 1:(lag+1)){ 
    //w_incr[t]  = w[lag+1-(t-1)];
  //}
 
  
  
  //w = exp(w);
  //w_sum = sum(w);
  //w = w/w_sum;
  //w = reverse(w_incr);
  
  //  for (site in 1:N_sites) {
  //    z[site] = rep_vector(0, N_years-lag);
  //    for (year in 1:(N_years-lag)) {
  //     for (l in 0:lag){
  //	//z[site][year] = z[site][year] + d[site][year-l] * w[l+1];
  //	z[site][year] = z[site][year] + dat[site][year+lag-l,N_covars+2] * w[l+1]; 
  //     }
  // }
  //}


  for (site in 1:N_sites) {
    z[site] = rep_vector(0, N_years);
    for (year in (lag+1):N_years) {
      for (l in 0:lag){
	//z[site][year] = z[site][year] + d[site][year-l] * w[l+1];
	z[site][year] = z[site][year] + dat[site][year-l, N_covars+2] * w[l+1];
	//z[site][year] = z[site][year] + dat[site][year+lag-l,N_covars+2] * w[l+1]; 
      }
    }
  }
  
  
  phi_dmem = 1/tau_dmem;
  

    f = rep_vector(0, N_years);
    for (year in 1:N_years) {
      for (prior_year in 1:(year-1)){
        if (dmem[year] == 1){
          f[year] = f[year] + spherical(year-prior_year, phi_dmem);
        }
      }
      //print(f[year]);
    }
  
  for (site in 1:N_sites) {
    alpha[site] = gamma0 + z[site][(lag+1):N_years] * gamma1 + f[(lag+1):N_years] * gamma2[site];

    //working
    //alpha[site][(lag+1):N_years] = gamma0 + z[site][(lag+1):N_years] * gamma1 + f[(lag+1):N_years] * gamma2;

  }
  

} 

model { 
  
  // priors
  //tau ~ uniform(1.0/15.0, 1); 
  tau_dmem ~ uniform(1.0/15.0, 1); 
  sigma_alpha ~ cauchy(1e-6, 5);
  sigma ~ cauchy(1e-6, 5);
  
  beta ~ normal(0, 10);
  
  gamma0 ~ normal(0, 2);
  gamma1 ~ normal(0, 1e5);
  gamma2 ~ normal(0, 1e5);

  // gamma0 ~ normal(0, 1e5);
  // gamma1 ~ normal(0, 1e5);
  // gamma2 ~ normal(0, 1e5);

  for (site in 1:N_sites){
    mu[site] ~ normal(0,2);
  }

  // MVN imputation
  for (site in 1:N_sites){
    for(year in 1:N_years){
      dat[site][year,] ~ multi_normal_cholesky(mu[site],Sigma);
    }
  }

  for (site in 1:N_sites){
    u[site] ~ normal(0, sigma_alpha);
  }
  
  eta ~ multi_normal(rep_vector(0, n_basis), sigma_eta^2 * S_inv);

  for (site in 1:N_sites){

    logY[site][(lag+1):N_years]  ~ normal(alpha[site] + dat[site][(lag+1):N_years,2:(N_covars+1)] * beta + u[site], sigma[site]);
    
    // dat[site][1:N_years,1]  ~ normal(alpha[site][1:N_years] + dat[site][1:N_years,2:(N_covars+1)] * beta + u[site][1:N_years], sigma[site]); 
    //Y[site][1:N_years]  ~ normal(alpha[site][1:N_years] + dat[site][1:N_years,2:(N_covars+1)] * beta + u[site][1:N_years], sigma[site]); 

    // working
    //Y[site][(lag+1):N_years]  ~ normal(alpha[site][(lag+1):N_years] + X[site][(lag+1):N_years] * beta + u[site][(lag+1):N_years], sigma[site]); 
    
  }
} 



