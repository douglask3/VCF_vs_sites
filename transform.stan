data{
int<lower = 0> n; // number of observations
vector[n] y; // vector of observed ssb
vector[n] x; // vector of recruits
}

//transformed data{
//vector[n] logit_y; // log recruitment logit_y = log(y);
//}//

parameters{
real<lower = 0> tau1; //screwness
real<lower = 0> tau2; //screwness
real alpha; // intercept
real<lower = 0> beta; // gradient
real<lower = 0> sigma; // recruitment standard deviation 
real<lower = 0.0, upper = 1 k; // for canopy cover conversion correction
real<lower = 1, upper = 100> rou; // for canopy cover conversion correction
}


model{
real m[n];
for (i in 1:n)
  m[i] = -log(1/(fmin(k * pow(pow(pow(1-x[i], tau2)/pow(x[i], tau1), beta) +1, -rou), 1)) -1);
//  m[i] = alpha +k * log(pow(x[i], tau1)/(pow(1-x[i], tau2))) * beta;
y ~ normal(m, sigma); //account for retransformation bias
//sigma ~ cauchy(0,2.5);
//alpha ~ normal(2*max_r, 0.2*max_r);
}
