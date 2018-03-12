data{
    int<lower=1> N;
    int<lower=1> k; //Number of users
    real rating[N];
    real expected_rating[N];
    int user_number[N];
}

parameters{
    vector[k] a;
    vector[k] b;
    real<lower=0> v;
    real mu_hyper_b;
    real mu_hyper_a;
    real sigma_hyper_b;
    real sigma_hyper_a;
}

model{
    vector[N] mu;
    v ~ lognormal( 0 , 10 );
    
    mu_hyper_b ~ normal(0,1);
    sigma_hyper_b ~ cauchy(0, 1);
    
    mu_hyper_a ~ normal(0,1);
    sigma_hyper_a ~ cauchy(0, 1);
    
    for (j in 1:k) {
      b[j] ~ normal( mu_hyper_b, sigma_hyper_b);
      a[j] ~ normal( mu_hyper_a, sigma_hyper_a);
    }

    for ( i in 1:N ) {
      mu[i] = inv_logit(a[user_number[i]] + b[user_number[i]] * logit(expected_rating[i]));
    }  
    rating ~ beta( (mu * v) , (1 - mu) * v );  
}
