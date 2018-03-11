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
}

model{
    vector[N] mu;
    v ~ lognormal( 0 , 10 );
    
    for (j in 1:k) {
      b[j] ~ normal( 1 , 1 );
      a[j] ~ normal( 0 , 1 );
    }

    for ( i in 1:N ) {
      mu[i] = inv_logit(a[user_number[i]] + b[user_number[i]] * logit(expected_rating[i]));
    }  
    rating ~ beta( (mu * v) , (1 - mu) * v );  
}
