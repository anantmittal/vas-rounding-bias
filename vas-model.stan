data{
    int<lower=1> N;
    real rating[N];
    real expected_rating[N];
}

parameters{
    real a;
    real b;
    real<lower=0> v;
}

model{
    vector[N] mu;
    v ~ lognormal( 0 , 10 );
    b ~ normal( 1 , 1 );
    a ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        mu[i] = inv_logit(a + b * logit(expected_rating[i]));
    }
    rating ~ beta( (mu * v) , (1 - mu) * v );
}
