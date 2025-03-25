
        data {
          int<lower=0> N; int<lower=0> P;
          matrix[N, P] X;
          vector[N] y;
        }
        parameters {
          real mu;
          vector[P] beta;
          real<lower=0> lambda;
          real<lower=0> sigma_e;
        }
        model {
          mu ~ normal(0, 5);
          beta ~ double_exponential(0, lambda);
          lambda ~ gamma(1, 1);
          sigma_e ~ normal(0, 1);
          y ~ normal(mu + X * beta, sigma_e);
        }
    