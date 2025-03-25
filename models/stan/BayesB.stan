
        data {
          int<lower=0> N; int<lower=0> P;
          matrix[N, P] X;
          vector[N] y;
        }
        parameters {
          real mu;
          vector[P] beta;
          real<lower=0> sigma_e;
          real<lower=0, upper=1> pi;
        }
        model {
          mu ~ normal(0, 5);
          pi ~ beta(1, 1);
          for (j in 1:P)
            target += log_mix(pi,
                              normal_lpdf(beta[j] | 0, 1),
                              normal_lpdf(beta[j] | 0, 0.01));
          sigma_e ~ normal(0, 1);
          y ~ normal(mu + X * beta, sigma_e);
        }
    