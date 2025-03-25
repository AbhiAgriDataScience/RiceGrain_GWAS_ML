
        data {
          int<lower=0> N; int<lower=0> P;
          matrix[N, P] X;
          vector[N] y;
        }
        parameters {
          real mu;
          vector[P] beta;
          vector<lower=0>[P] tau_sq;
          real<lower=0> sigma_e;
        }
        model {
          mu ~ normal(0, 5);
          tau_sq ~ inv_gamma(1, 1);
          for (j in 1:P)
            beta[j] ~ normal(0, sqrt(tau_sq[j]));
          sigma_e ~ normal(0, 1);
          y ~ normal(mu + X * beta, sigma_e);
        }
    