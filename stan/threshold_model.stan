data {
  int<lower=0> N;          // Number of clusters
  int<lower=0> J;          // Total number of observations
  int<lower=1,upper=N> cluster_idx[J];  // Cluster index for each observation
  real y[J];               // Observations
}
parameters {
  real mu[N];              // Mean expression for each cluster
  real<lower=0> sigma[N];  // Cluster-level standard deviations
  real T;                  // Global threshold for marking
  real<lower=0> alpha;     // Slope controlling sharpness of decision boundary
}
model {
  // Priors
  mu ~ normal(0.3699511, 0.4939401);      // Prior for cluster means
  sigma ~ cauchy(0.07308471, 0.1083554);    // Prior for cluster standard deviations
  T ~ normal(1, 0.5);      // Prior for threshold
  alpha ~ cauchy(0, 1);    // Prior for decision boundary slope

  // Likelihood
  for (j in 1:J) {
    y[j] ~ normal(mu[cluster_idx[j]], sigma[cluster_idx[j]]);
  }
}
generated quantities {
  real marked_prob[N];     // Probability of being marked for each cluster
  for (i in 1:N) {
    marked_prob[i] = inv_logit(alpha * (mu[i] - T));
  }
}
