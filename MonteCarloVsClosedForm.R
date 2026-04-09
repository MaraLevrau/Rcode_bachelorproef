library(ggplot2)
library(copula)
library(microbenchmark)

trials = 10000
n = 10

# We first generate a correlation matrix:
{
  Loadings <- matrix(runif(n^2, 0, 1), nrow = n)
  Symm <- Loadings %*% t(Loadings)
  D <- diag(1 / sqrt(diag(Symm)))
  
  R <- D %*% Symm %*% D # standardize to correlation matrix
  rm(Loadings, Symm, D)
}

gaussianCopula <- normalCopula(param = R[lower.tri(R)],
                               dim = n,
                               dispstr = "un")
Unif <- rCopula(trials, gaussianCopula)

# Get some interesting parameters for the lognormal distributions
μ = runif(n, min = -1, max = 1)
σ = runif(n, min =  0, max = 1)

multivariateDist <- mvdc(
  copula = gaussianCopula,
  margins = rep("lnorm", n),
  paramMargins = lapply(1:n, function(i) list(meanlog = μ[i], sdlog = σ[i]))
)

Y <- rMvdc(trials, multivariateDist)
S <- rowSums(Y)

β = runif(n, min = 0.5, max = 1.5)

{
  # theoretical covariance between Yᵢ and Yⱼ:
  tCov <- function(i, j) {
    exp(μ[i] + μ[j] + 0.5 * (σ[i]^2 + σ[j]^2)) * (exp(R[i, j] * σ[i] * σ[j]) - 1)
  }
  
  # Theoretical correlation between Y and Λ, making use of the fact that Cov is a bilinear form:
  VarΛ = sum(sapply(1:n, function(i) {
    sum(sapply(1:n, function(j) { β[i] * β[j] * tCov(i, j) }))
  }))
  
  ρ = sapply(1:n, function(i) {
    nom = sum(sapply(1:n, function(j) { β[j] * tCov(i, j) }))
    den = sqrt(tCov(i, i) * VarΛ)
    nom / den
  })
  
  rm(tCov, VarΛ)
}

df <- data.frame()
microbenchmarkControl <- list(warmup = 10, replications = 20)

for (p in seq(0.51, 0.995, by = 0.005)) {
  # since we chose the correlation matrix to be positive definite, any choice of
  # β's with ∀ i: βᵢ > 0 will lead to a positive correlation between Y and Λ.
  
  mc = microbenchmark(quantile(S, p), control = microbenchmarkControl)
  cf = microbenchmark(sum(sapply(1:n, function (i) {
    b = μ[i] + 0.5 * (1 - ρ[i]^2) * σ[i]^2
    exp(b + σ[i] * ρ[i] * qnorm(p))
  })), control = microbenchmarkControl)
  
  VaR_monteCarlo = quantile(S, p)
  VaR_closedForm = sum(sapply(1:n, function (i) {
    b = μ[i] + 0.5 * (1 - ρ[i]^2) * σ[i]^2
    exp(b + σ[i] * ρ[i] * qnorm(p))
  }))
  
  df <- rbind(df, list(
    pLevel = p, 
    MC = VaR_monteCarlo, 
    CF = VaR_closedForm,
    
    MCTimeMin = summary(mc)$min,
    CFTimeMin = summary(cf)$min,
    MCTimeAvg = summary(mc)$mean,
    CFTimeAvg = summary(cf)$mean,
    MCTimeMax = summary(mc)$max,
    CFTimeMax = summary(cf)$max
  ))
}

ggplot(df, aes(x = pLevel)) +
  geom_line(aes(y = MC, color = "Monte Carlo")) +
  geom_line(aes(y = CF, color = "Closed Form")) +
  labs(title = "Closed form expression vs. Monte Carlo simulation", x = "Probability level p", y = "VaR_p") +
  scale_color_manual(values = c("steelblue", "red")) +
  theme_minimal() +
  theme(legend.title = element_blank())

ggplot(subset(df, pLevel > 0.6), aes(x = pLevel)) +
  geom_ribbon(aes(ymin = MCTimeMin, ymax = MCTimeMax), fill = "blue", alpha = 0.2) +
  geom_ribbon(aes(ymin = CFTimeMin, ymax = CFTimeMax), fill = "red", alpha = 0.2) +
  geom_line(aes(y = MCTimeAvg), color = "darkblue", size = 1.5) +
  geom_line(aes(y = CFTimeAvg), color = "darkred", size = 1.5) +
  scale_color_manual(values = c("blue", "red"))
  labs(title = "Computation time for Monte Carlo vs. Closed Form", x = "Probability level p", y = "Time (μs)") +
  theme_minimal()

speedupFactor <- mean(df[,"MCTimeAvg"]) / mean(df[,"CFTimeAvg"])