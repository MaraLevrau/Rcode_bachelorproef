library(ggplot2)
library(copula)
library(microbenchmark)
library(gridExtra)
library(tidyr)
library(dplyr)

trials = 1e4
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

df <- data.frame()
microbenchmarkControl <- list(warmup = 20, replications = 5)


# Find some vector β such that Cor(Y, β · Y) is high; for simplicity we will do
# this using the builtin cor() function since that simplifies the code a bit

findBestβ <- function(minV, maxV, sieve, iterationsLeft) {
  stepSizes <- (maxV - minV) / sieve
  linspaces <- lapply(1:n, function(i) seq(from = minV[i], to = maxV[i], by = stepSizes[i]))
  searchSpace <- expand.grid(linspaces, KEEP.OUT.ATTRS = FALSE)
  
  bestβ = rep(NA, n)
  bestCor = -Inf
  
  for (i in 1:sieve^n) {
    β <- c(searchSpace[i, 1:n], use.names = FALSE, recursive = TRUE)
    correlation <- cor(Y, Y %*% β)
    euclNorm <- (t(correlation) %*% correlation)[1,1]
    if (euclNorm > bestCor) {
      bestβ <- β
      bestCor = euclNorm
    }
  }
  
  message(sprintf("best β with %d iterations left is %s", iterationsLeft, paste(bestβ, collapse = ", ")))
  
  if (iterationsLeft == 0) {
    return(bestβ)
  } else {
    return(findBestβ(bestβ - stepSizes / 2, bestβ + stepSizes / 2, sieve, iterationsLeft - 1))
  }
}

#  β = findBestβ(rep(1, n), rep(5, n), 2, 5)
β = runif(n, min = 1, max = 2)

# The following is debug code to draw a heatmap for correlations associated with β's 
# to check if there is some global maximum to aim for
# if this is not the case, the whole endeavour is futile
# because of the huge expand.grid, we will only do this for a small number of 
# coordinates and a coarse grid, but it should be enough to get an idea of the landscape
if (FALSE) {
  linspace <- seq(from = 0.05, to = 3, by = 0.1)
  searchSpace <- expand.grid(lapply(1:n, function(x) linspace))
  
  grid <- data.frame()
  
  for (i in 1:(length(linspace)^n)) {
    β <- c(searchSpace[i, 1:n], use.names = FALSE, recursive = TRUE)
    message(paste(β, collapse = ", "))
    correlation <- cor(Y, Y %*% β)
    euclNorm <- (t(correlation) %*% correlation)[1,1]
    
    coordValues <- as.list(β)
    names(coordValues) <- lapply(1:n, function (num) paste("coord_", num, sep=""))
    
    value <- c(val = euclNorm, coordValues)
    
    # store the position and value for the heatmap
    grid <- rbind(grid, value)
  }
  
  # view pairwise plots as heatmaps per coordinate pair in a grid
  plotList <- list()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      plotList[[length(plotList) + 1]] <- ggplot(grid, aes_string(x = paste("coord_", i, sep=""), y = paste("coord_", j, sep=""))) +
        geom_tile(aes(fill = val)) +
        scale_fill_gradient(low = "white", high = "red", oob = scales::squish, limits = c(quantile(grid[,"val"], 0.75), max(grid[,"val"]))) +
        labs(title = sprintf("coord_%d vs. coord_%d", i, j), x = sprintf("coord_%d", i), y = sprintf("coord_%d", j)) +
        theme_minimal() + 
        guides(fill="none")
    }
  }
  
  
  do.call("grid.arrange", c(plotList, ncol = n-1))
}

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

for (p in seq(0, 0.99, by = 0.01)) {
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
  
  message(sprintf("ran simulation for p = %f", p))
}

cfExceedsMc <- subset(df, CF > MC)
speedupFactor <- mean(df[,"MCTimeAvg"]) / mean(df[,"CFTimeAvg"])

ggplot(df, aes(x = pLevel)) +
  geom_line(aes(y = MC, color = "Monte Carlo")) +
  geom_line(aes(y = CF, color = "Closed Form")) +
  geom_point(data = cfExceedsMc, aes(y = CF)) + 
  labs(title = "Closed form expression vs. Monte Carlo simulation", x = "Probability level p", y = "VaR_p") +
  scale_color_manual(values = c("steelblue", "red")) +
  theme(legend.title = element_blank())

ggplot(subset(df, pLevel > 0.1), aes(x = pLevel)) +
  geom_ribbon(aes(ymin = MCTimeMin, ymax = MCTimeMax), fill = "blue", alpha = 0.2) +
  geom_ribbon(aes(ymin = CFTimeMin, ymax = CFTimeMax), fill = "red", alpha = 0.2) +
  geom_line(aes(y = MCTimeAvg), color = "darkblue", linewidth = 1.5) +
  geom_line(aes(y = CFTimeAvg), color = "darkred", linewidth = 1.5) +
  scale_color_manual(values = c("blue", "red")) + 
  labs(title = "Computation time for Monte Carlo vs. Closed Form", x = "Probability level p", y = "Time (μs)") +
  theme_minimal()


merged_times <- df %>%
  select(pLevel, MCTimeAvg, CFTimeAvg) %>%
  pivot_longer(cols = c(MCTimeAvg, CFTimeAvg), names_to = "Method", values_to = "Time")

ggplot(merged_times, aes(x = Time, fill=Method)) + 
  geom_histogram(binwidth=4, alpha=0.9, position="identity") + 
  theme_minimal()
