library(ggplot2)
library(copula)
library(microbenchmark)
library(gridExtra)
library(tidyr)
library(dplyr)

runAll <- function(
    trials, 
    n, 
    pValues = seq(0, 0.99, by = 0.01), 
    microbenchmarkControl = list(warmup = 20, replications = 50), 
    optimiseBeta = FALSE,
    inspectCorrelation = NULL) {
  message(sprintf("Running with %d trials and %d dimensions", trials, n))
  
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
  
  # since we chose the correlation matrix to be positive definite, any choice of
  # β's with ∀ i: βᵢ > 0 will lead to a positive correlation between Y and Λ.
  
  if (optimiseBeta) {
    β <- findBestβ(rep(1, n), rep(5, n), 3, 5)
  } else {
    β = runif(n, min = 1, max = 2)
  }
  
  # The following is debug code to draw a heatmap for correlations associated with β's 
  # to check if there is some global maximum to aim for
  # if this is not the case, the whole endeavour is futile
  # because of the huge expand.grid, we will only do this for a small number of 
  # coordinates and a coarse grid, but it should be enough to get an idea of the landscape
  if (!is.null(inspectCorrelation)) {
    linspace <- seq(from = inspectCorrelation$from, to = inspectCorrelation$to, by = inspectCorrelation$by)
    searchSpace <- expand.grid(lapply(1:n, function(x) linspace))
    
    gridDF <- data.frame()
    
    progress <- txtProgressBar(max=(length(linspace)^n), initial=0, style=3, char="#")
    
    for (i in 1:(length(linspace)^n)) {
      β <- c(searchSpace[i, 1:n], use.names = FALSE, recursive = TRUE)
      # message(paste(β, collapse = ", "))
      correlation <- cor(Y, Y %*% β)
      euclNorm <- (t(correlation) %*% correlation)[1,1]
      
      coordValues <- as.list(β)
      names(coordValues) <- lapply(1:n, function (num) paste("coord_", num, sep=""))
      
      value <- c(val = euclNorm, coordValues)
      
      # store the position and value for the heatmap
      gridDF <- rbind(gridDF, value)
      
      setTxtProgressBar(progress, i)
    }
    
    # view pairwise plots as heatmaps per coordinate pair in a grid
    plotList <- list()
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        plotList[[length(plotList) + 1]] <- ggplot(gridDF, aes_string(x = paste("coord_", i, sep=""), y = paste("coord_", j, sep=""))) +
          geom_tile(aes(fill = val)) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = quantile(gridDF[,"val"], 0.75)) +
          labs(title = sprintf("coord_%d vs. coord_%d", i, j), x = sprintf("coord_%d", i), y = sprintf("coord_%d", j)) +
          theme_minimal() + 
          guides(fill="none")
      }
    }
    
    
    do.call("grid.arrange", c(plotList, ncol = n-1))
    return(gridDF)
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
  
  if (length(pValues) == 1) {
    progress <- txtProgressBar(
      min = 0,
      max = max(pValues),
      initial = 0,
      style = 3,
      char = "#")
  } else {
    progress <- txtProgressBar(
      min = min(pValues), 
      max = max(pValues), 
      initial = min(pValues), 
      style = 3, 
      char = "#")
  }
  
  for (p in pValues) {
    mc = microbenchmark(quantile(S, p), unit = "us", control = microbenchmarkControl)
    cf = microbenchmark(sum(sapply(1:n, function (i) {
      b = μ[i] + 0.5 * (1 - ρ[i]^2) * σ[i]^2
      exp(b + σ[i] * ρ[i] * qnorm(p))
    })), unit = "us", control = microbenchmarkControl)
    
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
    
    setTxtProgressBar(progress, p)
  }
  
  message(sprintf("avg Monte Carlo time: %f μs, avg closed form time: %f μs", mean(df$MCTimeAvg), mean(df$CFTimeAvg)))
  
  rm(list = ls() %>% setdiff(c("df")))
  
  df
}

execution_times <- data.frame()
avgRuns <- 1
for (trials in seq(10e4, 10e6, by = 10e4)) {
  mcTimes <- numeric(avgRuns)
  cfTimes <- numeric(avgRuns)
  
  for (i in 1:avgRuns) {
    result <- runAll(trials = trials, n = 5, pValues = seq(0.1, 0.9, by = 0.1))
    mcTimes[i] <- mean(result$MCTimeAvg)
    cfTimes[i] <- mean(result$CFTimeAvg)
  }
  
  execution_times <- rbind(execution_times, list(
    trials = trials, 
    MCTimeMin = min(mcTimes),
    MCTimeAvg = mean(mcTimes),
    MCTimeMax = max(mcTimes),
    CFTimeMin = min(cfTimes),
    CFTimeAvg = mean(cfTimes),
    CFTimeMax = max(cfTimes)
  ))
}

speedup <- mutate(execution_times, Speedup = MCTimeAvg / CFTimeAvg)

ggplot(speedup, aes(x = trials)) +
  geom_line(aes(y = Speedup), color = "blue") +
  labs(title = "Speedup of Closed Form over Monte Carlo", x = "Number of trials", y = "Speedup (Closed Form Time / Monte Carlo Time)") +
  theme_minimal()

ggplot(execution_times, aes(x = trials)) +
  geom_line(aes(y = MCTimeAvg, colour = "Monte Carlo")) +
  geom_line(aes(y = CFTimeAvg, colour = "Closed Form")) +
  geom_ribbon(aes(ymin = MCTimeMin, ymax = MCTimeMax, fill = "Monte Carlo"), alpha = 0.2) +
  geom_ribbon(aes(ymin = CFTimeMin, ymax = CFTimeMax, fill = "Closed Form"), alpha = 0.2) +
  labs(title = "Average computation time for Monte Carlo vs. Closed Form (sample size 1e6)", x = "Number of trials", y = "Average Time (μs)") +
  theme_minimal() +
  theme(legend.title = element_blank())



# We can also look at individual calculation times in more detail

df <- runAll(
  trials = 1e4, 
  n = 10, 
  pValues = seq(0, 0.99, by = 0.01), 
  microbenchmarkControl = list(warmup = 20, replications = 100),
  optimiseBeta = FALSE)

ggplot(df, aes(x = pLevel)) +
  geom_line(aes(y = MC, color = "Monte Carlo")) +
  geom_line(aes(y = CF, color = "Closed Form")) +
# geom_point(data = cfExceedsMc, aes(y = CF)) +
  labs(x = "Probability level p", y = "VaR_p") +
  scale_color_manual(values = c("steelblue", "red")) +
  theme(legend.title = element_blank()) + 
  theme_minimal()

differences <- df %>%
  mutate(Difference = MC - CF) %>%
  mutate(DifferenceWhenCFExceeds = ifelse(Difference < 0, Difference, NA))

differences <- subset(differences, pLevel != 0)

ggplot(differences, aes(x = pLevel)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", color = "black") +
  geom_area(aes(y = DifferenceWhenCFExceeds), fill = "red", alpha = 0.2) + 
  geom_line(aes(y = Difference), color = "blue") +
  labs(title = "Difference between Closed Form and Monte Carlo VaR estimates", x = "Probability level p", y = "Difference (Closed Form - Monte Carlo)")

ggplot(subset(df, MCTimeMax < 10000), aes(x = pLevel)) +
  geom_ribbon(aes(ymin = MCTimeMin, ymax = MCTimeMax, fill = "Monte Carlo"), alpha = 0.2) +
  geom_ribbon(aes(ymin = CFTimeMin, ymax = CFTimeMax, fill = "Closed form"), alpha = 0.2) +
  geom_line(aes(y = MCTimeAvg, colour = "Monte Carlo"), linewidth = 1) +
  geom_line(aes(y = CFTimeAvg, colour = "Closed form"), linewidth = 1) +
  labs(title = "Computation time for Monte Carlo vs. Closed form (sample size 1e4)", x = "Probability level p", y = "Time (μs)") +
  guides(fill = guide_legend(title = "Method"))


merged_times <- df %>%
  select(pLevel, MCTimeAvg, CFTimeAvg) %>%
  pivot_longer(cols = c(MCTimeAvg, CFTimeAvg), names_to = "Method", values_to = "Time")

ggplot(subset(merged_times, Time < 400), aes(x = Time, fill=Method)) +
  geom_histogram(binwidth=4, alpha=0.9, position="identity") +
  theme_minimal()

# inspect correlation
gridDF <- runAll(
  trials = 50, 
  n = 3, 
  pValues = seq(0, 0.99, by = 0.01), 
  microbenchmarkControl = list(warmup = 20, replications = 100),
  inspectCorrelation = list(from = 0.5, to = 10, by = 0.15)
)
