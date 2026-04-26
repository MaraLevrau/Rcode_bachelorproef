library(ggplot2)
library(copula)
library(microbenchmark)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)
library(broman)

generateCorrelationMatrix <- function(n) {
  # We first generate a correlation matrix:
  Loadings <- matrix(runif(n^2, 0, 1), nrow = n)
  Symm <- Loadings %*% t(Loadings)
  D <- diag(1 / sqrt(diag(Symm)))
  
  D %*% Symm %*% D
}

runAll <- function(
    trials, 
    n, 
    correlationMatrix = NULL,
    βGiven = NULL,
    μRange = list(min = -1, max = 1),
    σRange = list(min = 0, max = 1),
    μGiven = NULL,
    σGiven = NULL,
    pValues = seq(0, 0.99, by = 0.01), 
    microbenchmarkControl = list(warmup = 20, replications = 50), 
    optimiseBeta = FALSE,
    inspectCorrelation = NULL) {
  message(sprintf("Running with %d trials and %d dimensions", trials, n))
  
  if (is.null(correlationMatrix)) {
    R <- generateCorrelationMatrix(n)
  } else {
    R <- correlationMatrix
  }
  
  gaussianCopula <- normalCopula(param = R[lower.tri(R)],
                                 dim = n,
                                 dispstr = "un")
  
  # Get some interesting parameters for the lognormal distributions
  if (is.null(μGiven)) {
    μ = runif(n, min = μRange$min, max = μRange$max)
  } else {
    μ = μGiven
  }
  
  if (is.null(σGiven)) {
    σ = runif(n, min = σRange$min, max = σRange$max)
  } else {
    σ = σGiven
  }
  
  multivariateDist <- mvdc(
    copula = gaussianCopula,
    margins = rep("lnorm", n),
    paramMargins = lapply(1:n, function(i) list(meanlog = μ[i], sdlog = σ[i]))
  )
  
  samplingTime <- summary(microbenchmark(
    rMvdc(trials, multivariateDist), 
    unit = "us", 
    control = microbenchmarkControl
  ))
  
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
  if (is.null(βGiven)) {
    if (optimiseBeta) {
      β <- findBestβ(rep(1, n), rep(5, n), 3, 5)
    } else {
      β <- runif(n, min = 1, max = 2)
    }
  } else {
    β <- βGiven
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
    
    quantileEstimation <- quantileSE(S, p)
    
    VaR_monteCarlo = quantileEstimation[[1]]
    VaR_closedForm = sum(sapply(1:n, function (i) {
      b = μ[i] + 0.5 * (1 - ρ[i]^2) * σ[i]^2
      exp(b + σ[i] * ρ[i] * qnorm(p))
    }))
    
    df <- rbind(df, list(
      pLevel = p, 
      MC = VaR_monteCarlo, 
      CF = VaR_closedForm,
      
      MCTimeMin = summary(mc)$min + samplingTime$min,
      CFTimeMin = summary(cf)$min,
      MCTimeAvg = summary(mc)$mean + samplingTime$mean,
      CFTimeAvg = summary(cf)$mean,
      MCTimeMax = summary(mc)$max + samplingTime$max,
      CFTimeMax = summary(cf)$max,
      
      MCError = quantileEstimation[[2]]
    ))
    
    setTxtProgressBar(progress, p)
  }
  
  message(sprintf("avg Monte Carlo time: %f μs, avg closed form time: %f μs", mean(df$MCTimeAvg), mean(df$CFTimeAvg)))
  
  rm(list = ls() %>% setdiff(c("df")))
  
  df
}

execution_times <- data.frame()
avgRuns <- 3
for (trials in seq(10e4, 10e6, by = 10e4)) {
  mcTimes <- numeric(avgRuns)
  cfTimes <- numeric(avgRuns)
  
  for (i in 1:avgRuns) {
    result <- runAll(trials = trials, n = 10, pValues = seq(0.1, 0.9, by = 0.1))
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

# split the plot into two subplots to better visualize the time ranges
pivoted <- execution_times %>%
  pivot_longer(
    cols = -trials,
    names_to = c("Method", ".value"),
    names_pattern = "(MC|CF)Time(Min|Avg|Max)"
  )

pivoted$Method <- recode(pivoted$Method,
                         MC = "Monte Carlo",
                         CF = "Closed Form")

ggplot(pivoted, aes(x = trials)) +
  geom_line(aes(y = Avg, colour = Method)) +
  geom_ribbon(aes(ymin = Min, ymax = Max, fill = Method), alpha = 0.2) +
  labs(title = "Average computation time for Closed form vs. Monte Carlo (sample size 1e6)", x = "Number of trials", y = "Average Time (μs)") +
  theme_minimal() +
  theme(legend.position = "none") + 
  facet_wrap(~Method, scales = "free_y")



# We can also look at individual calculation times in more detail

df <- runAll(
  trials = 2.5e2, 
  n = 10, 
  pValues = seq(0.01, 0.99, by = 0.005), 
  σRange = list(min = 0, max = 0.1),
  microbenchmarkControl = list(warmup = 0, replications = 0),
  optimiseBeta = FALSE)

ggplot(subset(df, pLevel > 0.8), aes(x = pLevel)) +
  geom_line(aes(y = MC, color = "Monte Carlo")) +
  geom_line(aes(y = CF, color = "Closed Form")) +
  
  geom_line(aes(y = MC - MCError, color = "Monte Carlo (lower bound)"), linetype = "dashed") +
  geom_line(aes(y = MC + MCError, color = "Monte Carlo (upper bound)"), linetype = "dashed") +
  
  geom_ribbon(aes(ymin = MC - 2 * MCError, ymax = MC + 2 * MCError, fill = "Monte Carlo"), alpha = 0.2) +
  geom_ribbon(aes(ymin = MC - MCError, ymax = MC + MCError, fill = "Monte Carlo"), alpha = 0.2) +
# geom_point(data = cfExceedsMc, aes(y = CF)) +
  labs(x = "Probability level p", y = "VaR_p") +
  scale_color_manual(values = c("steelblue", "red", "red", "red")) +
  theme(legend.title = element_blank()) + 
  theme_minimal()

stErrDeviation <- df %>%
  select(MC, CF, MCError, pLevel) %>%
  filter(MCError != 0) %>%
  mutate(StdErrDev = (CF - MC) / MCError)

maxDev = max(stErrDeviation$StdErrDev)  

ggplot(data = stErrDeviation, aes(x = pLevel)) +
  geom_line(aes(y = StdErrDev)) + 
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") + 
  geom_abline(intercept = -1, slope = 0, linetype = "dashed") +
  ylim(-maxDev, maxDev) + 
  geom_ribbon(aes(ymin=1, ymax=maxDev), fill="red", alpha=0.2) + 
  geom_ribbon(aes(ymin=-maxDev, ymax=-1), fill="red", alpha=0.2)
  

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

# by running a very large Monte Carlo simulation, we can get a very good
# estimate of the true VaR for a given set of parameters, and then we can
# compare the closed form estimate to that to see how well it performs for 
# different ranges of σ
σRangePlots <- list()
n = 20
μValues = runif(n, min = -2, max = 2)
pValues = seq(0.1, 0.9, by = 0.05)

rows = 3
cols = 4
σMin = 0.1
σMax = 5
diff = (σMax - σMin) / (rows * cols)

for (σRangeMax in seq(σMin, σMax - diff, by = diff)) {
  message(sprintf("Running for σ range max = %f", σRangeMax))
  σValues = runif(n, min = σRangeMax, max = σRangeMax + diff)
  
  accurateMCResult <- runAll(
    trials = 1e6, 
    n = n, 
    pValues = pValues, 
    microbenchmarkControl = list(warmup = 0, replications = 0),
    μGiven = μValues,
    σGiven = σValues
  )
  
  accurateMCResult <- select(accurateMCResult, pLevel, MC)
  
  CFResult <- runAll(
    trials = 1e4,
    n = n,
    pValues = pValues,
    microbenchmarkControl = list(warmup = 0, replications = 0),
    μGiven = μValues,
    σGiven = σValues
  )
  
  CFResult <- select(CFResult, pLevel, CF)
  
  comparison <- merge(accurateMCResult, CFResult, by = "pLevel") %>%
    mutate(σRangeMax = σRangeMax)
  
  plot <- ggplot(data = comparison, aes(x = pLevel)) + 
    geom_line(aes(y = MC, color = "Monte Carlo"), linewidth = 1) +
    geom_line(aes(y = CF, color = "Closed form"), linewidth = 1) + 
    labs(title = sprintf("σ ∈ [%.2f, %.2f]", σRangeMax, σRangeMax + diff), y = NULL, x = NULL)
  
  σRangePlots[[length(σRangePlots) + 1]] <- plot
}

wrap_plots(σRangePlots, ncol = 4, guides = "collect") + 
  plot_annotation(title = "Comparison of Monte Carlo and Closed Form VaR estimates for different σ ranges") &
  theme(legend.position = "bottom", legend.title = element_blank())

## Let's inspect the way β's affect the predictions

βComparisonDf <- data.frame()

R <- generateCorrelationMatrix(3)

μ = runif(3, min = -1, max = 1)
σ = runif(3, min =  0, max = 1)

pValues <- seq(0.01, 0.99, by = 0.001)
betas <- tibble(names= c("β1", "β2", "β3"), comp1 = c(5, 1, 1), comp2 = c(1, 5, 1), comp3 = c(1, 1, 5))

gaussianCopula <- normalCopula(param = R[lower.tri(R)], dim = 3, dispstr = "un")

multivariateDist <- mvdc(
  copula = gaussianCopula,
  margins = rep("lnorm", 3),
  paramMargins = lapply(1:3, function(i) list(meanlog = μ[i], sdlog = σ[i]))
)

Y <- rMvdc(100, multivariateDist)

for (i in 1:3) {
  β = c(betas$comp1[i], betas$comp2[i], betas$comp3[i])
  
  {
    # theoretical covariance between Yᵢ and Yⱼ:
    tCov <- function(i, j) {
      exp(μ[i] + μ[j] + 0.5 * (σ[i]^2 + σ[j]^2)) * (exp(R[i, j] * σ[i] * σ[j]) - 1)
    }
    
    # Theoretical correlation between Y and Λ, making use of the fact that Cov is a bilinear form:
    VarΛ = sum(sapply(1:3, function(i) {
      sum(sapply(1:3, function(j) { β[i] * β[j] * tCov(i, j) }))
    }))
    
    ρ = sapply(1:3, function(i) {
      nom = sum(sapply(1:3, function(j) { β[j] * tCov(i, j) }))
      den = sqrt(tCov(i, i) * VarΛ)
      nom / den
    })
    
    rm(tCov, VarΛ)
  }
  
  # covariance size
  euclNorm <- sqrt(sum(ρ^2))
  message(sprintf("beta %s : cov %f", betas$names[i], euclNorm))
  result <- data.frame()
  for (p in pValues) {
    VaR_closedForm = sum(sapply(1:3, function (i) {
      b = μ[i] + 0.5 * (1 - ρ[i]^2) * σ[i]^2
      exp(b + σ[i] * ρ[i] * qnorm(p))
    }))
    
    result <- rbind(result, list(
      pLevel = p, 
      CF = VaR_closedForm
    ))
    
  }
  
  βComparisonDf <- rbind(βComparisonDf, result)
}

MCValues <- runAll(
  trials = 1e6, 
  n = 3, 
  pValues = pValues, 
  microbenchmarkControl = list(warmup = 0, replications = 0),
  μGiven = μ,
  σGiven = σ,
  βGiven = β,
  correlationMatrix = R
) %>% select(pLevel, MC)

βComparisonDfMerged <- βComparisonDf %>%
  pivot_wider(names_from = βName, values_from = c(CF)) %>%
  merge(MCValues, by = "pLevel") %>%
  mutate(Δβ1 = MC - β1, Δβ2 = MC - β2, Δβ3 = MC - β3) %>%
  select(pLevel, MC, Δβ1, Δβ2, Δβ3)

sum(abs(βComparisonDfMerged$Δβ1))
sum(abs(βComparisonDfMerged$Δβ2))
sum(abs(βComparisonDfMerged$Δβ3))

ggplot(βComparisonDfMerged, aes(x = pLevel)) +
  geom_line(aes(y = Δβ1, color = "β₁"), linewidth = 1) + 
  geom_line(aes(y = Δβ2, color = "β₂"), linewidth = 1) + 
  geom_line(aes(y = Δβ3, color = "β₃"), linewidth = 1) + 
  geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed", color = "black") +
  labs(title = "Comparison of Monte Carlo and Closed Form VaR estimates for different β's", y = "VaR", x = "Probability level p") +
  theme_minimal() +
  theme(legend.title = element_blank())
