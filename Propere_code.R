## "Propere" code voor bachelorproef
library(mvtnorm)
library(ggplot2)
library(ggExtra)
library(sn)
library(twopiece)
library(plotly)
library(copula)
library(PerformanceAnalytics)
library(tidyr)


## Different copulas with lognormal marginal distributions


n  <- 5000
m <- 2

# Frank's copula
alpha <- 10


frankCop <- frankCopula(alpha, m)
U1 <- rCopula(n, frankCop)

X1 <- qlnorm(U1, meanlog = 0, sdlog = 1)

df1 <- as.data.frame(X1)
colnames(df1) <- c("x","y")

p1 <- ggplot(df1, aes(x, y)) + 
  geom_point() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  geom_point(size = .5) + 
  xlab("X") + ylab("Y")

p1


# Gaussian copula
correlation <- cor(U1)
rho <- correlation[1,2]

gaussianCop <- normalCopula(rho, m, dispstr = "ex")
U2 <- rCopula(n, gaussianCop)

X2 <- qlnorm(U2, meanlog = 0, sdlog = 1)


df2 <- as.data.frame(X2)
colnames(df2) <- c("x","y")


p2 <- ggplot(df2, aes(x, y)) + 
  geom_point() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  geom_point(size = .5) +
  xlab("X") + ylab("Y")

p2


## Aggregate loss with different copula's 

n  <- 10000
m <- 100

# Frank's copula
alpha    <- 10

frankCop <- frankCopula(alpha, m)
U1 <- rCopula(n, frankCop)

X1 <- qlnorm(U1, meanlog = 0, sdlog = 1)
S1 <- rowSums(X1)

ggplot(data.frame(S1), aes(S1)) +
  geom_histogram(bins = 60, fill = "darkred", color = "white") +
  labs(x = "Total Loss S",
       y = "Frequency")


var1 <- quantile(S1, 0.99)
tvar1 <- mean(S1[S1 > var1])

# Gaussian copula

corr <- cor(U1)

gaussianCop <- normalCopula(param = P2p(corr), dim = m, dispstr = "un")
U2 <- rCopula(n, gaussianCop)

X2 <- qlnorm(U2, meanlog = 0, sdlog = 1)
S2 <- rowSums(X2)

ggplot(data.frame(S2), aes(S2)) +
  geom_histogram(bins = 60, fill = "darkred", color = "white") +
  labs(x = "Total Loss S",y = "Frequency")

var2 <- quantile(S2, 0.99)
tvar2 <- mean(S2[S2 > var2])



## Value at risk and Tail value at risk

n  <- 10000
m <- 10

rho <- seq(0, 0.99, 0.05)
results <- data.frame(rho = rho, VaR = NA, TVaR = NA)

for (r in 1:length(rho)){
  corr_matrix <- matrix(rho[r], m, n)
  diag(corr_matrix) <- 1
  
  gaussianCop <- normalCopula(rho[r], m, dispstr = "ex")
  U <- rCopula(n, gaussianCop)
  X <- qlnorm(U, meanlog = 0, sdlog = 1)
  S <- rowSums(X)
  
  var <- quantile(S, 0.99)
  results$VaR[r] <- var
  results$TVaR[r] <- mean(S[S > var])
  
}

p1 <- ggplot(results, aes(x = rho, y = VaR)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(
    x = expression(Correlation),
    y = "Value at risk"
  ) +
  theme_minimal()
p1

p2 <- ggplot(results, aes(x = rho, y = TVaR)) +
  geom_line(color = "darkred", linewidth = 1) +
  labs(
    x = expression(Correlation),
    y = "Tail value at risk"
  ) +
  theme_minimal()
p2


