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
library(dplyr)


## Different copulas with lognormal marginal distributions


n  <- 5000
m <- 2

# Frank's copula
alpha <- 55
frankCop <- frankCopula(alpha, m)
U1 <- rCopula(n, frankCop)
X1 <- qnorm(U1, mean = 0, sd = 1)
df1 <- as.data.frame(X1)
colnames(df1) <- c("x","y")

# Gaussian copula
correlation <- cor(U1)
rho <- correlation[1,2]

gaussianCop <- normalCopula(rho, m, dispstr = "ex")
U2 <- rCopula(n, gaussianCop)
X2 <- qnorm(U2, mean = 0, sd = 1)
df2 <- as.data.frame(X2)
colnames(df2) <- c("x","y")

# we voegen df1 en df2 samen
df1$Copula <- "Frank's"
df2$Copula <- "Gaussian"
df <- rbind(df1, df2)
ggplot(df, aes(x, y, color = Copula)) + 
  geom_point(size = .5) + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("X") + ylab("Y") +
  scale_color_manual(values = c("darkred", "steelblue")) +
  labs(color = "Copula type") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))


## Aggregate loss with different copula's 

n  <- 10000
m <- 100

rho <- 0.7

gaussianCop <- normalCopula(rho, m, dispstr = "ex")
U2 <- rCopula(n, gaussianCop)

X2 <- qlnorm(U2, meanlog = 0, sdlog = 1)
S2 <- rowSums(X2)

ggplot(data.frame(S2), aes(S2)) +
  geom_histogram(bins = 60, fill = "darkred", color = "white") +
  labs(x = "Total Loss S",y = "Frequency")

var2 <- quantile(S2, 0.99)
tvar2 <- mean(S2[S2 > var2])


## Value at risk and Tail value at risk

n  <- 25000
m <- 10
repetitions <- 50

rho <- seq(0, 0.99, 0.05)
results <- data.frame()

for (r in 1:length(rho)){
  for (repetition in 1:repetitions) {
    gaussianCop <- normalCopula(rho[r], m, dispstr = "ex")
    U <- rCopula(n, gaussianCop)
    X <- qlnorm(U, meanlog = 0, sdlog = 1)
    S <- rowSums(X)
    
    var <- quantile(S, 0.99)
    tvar <- mean(S[S > var])
    
    results <- rbind(results, data.frame(rho = rho[r], repetition = repetition, VaR = var, TVaR = tvar))
  }
}

summary_df <- results %>%
  group_by(rho) %>%
  summarise(min_VaR = min(VaR), max_VaR = max(VaR), mean_VaR = mean(VaR),
            min_TVaR = min(TVaR), max_TVaR = max(TVaR), mean_TVaR = mean(TVaR))

p1 <- ggplot() +
  geom_line(data=results, aes(x=rho, y=VaR, group=repetition), alpha=0.3, color="gray") +
  geom_ribbon(data=summary_df, aes(x=rho, ymin=min_VaR, ymax=max_VaR), alpha=0.3, fill="lightblue") +
  geom_line(data=summary_df, aes(x=rho, y=mean_VaR), color="black", linewidth=1.5) +
  labs(
    x = expression(Correlation),
    y = "Value at risk (q = 0.99)"
  ) +
  theme_minimal(base_size = 16)
p1

p2 <- ggplot() +
  geom_line(data=results, aes(x=rho, y=TVaR, group=repetition), alpha=0.3, color="gray") +
  geom_ribbon(data=summary_df, aes(x=rho, ymin=min_TVaR, ymax=max_TVaR), alpha=0.3, fill="lightcoral") +
  geom_line(data=summary_df, aes(x=rho, y=mean_TVaR), color="black", linewidth=1.5) +
  labs(
    x = expression(Correlation),
    y = "Tail value at risk (q = 0.99)"
  ) +
  theme_minimal(base_size = 16)
p2
