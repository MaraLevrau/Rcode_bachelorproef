# We want to investigate how the runtime of the `quantile` function scales with input size.

library(microbenchmark)
library(ggplot2)
library(dplyr)

times <- data.frame(
  n = integer(),
  time = numeric()
)

progress <- txtProgressBar(
  min = 1e4,
  max = 1e6,
  initial = 1e4,
  style = 3,
  char = "#")

for (n in seq(1e4, 1e6, by = 1e4)) {
  sample <- runif(n)
  benchmark <- microbenchmark(quantile(sample, 0.5), times = 1000)
  times <- rbind(times, data.frame(n = n, time = mean(benchmark$time)))
  setTxtProgressBar(progress, n)
}

ggplot(times, aes(x = n, y = time)) +
  geom_line() +
  geom_point() +
  labs(x = "Input Size (n)", y = "Average Time (nanoseconds)") +
  theme_minimal()

outliers = c()
times_without_outliers = times[seq(nrow(times)) %>% setdiff(outliers),]

# fit a linear model to the data without outliers
model <- lm(time ~ n, data = times_without_outliers)

# plot the regression line
ggplot(times_without_outliers, aes(x = n, y = time)) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = coef(model)[1], slope = coef(model)[2], color = "red") +
  labs(title = "Average execution time of the `quantile` method", x = "Input Size (n)", y = "Average `quantile` execution time over 1000 runs (nanoseconds)")

message(sprintf("quantile time scales approximately linearly as y = ax + b, where a = %f and b = %f", coef(model)[2], coef(model)[1]))
