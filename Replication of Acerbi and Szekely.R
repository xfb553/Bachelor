Replikation af artiklen "Backtesting Expected Shortfall" fra Acerbi & Szekely, 2014

Pakker
```{r}
library(esreg)
library(stats)
library(ggplot2)
library(dplyr)
```

VaR og ES
```{r}
VaR <- function(alpha = 0.975, sigma2 = 1, mu = 0, shift = 0, scale = 1, df = NULL, type = "Normal", Norm = FALSE) {
  if (!Norm) {
    if (type == "Normal") {
      VaR <- mu + sqrt(sigma2) * qnorm(1 - alpha)
    } else if (type == "t") {
      VaR <- qt(1 - alpha, df = df) * scale+shift
    }
  } else {
    if (type == "Normal") {
      VaR <- qnorm(1 - alpha)
    } else if (type == "t") {
      VaR <- sqrt((df - 2) / df) * qt(1 - alpha, df = df)*scale+shift
    }
  }
  return(-VaR)
}

ES <- function(alpha = 0.975, sigma2 = 1, mu = 0, shift = 0, scale = 1, df = NULL, type = "Normal", Norm = FALSE) {
  if (!Norm) {
    if (type == "Normal") {
      ES <- (mu - sqrt(sigma2) / (1 - alpha) * dnorm(qnorm(1 - alpha))) * scale + shift
    } else if (type == "t") {
      x <- qt(1 - alpha, df = df)
      ES <- -(dt(x, df) / (1 - alpha) * (df + x^2) / (df - 1)) * scale + shift
    }
  } else {
    if (type == "Normal") {
      ES <- (mu - sqrt(sigma2) / (1 - alpha) * dnorm(qnorm(1 - alpha))) * scale + shift
    } else if (type == "t") {
      x <- qt(1 - alpha, df = df)
      ES <- -sqrt((df - 2) / df) * (dt(x, df) / (1 - alpha) * (df + x^2) / (df - 1))*scale+shift

    }
  }
  return(-ES)
}
```

Shift og scale
```{r}
shift <- function(alpha, df, dfnull, Norm=FALSE) {
  VaR1 <- VaR(1-alpha, df = dfnull, Norm = Norm, type = "t")
  VaR2 <- VaR(1-alpha, df = df, Norm = Norm, type = "t")
  res <- (VaR1 - VaR2)
  return(res)
}

scale <- function(signi = 0.05, df = 100, type = "Normal", Norm=FALSE) {
  res <- ES(alpha = 0.975, sigma2 = 1, mu = 0, df = df, type = type, Norm=Norm) / 
    ES(alpha = 1 - signi, sigma2 = 1, mu = 0, df = df, type = type, Norm=Norm)
  return(res)
}
```

Test 1
```{r}

Z_1 <- function(quan, sigma2, mu, T, n,scale=1, shift = 0, df = NULL, dfnull = 100, Norm=TRUE, type = "Normal") {
  Z_1_list <- c()
  var <- VaR(alpha = 0.975, sigma2 = sigma2, mu = mu, df = dfnull, type = type)
  es <- ES(alpha = 0.975, sigma2 = sigma2, mu = mu, df = dfnull, type = type)

  for (i in 1:n) {
    if(!Norm){
      X_t <- rt(T, df = df)+shift
    } else {
      X_t<-rt(T,df=df)*sqrt((df-2)/df)+shift
    }
    q <- 0
    P <- var + X_t
    N_t <- sum(P < 0)

    for (j in 1:T) {
      if (X_t[j] + var < 0) {
        q <- q + X_t[j] / es
      }
    }
    Z_1 <- q / N_t + 1
    Z_1_list <- c(Z_1_list, Z_1)
  }

  Z_1_list <- Z_1_list[!is.na(Z_1_list)]
  quant <- quantile(Z_1_list, quan)

  return(list(quant = quant, Z_1_list = Z_1_list, mean = mean(Z_1_list)))
}
```

Test 2
```{r}
Z_2 <- function(quan, sigma2, mu, T, n, shift = 0, scale = 1, df = NULL, dfnull = 100, Norm = FALSE, type = "Normal") {
  Z_2_list <- c()
  var <- VaR(alpha = 0.975, sigma2 = sigma2,mu = mu, df = dfnull, type = type, Norm = Norm)
  es <- ES(alpha = 0.975, sigma2 = sigma2,mu = mu, df = dfnull, type = type, Norm = Norm)

  for (i in 1:n) {
    if (!Norm) {
      X_t <- rt(T, df = df) * scale+shift
    } else {
      X_t <- rt(T, df = df) * sqrt((df - 2) / df)+shift
    }

    q <- 0
    Talpha <- T * 0.025

    for (j in 1:T) {
      if (X_t[j] + var < 0) {
        q <- q + X_t[j] / es
      }
    }
    Z_2 <- q / Talpha + 1
    Z_2_list <- c(Z_2_list, Z_2)
  }

  Z_2_list <- Z_2_list[!is.na(Z_2_list)]
  quant <- quantile(Z_2_list, quan)
  return(list(quant = quant, Z_2_list = Z_2_list, mean = mean(Z_2_list)))
}
```

Test 3
```{r}
Z_3 <- function(quan, T, n, df, dfnull, shift = 0, scale = 1, Norm = FALSE) {
  Z_3_list <- numeric(n)
  TAlpha <- round(T * 0.025)
  
  if (!Norm) {
    integral <- -T / TAlpha * integrate(
      function(p) pbeta(1 - p, shape1 = T - TAlpha, shape2 = TAlpha) * qt(p, df = dfnull),
      lower = 0, upper = 1
    )$value
  } else {
    integral <- -T / TAlpha * integrate(
      function(p) pbeta(1 - p, shape1 = T - TAlpha, shape2 = TAlpha) * qt(p, df = dfnull, ncp = 0) * sqrt((dfnull - 2) / dfnull),
      lower = 0, upper = 1
    )$value
  }
  
  for (j in 1:n) {
    U <- runif(T)
    PU <- if (!Norm) {
      qt(U, df = df) * scale+shift
    } else {
      qt(U, df = df) * sqrt((df - 2) / df)+shift
    }
    
    PU <- sort(PU)
    q <- sum(PU[1:TAlpha])
    ES_hat <- -q / TAlpha
    
    Z_3_list[j] <- -ES_hat / integral + 1
  }
  
  Z_3_list <- Z_3_list[!is.na(Z_3_list)]
  
  quant <- quantile(Z_3_list, quan, na.rm = TRUE)
  mean_Z_3 <- mean(Z_3_list)
  
  return(list(quant = quant, Z_3_list = Z_3_list, mean_Z_3 = mean_Z_3))
}
```

VaR backtest
```{r}
VaR_backtest <- function(sigma2, mu, T, n, shift = 0, scale = 1, df = NULL, dfnull = 100, type = 'Normal', Norm = FALSE, alpha=0.99) {
  var <- VaR(alpha = alpha, sigma2 = sigma2, mu = mu, df = dfnull, type = type, Norm = Norm)
  var_back <- numeric(n)
  
  if (!Norm) {
    for (j in 1:n) {
      X_t <- rt(T, df = df) * scale+shift
      exceedances <- sum(X_t + var < 0)
      var_back[j] <- exceedances
    }
  } else {
    for (j in 1:n) {
      X_t <- rt(T, df = df) * sqrt((df - 2) / df)+shift
      exceedances <- sum(X_t + var < 0)
      var_back[j] <- exceedances
    }
  }
  
  var_back <- var_back[!is.na(var_back)]
  return(list(var_back = var_back, mean_exceedances = mean(var_back)))
}

```

Power
```{r}
Power <- function(data_null, data_alternatives, significance_level, alternative_names = NULL) {
  sorted_null <- sort(data_null)
  crit_value <- quantile(sorted_null, probs = significance_level)

  powers <- list(Significance_Level = significance_level)
  for (i in seq_along(data_alternatives)) {
    power <- mean(data_alternatives[[i]] <= crit_value)
    powers[[alternative_names[i]]] <- power
  }

  return(as.data.frame(powers))
}
```

Til print af tabeller
```{r}
compute_table <- function(data_null, data_alt, significance_levels, df_value, source_name, power_labels) {
  lapply(significance_levels, function(sl) {
    Power(data_null, data_alt, sl, power_labels)
  }) %>%
    bind_rows() %>%
    mutate(H0_df = df_value, Source = source_name)
}
```

Til plots
```{r}
smoothed_cdf <- function(data, n_points = 100000) {
  density_data <- density(data, n = n_points, adjust = 3, from = min(data), to = max(data))
  x <- density_data$x
  y <- cumsum(density_data$y) / sum(density_data$y)
  return(data.frame(value = x, cdf = y))
}

plot_cdfs <- function(null_data, alternatives, critical_value1, critical_value2, pct_5 = NULL, pct_10 = NULL, alternative_names, plot_title = "Comparison of CDF (Null) and 1 - CDF (Alternatives)", xlim=NULL, ylim=NULL) {
  cdf_null <- smoothed_cdf(null_data)
  cdf_null$group <- "Null Hypothesis (CDF)"
  
  cdf_alternatives <- do.call(rbind, lapply(1:length(alternatives), function(i) {
    cdf_data <- smoothed_cdf(alternatives[[i]])
    cdf_data$cdf <- 1 - cdf_data$cdf
    cdf_data$group <- alternative_names[i]
    return(cdf_data)
  }))
  
  combined_data <- rbind(cdf_null, cdf_alternatives)
  
  p <- ggplot(combined_data, aes(x = value, y = cdf, color = group)) +
    geom_line(size = 1) +
    geom_vline(xintercept = critical_value1, linetype = "dashed", color = "black") +
    geom_vline(xintercept = critical_value2, linetype = "dashed", color = "black")+
    geom_vline(xintercept = pct_5, linetype = "dashed", color = "green")+
    geom_vline(xintercept = pct_10, linetype = "dashed", color = "green")+
    labs(
    title = plot_title,
    x = "Value",
    y = "Probability",
    color = "Group"
  ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  if (!is.null(xlim)) {
    p <- p + scale_x_continuous(limits = xlim)
  }
  
  if (!is.null(ylim)) {
    p <- p + scale_y_continuous(limits = ylim)
  }
  
  print(p)
}
```

Standard antagelser
```{r}
T <- 250
n <- 100000
alpha <- 0.975
beta <- 0.99
```

Tabel 1 del 1:
```{r}
T <- 250 # Antal dage
n <- 100000
alpha <- 0.975 # Niveau for Expected Shortfall (97.5%)
beta <- 0.99 # Niveau for Value at Risk (1%)
nu <- c(5, 100) # Frihedsgrader for Student-t-fordeling
alpha_prime <- c(0.975, 0.95, 0.9) # Alpha'-værdier til ES

results <- data.frame()

for (df in nu) {
  
  for (ap in alpha_prime) {
    
    gamma <- ES(alpha, df = df, type = 't', Norm = FALSE)/ES(ap, df = df, type = 't', Norm = FALSE)  
    
    ES_new <- ES(alpha, df = df, type = 't', Norm = FALSE) * gamma
    
    VaR_new1pct <- VaR(alpha=beta, df = df, type = 't', Norm = FALSE) * gamma
    
    results <- rbind(results, data.frame(
      nu = df,
      alpha_prime = ap,
      gamma = gamma,
      VaR_new = VaR_new1pct,
      ES_new = ES_new
    ))
  }
}

print(results)
```
