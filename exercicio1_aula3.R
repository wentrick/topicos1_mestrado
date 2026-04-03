###############################################################################
## Exercício 1 - Aula 3: Métodos Multivariados e Machine Learning
## Distribuição de Morgenstern bivariada com marginais Gamma (Cópula FGM)
## Felipe Quintino (UnB) - 01/2026
###############################################################################

# ============================================================================
# ITEM (a): Função de densidade correspondente à CDF do slide 15
# ============================================================================

# A CDF bivariada (slide 14, eq. 2) é:
#   F_{X,Y}(x,y) = F_X(x)*F_Y(y)*[1 + rho*(1 - F_X(x))*(1 - F_Y(y))]
#
# onde F_X ~ Gamma(alpha, theta) e F_Y ~ Gamma(beta, theta), -1 <= rho <= 1.
#
# A cópula bivariada (slide 15) é a FGM:
#   C(u1, u2) = u1*u2*[1 + rho*(1 - u1)*(1 - u2)]
#
# A PDF bivariada é obtida pela derivada parcial mista:
#   f_{X,Y}(x,y) = d^2 F / (dx dy)
#
# Resultado:
#   f_{X,Y}(x,y) = f_X(x)*f_Y(y)*[1 + rho*(1 - 2*F_X(x))*(1 - 2*F_Y(y))]
#
# Demonstração:
#   F = u*v*[1 + rho*(1-u)*(1-v)]  onde u=F_X(x), v=F_Y(y)
#   dF/du = v*[1 + rho*(1-2u)*(1-v)]
#   d^2F/(du dv) = 1 + rho*(1-2u)*(1-2v)
#   Pela regra da cadeia: f(x,y) = f_X(x)*f_Y(y)*[1 + rho*(1-2u)*(1-2v)]

# --- Funções do modelo ---

# CDF bivariada
cdf_morgenstern <- function(x, y, alpha, beta, theta, rho) {
  u <- pgamma(x, shape = alpha, rate = 1/theta)
  v <- pgamma(y, shape = beta, rate = 1/theta)
  u * v * (1 + rho * (1 - u) * (1 - v))
}

# PDF bivariada
pdf_morgenstern <- function(x, y, alpha, beta, theta, rho) {
  fx <- dgamma(x, shape = alpha, rate = 1/theta)
  fy <- dgamma(y, shape = beta, rate = 1/theta)
  Fx <- pgamma(x, shape = alpha, rate = 1/theta)
  Fy <- pgamma(y, shape = beta, rate = 1/theta)
  fx * fy * (1 + rho * (1 - 2 * Fx) * (1 - 2 * Fy))
}

# --- Plots da PDF bivariada ---

# Configuração de grade
n_grid <- 80

# Parâmetros de exemplo
params_list <- list(
  list(alpha = 2, beta = 3, theta = 1, rho = 0.8, title = expression(paste(alpha, "=2, ", beta, "=3, ", theta, "=1, ", rho, "=0.8"))),
  list(alpha = 2, beta = 3, theta = 1, rho = -0.8, title = expression(paste(alpha, "=2, ", beta, "=3, ", theta, "=1, ", rho, "=-0.8"))),
  list(alpha = 5, beta = 5, theta = 0.5, rho = 0.5, title = expression(paste(alpha, "=5, ", beta, "=5, ", theta, "=0.5, ", rho, "=0.5"))),
  list(alpha = 1, beta = 2, theta = 2, rho = 0, title = expression(paste(alpha, "=1, ", beta, "=2, ", theta, "=2, ", rho, "=0 (indep.)")))
)

# Plot: contornos da PDF bivariada
pdf("plots_item_a.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (p in params_list) {
  xmax <- qgamma(0.995, shape = p$alpha, rate = 1/p$theta)
  ymax <- qgamma(0.995, shape = p$beta, rate = 1/p$theta)
  
  x_seq <- seq(0.01, xmax, length.out = n_grid)
  y_seq <- seq(0.01, ymax, length.out = n_grid)
  
  z_mat <- outer(x_seq, y_seq, function(x, y) {
    mapply(pdf_morgenstern, x, y,
           MoreArgs = list(alpha = p$alpha, beta = p$beta,
                           theta = p$theta, rho = p$rho))
  })
  
  contour(x_seq, y_seq, z_mat, nlevels = 20,
          xlab = "x", ylab = "y",
          main = p$title, col = heat.colors(20))
  title(sub = "Contornos da PDF bivariada", cex.sub = 0.8)
}
dev.off()

# Plot 3D (perspectiva)
pdf("plots_item_a_3d.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))

for (p in params_list) {
  xmax <- qgamma(0.995, shape = p$alpha, rate = 1/p$theta)
  ymax <- qgamma(0.995, shape = p$beta, rate = 1/p$theta)
  
  x_seq <- seq(0.01, xmax, length.out = n_grid)
  y_seq <- seq(0.01, ymax, length.out = n_grid)
  
  z_mat <- outer(x_seq, y_seq, function(x, y) {
    mapply(pdf_morgenstern, x, y,
           MoreArgs = list(alpha = p$alpha, beta = p$beta,
                           theta = p$theta, rho = p$rho))
  })
  
  persp(x_seq, y_seq, z_mat, theta = 30, phi = 25,
        xlab = "x", ylab = "y", zlab = "f(x,y)",
        main = p$title, col = "lightblue", shade = 0.3,
        ticktype = "detailed", cex.main = 0.9)
}
dev.off()

cat("Item (a): Plots salvos em plots_item_a.pdf e plots_item_a_3d.pdf\n")


# ============================================================================
# ITEM (b): Estimadores de Máxima Verossimilhança (EMV)
# ============================================================================

# A log-verossimilhança para uma amostra (x1,y1),...,(xn,yn) é:
#
#   l(alpha, beta, theta, rho) = sum_i { log f_X(xi) + log f_Y(yi)
#                                  + log[1 + rho*(1 - 2*F_X(xi))*(1 - 2*F_Y(yi))] }
#
# Não é possível obter os EMVs explicitamente pois:
# 1) A Gamma não possui EMV em forma fechada para o parâmetro de forma (alpha, beta)
# 2) O termo da cópula FGM acopla todos os parâmetros de forma não-linear
#
# Estratégias numéricas implementadas:
#   1) optim com método Nelder-Mead (derivada-free)
#   2) optim com método L-BFGS-B (quasi-Newton com restrições de box)
#   3) optim com método BFGS (quasi-Newton sem restrições, via reparametrização)

# Função de log-verossimilhança negativa
neg_loglik <- function(par, x, y) {
  alpha <- par[1]
  beta  <- par[2]
  theta <- par[3]
  rho   <- par[4]
  
  # Restrições de domínio
  if (alpha <= 0 || beta <= 0 || theta <= 0 || abs(rho) > 1) return(1e10)
  
  fx <- dgamma(x, shape = alpha, rate = 1/theta, log = TRUE)
  fy <- dgamma(y, shape = beta, rate = 1/theta, log = TRUE)
  Fx <- pgamma(x, shape = alpha, rate = 1/theta)
  Fy <- pgamma(y, shape = beta, rate = 1/theta)
  
  copula_term <- 1 + rho * (1 - 2 * Fx) * (1 - 2 * Fy)
  
  # Verificar se o termo da cópula é positivo (condição de densidade)
  if (any(copula_term <= 0)) return(1e10)
  
  ll <- sum(fx + fy + log(copula_term))
  return(-ll)
}

# Reparametrização para otimização irrestrita (método BFGS)
# alpha = exp(phi1), beta = exp(phi2), theta = exp(phi3), rho = tanh(phi4)
neg_loglik_repar <- function(phi, x, y) {
  par <- c(exp(phi[1]), exp(phi[2]), exp(phi[3]), tanh(phi[4]))
  neg_loglik(par, x, y)
}

# Função para estimar com múltiplas estratégias
fit_morgenstern <- function(x, y, verbose = TRUE) {
  
  # Valores iniciais baseados nos momentos marginais
  mean_x <- mean(x); var_x <- var(x)
  mean_y <- mean(y); var_y <- var(y)
  
  theta0 <- (var_x / mean_x + var_y / mean_y) / 2
  alpha0 <- mean_x / theta0
  beta0  <- mean_y / theta0
  rho0   <- 0.1 * sign(cor(x, y))
  
  par0 <- c(alpha0, beta0, theta0, rho0)
  
  results <- list()
  
  # --- Estratégia 1: Nelder-Mead ---
  t1 <- system.time({
    fit1 <- tryCatch(
      optim(par0, neg_loglik, x = x, y = y, method = "Nelder-Mead",
            control = list(maxit = 10000, reltol = 1e-10)),
      error = function(e) NULL
    )
  })
  if (!is.null(fit1)) {
    results[["Nelder-Mead"]] <- list(
      par = fit1$par, loglik = -fit1$value,
      convergence = fit1$convergence, time = t1["elapsed"]
    )
  }
  
  # --- Estratégia 2: L-BFGS-B (com box constraints) ---
  t2 <- system.time({
    fit2 <- tryCatch(
      optim(par0, neg_loglik, x = x, y = y, method = "L-BFGS-B",
            lower = c(1e-4, 1e-4, 1e-4, -0.999),
            upper = c(Inf, Inf, Inf, 0.999),
            control = list(maxit = 10000)),
      error = function(e) NULL
    )
  })
  if (!is.null(fit2)) {
    results[["L-BFGS-B"]] <- list(
      par = fit2$par, loglik = -fit2$value,
      convergence = fit2$convergence, time = t2["elapsed"]
    )
  }
  
  # --- Estratégia 3: BFGS com reparametrização ---
  phi0 <- c(log(alpha0), log(beta0), log(theta0), atanh(max(-0.99, min(0.99, rho0))))
  
  t3 <- system.time({
    fit3 <- tryCatch(
      optim(phi0, neg_loglik_repar, x = x, y = y, method = "BFGS",
            control = list(maxit = 10000, reltol = 1e-10)),
      error = function(e) NULL
    )
  })
  if (!is.null(fit3)) {
    par3 <- c(exp(fit3$par[1]), exp(fit3$par[2]), exp(fit3$par[3]), tanh(fit3$par[4]))
    results[["BFGS-repar"]] <- list(
      par = par3, loglik = -fit3$value,
      convergence = fit3$convergence, time = t3["elapsed"]
    )
  }
  
  # Comparação
  if (verbose) {
    cat("\n=== Comparação de Estratégias de Otimização ===\n")
    cat(sprintf("%-15s %-40s %-15s %-10s %-10s\n",
                "Método", "Parâmetros (alpha, beta, theta, rho)",
                "log-lik", "Conv.", "Tempo(s)"))
    cat(paste(rep("-", 90), collapse = ""), "\n")
    for (name in names(results)) {
      r <- results[[name]]
      cat(sprintf("%-15s (%.4f, %.4f, %.4f, %.4f) %15.4f %5d %10.4f\n",
                  name, r$par[1], r$par[2], r$par[3], r$par[4],
                  r$loglik, r$convergence, r$time))
    }
  }
  
  # Retornar o melhor resultado
  best <- which.max(sapply(results, function(r) r$loglik))
  return(list(best = results[[best]], all = results, best_method = names(results)[best]))
}

# --- Teste com dados simulados para validar ---
set.seed(42)
n_sim <- 500
alpha_true <- 3; beta_true <- 2; theta_true <- 1; rho_true <- 0.6

# Geração via cópula FGM (algoritmo do slide 13)
u2_sim <- runif(n_sim)
v_sim  <- runif(n_sim)
a_sim  <- rho_true * (1 - 2 * u2_sim)
# Inversa condicional (slide 16)
u1_sim <- ((1 + a_sim) - sqrt((1 + a_sim)^2 - 4 * a_sim * v_sim)) / (2 * a_sim)
# Tratar caso rho = 0 (a = 0)
u1_sim[a_sim == 0] <- v_sim[a_sim == 0]

x_sim <- qgamma(u1_sim, shape = alpha_true, rate = 1/theta_true)
y_sim <- qgamma(u2_sim, shape = beta_true, rate = 1/theta_true)

cat("\n--- Validação com dados simulados ---\n")
cat(sprintf("Parâmetros verdadeiros: alpha=%.1f, beta=%.1f, theta=%.1f, rho=%.1f\n",
            alpha_true, beta_true, theta_true, rho_true))
fit_sim <- fit_morgenstern(x_sim, y_sim)
cat(sprintf("Melhor método: %s\n", fit_sim$best_method))


# ============================================================================
# ITEM (c): Aplicação aos dados de income-consumption
# ============================================================================

cat("\n\n=== ITEM (c): Dados de income-consumption ===\n")

# Carregar dados
rfam  <- read.csv("rfam08.csv")    # Y = renda
risfam <- read.csv("risfam08.csv") # C = consumo

# Merge pelo identificador da família
dados <- merge(rfam[, c("nquest", "Y")], risfam[, c("nquest", "C")], by = "nquest")
colnames(dados) <- c("nquest", "income", "consumption")

# Remover entradas com valores negativos ou nulos
dados <- dados[dados$income > 0 & dados$consumption > 0, ]
cat(sprintf("Número de famílias após limpeza: %d\n", nrow(dados)))

# Estatísticas descritivas
cat("\n--- Estatísticas descritivas ---\n")
cat(sprintf("Renda:    min=%.2f, mediana=%.2f, média=%.2f, max=%.2f, sd=%.2f\n",
            min(dados$income), median(dados$income), mean(dados$income),
            max(dados$income), sd(dados$income)))
cat(sprintf("Consumo:  min=%.2f, mediana=%.2f, média=%.2f, max=%.2f, sd=%.2f\n",
            min(dados$consumption), median(dados$consumption), mean(dados$consumption),
            max(dados$consumption), sd(dados$consumption)))
cat(sprintf("Correlação de Pearson: %.4f\n", cor(dados$income, dados$consumption)))

# Escalonar os dados para evitar problemas numéricos (dividir por 1000)
x_data <- dados$consumption / 1000
y_data <- dados$income / 1000

cat(sprintf("\nDados escalados (/ 1000):\n"))
cat(sprintf("Consumo/1000: média=%.4f, sd=%.4f\n", mean(x_data), sd(x_data)))
cat(sprintf("Renda/1000:   média=%.4f, sd=%.4f\n", mean(y_data), sd(y_data)))

# Ajuste do modelo FGM-Gamma
cat("\n--- Ajustando modelo de Morgenstern bivariado com marginais Gamma ---\n")
fit_data <- fit_morgenstern(x_data, y_data)

# Parâmetros estimados
par_hat <- fit_data$best$par
alpha_hat <- par_hat[1]
beta_hat  <- par_hat[2]
theta_hat <- par_hat[3]
rho_hat   <- par_hat[4]

cat(sprintf("\nParâmetros estimados (melhor método: %s):\n", fit_data$best_method))
cat(sprintf("  alpha = %.6f\n", alpha_hat))
cat(sprintf("  beta  = %.6f\n", beta_hat))
cat(sprintf("  theta = %.6f\n", theta_hat))
cat(sprintf("  rho   = %.6f\n", rho_hat))
cat(sprintf("  log-verossimilhança = %.4f\n", fit_data$best$loglik))

# --- Critérios de informação ---
n <- nrow(dados)
k <- 4  # número de parâmetros
ll_max <- fit_data$best$loglik

AIC_val <- -2 * ll_max + 2 * k
BIC_val <- -2 * ll_max + k * log(n)

cat(sprintf("\n--- Critérios de Informação ---\n"))
cat(sprintf("  AIC = %.4f\n", AIC_val))
cat(sprintf("  BIC = %.4f\n", BIC_val))

# --- Análise da qualidade do ajuste ---

# 1) Verificação das marginais: QQ-plot contra Gamma ajustada
pdf("plots_item_c.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))

# QQ-plot marginal X (consumo)
q_x <- sort(x_data)
p_x <- (1:n) / (n + 1)
q_x_theo <- qgamma(p_x, shape = alpha_hat, rate = 1/theta_hat)
plot(q_x_theo, q_x, pch = ".", col = "blue",
     xlab = "Quantis teóricos Gamma", ylab = "Quantis amostrais",
     main = "QQ-plot: Consumo vs Gamma")
abline(0, 1, col = "red", lwd = 2)

# QQ-plot marginal Y (renda)
q_y <- sort(y_data)
q_y_theo <- qgamma(p_x, shape = beta_hat, rate = 1/theta_hat)
plot(q_y_theo, q_y, pch = ".", col = "blue",
     xlab = "Quantis teóricos Gamma", ylab = "Quantis amostrais",
     main = "QQ-plot: Renda vs Gamma")
abline(0, 1, col = "red", lwd = 2)

# 2) Histograma marginal do consumo com densidade ajustada
hist(x_data, breaks = 60, prob = TRUE, col = "lightgray",
     main = "Histograma: Consumo", xlab = "Consumo/1000", xlim = c(0, quantile(x_data, 0.99)))
curve(dgamma(x, shape = alpha_hat, rate = 1/theta_hat), add = TRUE, col = "red", lwd = 2)
legend("topright", "Gamma ajustada", col = "red", lwd = 2, cex = 0.8)

# 3) Histograma marginal da renda com densidade ajustada
hist(y_data, breaks = 60, prob = TRUE, col = "lightgray",
     main = "Histograma: Renda", xlab = "Renda/1000", xlim = c(0, quantile(y_data, 0.99)))
curve(dgamma(x, shape = beta_hat, rate = 1/theta_hat), add = TRUE, col = "blue", lwd = 2)
legend("topright", "Gamma ajustada", col = "blue", lwd = 2, cex = 0.8)

# 4) Scatterplot com contornos do modelo ajustado
x_range <- seq(0.01, quantile(x_data, 0.98), length.out = 80)
y_range <- seq(0.01, quantile(y_data, 0.98), length.out = 80)

z_fit <- outer(x_range, y_range, function(x, y) {
  mapply(pdf_morgenstern, x, y,
         MoreArgs = list(alpha = alpha_hat, beta = beta_hat,
                         theta = theta_hat, rho = rho_hat))
})

plot(x_data, y_data, pch = ".", col = rgb(0, 0, 0, 0.1),
     xlab = "Consumo/1000", ylab = "Renda/1000",
     main = "Dados e contornos do modelo ajustado",
     xlim = range(x_range), ylim = range(y_range))
contour(x_range, y_range, z_fit, nlevels = 15, add = TRUE, col = "red", lwd = 1.5)

# 5) CDF empírica bivariada vs teórica (avaliada em pontos da grade)
# Teste de Kolmogorov-Smirnov para cada marginal
ks_x <- ks.test(x_data, "pgamma", shape = alpha_hat, rate = 1/theta_hat)
ks_y <- ks.test(y_data, "pgamma", shape = beta_hat, rate = 1/theta_hat)

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.8, "Testes de Kolmogorov-Smirnov", cex = 1.3, font = 2)
text(0.5, 0.65, sprintf("Consumo: D=%.4f, p-valor=%.4f", ks_x$statistic, ks_x$p.value), cex = 1.1)
text(0.5, 0.50, sprintf("Renda:   D=%.4f, p-valor=%.4f", ks_y$statistic, ks_y$p.value), cex = 1.1)
text(0.5, 0.30, sprintf("rho estimado = %.4f", rho_hat), cex = 1.1)
text(0.5, 0.15, paste("Nota: Cópula FGM tem limitação |rho| <= 1,\n",
                       "mas a correlação de Kendall fica em [-2/9, 2/9]."),
     cex = 0.9)

dev.off()

cat("\nPlots salvos em plots_item_c.pdf\n")

# --- Discussão sobre a adequação do modelo ---
cat("\n--- Discussão sobre a adequação do modelo ---\n")
cat(sprintf("Correlação de Pearson nos dados: %.4f\n", cor(x_data, y_data)))
cat(sprintf("Correlação de Kendall nos dados: %.4f\n", cor(x_data, y_data, method = "kendall")))
cat(sprintf("\nLimitação da cópula FGM:\n"))
cat(sprintf("  A cópula FGM captura apenas dependência fraca.\n"))
cat(sprintf("  O tau de Kendall teórico fica em [-2/9, 2/9] ≈ [-0.222, 0.222].\n"))
cat(sprintf("  Se a dependência nos dados for mais forte, o modelo pode não ser adequado.\n"))

tau_data <- cor(x_data, y_data, method = "kendall")
if (abs(tau_data) > 2/9) {
  cat(sprintf("\n  ATENÇÃO: tau de Kendall dos dados (%.4f) excede o limite da FGM (2/9 ≈ 0.222).\n", tau_data))
  cat(sprintf("  O modelo FGM-Gamma pode NÃO ser adequado para estes dados.\n"))
  cat(sprintf("  Recomenda-se considerar cópulas com dependência mais forte\n"))
  cat(sprintf("  (e.g., Clayton, Frank, Gumbel) ou o modelo Unit-Fréchet do artigo.\n"))
} else {
  cat(sprintf("\n  O tau de Kendall dos dados (%.4f) está dentro do limite da FGM.\n", tau_data))
  cat(sprintf("  O modelo pode ser adequado.\n"))
}

cat("\n=== Fim do Exercício 1 ===\n")
