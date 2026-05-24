
###############################################################################
# Prova 1 - Questão 1: Distribuição Bivariada com Marginais Fréchet
# Disciplina: Métodos Multivariados e Machine Learning
###############################################################################

# Pacotes necessários
if (!require(plotly)) install.packages("plotly")
if (!require(copula)) install.packages("copula")
if (!require(VineCopula)) install.packages("VineCopula")
if (!require(e1071)) install.packages("e1071")
if (!require(numDeriv)) install.packages("numDeriv")
if (!require(stats4)) install.packages("stats4")
if (!require(knitr)) install.packages("knitr")

library(plotly)
library(copula)
library(e1071)
library(numDeriv)

###############################################################################
# (a) Função de densidade e plots
###############################################################################

# CDF Bivariada Fréchet
CDF_biv <- function(x, y, sigma1, sigma2, alpha, rho) {
  x1 <- x / sigma1
  x2 <- y / sigma2
  exp(-(x1^(-alpha)) - (x2^(-alpha)) + rho * ((x1^alpha) + (x2^alpha))^(-1))
}

# Marginal Fréchet CDF
pFrechet <- function(x, sigma, alpha) {
  ifelse(x <= 0, 0, exp(-(x / sigma)^(-alpha)))
}

# Marginal Fréchet PDF
dFrechet <- function(x, sigma, alpha) {
  ifelse(x <= 0, 0,
         (alpha / sigma) * (x / sigma)^(-alpha - 1) * exp(-(x / sigma)^(-alpha)))
}

# PDF Bivariada (Derivação Analítica)
# f(x,y) = F(x,y) * [g_x * g_y + g_xy]
# onde g(x,y) = ln F(x,y)
#
# Derivação:
# Definindo a = (x/s1)^alpha, b = (y/s2)^alpha, S = a + b:
# g(x,y) = -1/a - 1/b + rho/S
# g_x = (alpha/x) * [1/a - rho*a/S^2]
# g_y = (alpha/y) * [1/b - rho*b/S^2]
# g_xy = 2*rho*alpha^2 * a*b / (x*y*S^3)
#
# f(x,y) = F * (g_x * g_y + g_xy)

PDF_biv <- function(x, y, sigma1, sigma2, alpha, rho) {
  x1 <- x / sigma1
  x2 <- y / sigma2
  a <- x1^alpha
  b <- x2^alpha
  S <- a + b
  
  Fxy <- CDF_biv(x, y, sigma1, sigma2, alpha, rho)
  
  gx <- (alpha / x) * (1/a - rho * a / S^2)
  gy <- (alpha / y) * (1/b - rho * b / S^2)
  gxy <- 2 * rho * alpha^2 * a * b / (x * y * S^3)
  
  f <- Fxy * (gx * gy + gxy)
  return(f)
}

# Verificação numérica da PDF (derivada numérica de 2a ordem)
PDF_biv_num <- function(x, y, sigma1, sigma2, alpha, rho) {
  h <- 1e-5
  (CDF_biv(x + h, y + h, sigma1, sigma2, alpha, rho) -
     CDF_biv(x + h, y, sigma1, sigma2, alpha, rho) -
     CDF_biv(x, y + h, sigma1, sigma2, alpha, rho) +
     CDF_biv(x, y, sigma1, sigma2, alpha, rho)) / (h^2)
}

# ---- Plots da densidade ----
cat("=== (a) Gerando plots da densidade ===\n")

# Parâmetros para os plots
params_list <- list(
  list(alpha = 2, sigma1 = 1, sigma2 = 1, rho = 0,   label = "rho = 0 (independência)"),
  list(alpha = 2, sigma1 = 1, sigma2 = 1, rho = 0.5, label = "rho = 0.5"),
  list(alpha = 2, sigma1 = 1, sigma2 = 1, rho = 0.9, label = "rho = 0.9")
)

x_vals <- seq(0.01, 5, length.out = 80)
y_vals <- seq(0.01, 5, length.out = 80)

# Gerando plots para cada valor de rho
for (par_set in params_list) {
  z <- outer(x_vals, y_vals, Vectorize(function(x, y)
    PDF_biv(x, y, par_set$sigma1, par_set$sigma2, par_set$alpha, par_set$rho)))
  
  # Plot de superfície 3D
  p3d <- plot_ly(x = x_vals, y = y_vals, z = z, type = "surface",
                 colorscale = "Viridis") %>%
    layout(title = paste("Densidade Conjunta -", par_set$label),
           scene = list(xaxis = list(title = "x"),
                        yaxis = list(title = "y"),
                        zaxis = list(title = "f(x,y)")))
  print(p3d)
  
  # Plot de contorno
  p_cont <- plot_ly(x = x_vals, y = y_vals, z = z, type = "contour",
                    colorscale = "Earth") %>%
    layout(title = paste("Contornos -", par_set$label),
           xaxis = list(title = "x"), yaxis = list(title = "y"))
  print(p_cont)
}

# Verificação: PDF analítica vs numérica
cat("\nVerificação PDF analítica vs numérica (alpha=2, sigma1=sigma2=1, rho=0.5):\n")
test_pts <- expand.grid(x = c(0.5, 1, 2), y = c(0.5, 1, 2))
for (i in 1:nrow(test_pts)) {
  fa <- PDF_biv(test_pts$x[i], test_pts$y[i], 1, 1, 2, 0.5)
  fn <- PDF_biv_num(test_pts$x[i], test_pts$y[i], 1, 1, 2, 0.5)
  cat(sprintf("  (%.1f, %.1f): analítica = %.6f, numérica = %.6f, diff = %.2e\n",
              test_pts$x[i], test_pts$y[i], fa, fn, abs(fa - fn)))
}

###############################################################################
# (b) Estimadores de Máxima Verossimilhança
###############################################################################

cat("\n=== (b) Estimação por MLE ===\n")

# Log-verossimilhança negativa
neg_loglik <- function(par, x, y) {
  sigma1 <- par[1]
  sigma2 <- par[2]
  alpha  <- par[3]
  rho    <- par[4]
  
  # Restrições nos parâmetros
  if (sigma1 <= 0 || sigma2 <= 0 || alpha <= 0 || rho < 0 || rho > 1)
    return(1e20)
  
  ll <- tryCatch({
    f <- PDF_biv(x, y, sigma1, sigma2, alpha, rho)
    f[f <= 0] <- .Machine$double.xmin
    sum(log(f))
  }, error = function(e) return(-1e20))
  
  if (!is.finite(ll)) return(1e20)
  return(-ll)
}

# Estratégia 1: optim com L-BFGS-B
MLE_LBFGSB <- function(x, y, start = c(1, 1, 1, 0.5)) {
  t0 <- proc.time()
  fit <- optim(par = start,
               fn = neg_loglik,
               x = x, y = y,
               method = "L-BFGS-B",
               lower = c(1e-4, 1e-4, 1e-4, 0),
               upper = c(Inf, Inf, Inf, 1),
               control = list(maxit = 5000))
  tempo <- (proc.time() - t0)[3]
  return(list(par = fit$par, loglik = -fit$value, convergence = fit$convergence,
              method = "L-BFGS-B", time = tempo))
}

# Estratégia 2: optim com Nelder-Mead (reparametrizado)
MLE_NM <- function(x, y, start = c(1, 1, 1, 0.5)) {
  # Reparametrização: sigma = exp(theta), rho = 1/(1+exp(-phi))
  neg_loglik_reparam <- function(theta, x, y) {
    sigma1 <- exp(theta[1])
    sigma2 <- exp(theta[2])
    alpha  <- exp(theta[3])
    rho    <- 1 / (1 + exp(-theta[4]))
    neg_loglik(c(sigma1, sigma2, alpha, rho), x, y)
  }
  
  theta0 <- c(log(start[1]), log(start[2]), log(start[3]),
               log(start[4] / (1 - start[4] + 1e-10)))
  
  t0 <- proc.time()
  fit <- optim(par = theta0,
               fn = neg_loglik_reparam,
               x = x, y = y,
               method = "Nelder-Mead",
               control = list(maxit = 10000))
  tempo <- (proc.time() - t0)[3]
  
  par_est <- c(exp(fit$par[1]), exp(fit$par[2]), exp(fit$par[3]),
               1 / (1 + exp(-fit$par[4])))
  return(list(par = par_est, loglik = -fit$value, convergence = fit$convergence,
              method = "Nelder-Mead", time = tempo))
}

# Estratégia 3: nlminb
MLE_nlminb <- function(x, y, start = c(1, 1, 1, 0.5)) {
  t0 <- proc.time()
  fit <- nlminb(start = start,
                objective = neg_loglik,
                x = x, y = y,
                lower = c(1e-4, 1e-4, 1e-4, 0),
                upper = c(Inf, Inf, Inf, 1),
                control = list(iter.max = 5000))
  tempo <- (proc.time() - t0)[3]
  return(list(par = fit$par, loglik = -fit$objective, convergence = fit$convergence,
              method = "nlminb", time = tempo))
}

# Estratégia 4: PSO (Particle Swarm Optimization) via grid search + optim
MLE_grid_optim <- function(x, y) {
  # Busca em grid para bons valores iniciais
  sigma_grid <- c(0.5, 1, 2, 5)
  alpha_grid <- c(0.5, 1, 2, 3)
  rho_grid   <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  
  best_ll <- Inf
  best_start <- c(1, 1, 1, 0.5)
  
  for (s1 in sigma_grid)
    for (s2 in sigma_grid)
      for (a in alpha_grid)
        for (r in rho_grid) {
          ll <- tryCatch(neg_loglik(c(s1, s2, a, r), x, y),
                         error = function(e) Inf)
          if (ll < best_ll) {
            best_ll <- ll
            best_start <- c(s1, s2, a, r)
          }
        }
  
  t0 <- proc.time()
  fit <- optim(par = best_start,
               fn = neg_loglik,
               x = x, y = y,
               method = "L-BFGS-B",
               lower = c(1e-4, 1e-4, 1e-4, 0),
               upper = c(Inf, Inf, Inf, 1),
               control = list(maxit = 5000))
  tempo <- (proc.time() - t0)[3]
  
  return(list(par = fit$par, loglik = -fit$value, convergence = fit$convergence,
              method = "Grid+L-BFGS-B", time = tempo, start = best_start))
}

###############################################################################
# (c) Aplicação aos dados income-consumption
###############################################################################

cat("\n=== (c) Aplicação aos dados income-consumption ===\n")

# Leitura e preparação dos dados
rfam08   <- read.csv("rfam08.csv")
risfam08 <- read.csv("risfam08.csv")

# Converter se necessário
rfam08$Y   <- as.numeric(rfam08$Y)
risfam08$C <- as.numeric(risfam08$C)

# Renomear Y em rfam08 para evitar conflito
names(rfam08)[names(rfam08) == "Y"] <- "Y_fam"

# Merge
dfam <- merge(rfam08, risfam08, by = "nquest")

# Remover valores não-positivos
dfam2 <- dfam[dfam$C > 0 & dfam$Y_fam > 0, c("Y_fam", "C")]
colnames(dfam2) <- c("income", "consumption")

cat("Dimensão dos dados (após limpeza):", nrow(dfam2), "observações\n")

# Estatísticas descritivas
cat("\nEstatísticas descritivas:\n")
cat("Income - Média:", mean(dfam2$income), " Mediana:", median(dfam2$income), "\n")
cat("Consumption - Média:", mean(dfam2$consumption), " Mediana:", median(dfam2$consumption), "\n")
cat("Correlação:", cor(dfam2$income, dfam2$consumption), "\n")

# Scatterplot
plot(dfam2$income, dfam2$consumption,
     xlab = "Income (Y)", ylab = "Consumption (C)",
     main = "Income vs Consumption", pch = 16, cex = 0.3, col = "steelblue")

# Escalamento para estabilidade numérica (dividir pela mediana)
scale_x <- median(dfam2$income)
scale_y <- median(dfam2$consumption)
x_data <- dfam2$income / scale_x
y_data <- dfam2$consumption / scale_y

cat("\nFatores de escala: income /", scale_x, ", consumption /", scale_y, "\n")

# Estimação MLE - Modelo Completo (4 parâmetros: sigma1, sigma2, alpha, rho)
cat("\n--- Estimação MLE: Modelo Completo ---\n")

# Testar múltiplas estratégias
cat("Estratégia 1: L-BFGS-B\n")
fit1 <- tryCatch(MLE_LBFGSB(x_data, y_data, start = c(1, 1, 1, 0.5)),
                 error = function(e) list(par = NA, loglik = NA, time = NA, method = "L-BFGS-B"))
cat("  Par:", fit1$par, "\n  LogLik:", fit1$loglik, "\n  Tempo:", fit1$time, "s\n")

cat("Estratégia 2: Nelder-Mead\n")
fit2 <- tryCatch(MLE_NM(x_data, y_data, start = c(1, 1, 1, 0.5)),
                 error = function(e) list(par = NA, loglik = NA, time = NA, method = "Nelder-Mead"))
cat("  Par:", fit2$par, "\n  LogLik:", fit2$loglik, "\n  Tempo:", fit2$time, "s\n")

cat("Estratégia 3: nlminb\n")
fit3 <- tryCatch(MLE_nlminb(x_data, y_data, start = c(1, 1, 1, 0.5)),
                 error = function(e) list(par = NA, loglik = NA, time = NA, method = "nlminb"))
cat("  Par:", fit3$par, "\n  LogLik:", fit3$loglik, "\n  Tempo:", fit3$time, "s\n")

cat("Estratégia 4: Grid Search + L-BFGS-B\n")
fit4 <- tryCatch(MLE_grid_optim(x_data, y_data),
                 error = function(e) list(par = NA, loglik = NA, time = NA, method = "Grid+L-BFGS-B"))
cat("  Par:", fit4$par, "\n  LogLik:", fit4$loglik, "\n  Tempo:", fit4$time, "s\n")

# Selecionar o melhor ajuste
fits <- list(fit1, fit2, fit3, fit4)
logliks <- sapply(fits, function(f) ifelse(is.na(f$loglik[1]), -Inf, f$loglik))
best_fit <- fits[[which.max(logliks)]]
cat("\nMelhor ajuste:", best_fit$method, "\n")
cat("  sigma1 =", best_fit$par[1], "\n")
cat("  sigma2 =", best_fit$par[2], "\n")
cat("  alpha  =", best_fit$par[3], "\n")
cat("  rho    =", best_fit$par[4], "\n")
cat("  LogLik =", best_fit$loglik, "\n")

# Parâmetros na escala original
cat("\nParâmetros na escala original:\n")
cat("  sigma1_orig =", best_fit$par[1] * scale_x, "\n")
cat("  sigma2_orig =", best_fit$par[2] * scale_y, "\n")

# Critérios de informação - modelo completo
n <- nrow(dfam2)
k_full <- 4
AIC_full <- -2 * best_fit$loglik + 2 * k_full
BIC_full <- -2 * best_fit$loglik + k_full * log(n)
cat("\n  AIC =", AIC_full, "\n")
cat("  BIC =", BIC_full, "\n")

# ----- Qualidade do ajuste -----

# 1. Teste de Kolmogorov-Smirnov para as marginais
cat("\n--- Testes de qualidade do ajuste ---\n")

# Marginal X (income)
ks_x <- ks.test(x_data, pFrechet, sigma = best_fit$par[1], alpha = best_fit$par[3])
cat("KS test marginal X (income): p-value =", ks_x$p.value, "\n")

# Marginal Y (consumption)
ks_y <- ks.test(y_data, pFrechet, sigma = best_fit$par[2], alpha = best_fit$par[3])
cat("KS test marginal Y (consumption): p-value =", ks_y$p.value, "\n")

# 2. QQ-plots marginais
par(mfrow = c(2, 2))

# QQ-plot marginal X
u_emp_x <- rank(x_data) / (n + 1)
q_teo_x <- (-log(sort(u_emp_x)))^(-1/best_fit$par[3]) * best_fit$par[1]
plot(sort(x_data), q_teo_x, main = "QQ-Plot Marginal X (Income)",
     xlab = "Quantis empíricos", ylab = "Quantis teóricos",
     pch = 16, cex = 0.3, col = "steelblue")
abline(0, 1, col = "red", lwd = 2)

# QQ-plot marginal Y
u_emp_y <- rank(y_data) / (n + 1)
q_teo_y <- (-log(sort(u_emp_y)))^(-1/best_fit$par[3]) * best_fit$par[2]
plot(sort(y_data), q_teo_y, main = "QQ-Plot Marginal Y (Consumption)",
     xlab = "Quantis empíricos", ylab = "Quantis teóricos",
     pch = 16, cex = 0.3, col = "steelblue")
abline(0, 1, col = "red", lwd = 2)

# 3. Contornos da densidade ajustada vs dados
x_grid <- seq(min(x_data) * 0.8, quantile(x_data, 0.95), length.out = 60)
y_grid <- seq(min(y_data) * 0.8, quantile(y_data, 0.95), length.out = 60)
z_fitted <- outer(x_grid, y_grid, Vectorize(function(x, y)
  PDF_biv(x, y, best_fit$par[1], best_fit$par[2], best_fit$par[3], best_fit$par[4])))

plot(x_data, y_data, pch = 16, cex = 0.2, col = "gray70",
     xlab = "Income (escalonado)", ylab = "Consumption (escalonado)",
     main = "Contornos da densidade ajustada")
contour(x_grid, y_grid, z_fitted, add = TRUE, col = "red", nlevels = 15)

par(mfrow = c(1, 1))

###############################################################################
# (d) Modelo reduzido: sigma1 = sigma2
###############################################################################

cat("\n=== (d) Modelo reduzido: sigma1 = sigma2 ===\n")

# Log-verossimilhança negativa com restrição sigma1 = sigma2
neg_loglik_reduced <- function(par, x, y) {
  sigma <- par[1]
  alpha <- par[2]
  rho   <- par[3]
  neg_loglik(c(sigma, sigma, alpha, rho), x, y)
}

# Ajuste do modelo reduzido
fit_red <- nlminb(start = c(1, 1, 0.5),
                  objective = neg_loglik_reduced,
                  x = x_data, y = y_data,
                  lower = c(1e-4, 1e-4, 0),
                  upper = c(Inf, Inf, 1),
                  control = list(iter.max = 5000))

cat("Modelo Reduzido (sigma1 = sigma2):\n")
cat("  sigma =", fit_red$par[1], "\n")
cat("  alpha =", fit_red$par[2], "\n")
cat("  rho   =", fit_red$par[3], "\n")
cat("  LogLik =", -fit_red$objective, "\n")

k_red <- 3
AIC_red <- -2 * (-fit_red$objective) + 2 * k_red
BIC_red <- -2 * (-fit_red$objective) + k_red * log(n)

cat("  AIC =", AIC_red, "\n")
cat("  BIC =", BIC_red, "\n")

# Comparação
cat("\n--- Comparação dos modelos ---\n")
cat(sprintf("%-20s %10s %10s %12s\n", "Modelo", "AIC", "BIC", "LogLik"))
cat(sprintf("%-20s %10.2f %10.2f %12.2f\n", "Completo (4 par)",
            AIC_full, BIC_full, best_fit$loglik))
cat(sprintf("%-20s %10.2f %10.2f %12.2f\n", "Reduzido (3 par)",
            AIC_red, BIC_red, -fit_red$objective))

# Teste de razão de verossimilhanças
LRT <- 2 * (best_fit$loglik - (-fit_red$objective))
p_lrt <- pchisq(LRT, df = 1, lower.tail = FALSE)
cat("\nTeste de Razão de Verossimilhanças:\n")
cat("  LRT =", LRT, "\n")
cat("  p-value =", p_lrt, "\n")
cat("  Conclusão:", ifelse(p_lrt < 0.05,
    "Rejeita H0 -> modelos diferem significativamente, preferir o completo.",
    "Não rejeita H0 -> modelo reduzido é adequado."), "\n")

###############################################################################
# (e) Comparação com cópulas concorrentes
###############################################################################

cat("\n=== (e) Cópulas concorrentes ===\n")

# Transformar dados para uniformes usando as marginais Fréchet estimadas
# Usando os parâmetros do modelo completo
u_data <- pFrechet(x_data, sigma = best_fit$par[1], alpha = best_fit$par[3])
v_data <- pFrechet(y_data, sigma = best_fit$par[2], alpha = best_fit$par[3])

# Verificar se estão em (0,1)
u_data <- pmin(pmax(u_data, 1e-10), 1 - 1e-10)
v_data <- pmin(pmax(v_data, 1e-10), 1 - 1e-10)

pseudo_obs <- cbind(u_data, v_data)

cat("Tau de Kendall empírico:", cor(u_data, v_data, method = "kendall"), "\n")

# --- Ajuste das cópulas ---

# 1. Gumbel-Hougaard
cat("\n1. Cópula de Gumbel-Hougaard:\n")
cop_gumbel <- gumbelCopula(dim = 2)
fit_gumbel <- tryCatch(
  fitCopula(cop_gumbel, pseudo_obs, method = "mpl"),
  error = function(e) {
    cat("  Erro no ajuste, tentando método itau...\n")
    fitCopula(cop_gumbel, pseudo_obs, method = "itau")
  }
)
cat("  theta =", coef(fit_gumbel), "\n")
cat("  LogLik =", logLik(fit_gumbel), "\n")

# 2. Frank
cat("\n2. Cópula de Frank:\n")
cop_frank <- frankCopula(dim = 2)
fit_frank <- tryCatch(
  fitCopula(cop_frank, pseudo_obs, method = "mpl"),
  error = function(e) fitCopula(cop_frank, pseudo_obs, method = "itau"))
cat("  theta =", coef(fit_frank), "\n")
cat("  LogLik =", logLik(fit_frank), "\n")

# 3. Clayton
cat("\n3. Cópula de Clayton:\n")
cop_clayton <- claytonCopula(dim = 2)
fit_clayton <- tryCatch(
  fitCopula(cop_clayton, pseudo_obs, method = "mpl"),
  error = function(e) fitCopula(cop_clayton, pseudo_obs, method = "itau"))
cat("  theta =", coef(fit_clayton), "\n")
cat("  LogLik =", logLik(fit_clayton), "\n")

# 4. Farlie-Gumbel-Morgenstern (FGM)
cat("\n4. Cópula de Farlie-Gumbel-Morgenstern:\n")
cop_fgm <- fgmCopula(dim = 2)
fit_fgm <- tryCatch(
  fitCopula(cop_fgm, pseudo_obs, method = "mpl"),
  error = function(e) fitCopula(cop_fgm, pseudo_obs, method = "itau"))
cat("  theta =", coef(fit_fgm), "\n")
cat("  LogLik =", logLik(fit_fgm), "\n")

# --- Log-verossimilhança do modelo original (cópula implícita) ---
ll_original <- best_fit$loglik
# Subtrair as marginais para obter a log-lik da cópula apenas
ll_marginal_x <- sum(log(dFrechet(x_data, best_fit$par[1], best_fit$par[3])))
ll_marginal_y <- sum(log(dFrechet(y_data, best_fit$par[2], best_fit$par[3])))
ll_copula_original <- ll_original - ll_marginal_x - ll_marginal_y

cat("\n--- Comparação das cópulas ---\n")
cat("Log-verossimilhança da dependência (cópula) apenas:\n")

# Montar tabela comparativa
copula_names <- c("Modelo Original (Fréchet biv.)", "Gumbel-Hougaard", "Frank",
                  "Clayton", "FGM")
copula_lls <- c(ll_copula_original, as.numeric(logLik(fit_gumbel)),
                as.numeric(logLik(fit_frank)), as.numeric(logLik(fit_clayton)),
                as.numeric(logLik(fit_fgm)))
copula_npars <- c(1, 1, 1, 1, 1)  # cada cópula tem 1 parâmetro de dependência
copula_aics <- -2 * copula_lls + 2 * copula_npars
copula_bics <- -2 * copula_lls + copula_npars * log(n)

cat(sprintf("\n%-30s %12s %10s %10s\n", "Cópula", "LogLik", "AIC", "BIC"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (i in seq_along(copula_names)) {
  cat(sprintf("%-30s %12.2f %10.2f %10.2f\n",
              copula_names[i], copula_lls[i], copula_aics[i], copula_bics[i]))
}

# --- Testes de qualidade do ajuste (GoF) para as cópulas ---
cat("\n--- Testes de bondade do ajuste (Cramér-von Mises) ---\n")

gof_gumbel  <- tryCatch(gofCopula(cop_gumbel, pseudo_obs, N = 100, method = "Sn"),
                         error = function(e) list(p.value = NA))
gof_frank   <- tryCatch(gofCopula(cop_frank, pseudo_obs, N = 100, method = "Sn"),
                         error = function(e) list(p.value = NA))
gof_clayton <- tryCatch(gofCopula(cop_clayton, pseudo_obs, N = 100, method = "Sn"),
                         error = function(e) list(p.value = NA))
gof_fgm     <- tryCatch(gofCopula(cop_fgm, pseudo_obs, N = 100, method = "Sn"),
                         error = function(e) list(p.value = NA))

cat(sprintf("  Gumbel-Hougaard: p-value = %s\n",
            ifelse(is.na(gof_gumbel$p.value), "NA", sprintf("%.4f", gof_gumbel$p.value))))
cat(sprintf("  Frank:           p-value = %s\n",
            ifelse(is.na(gof_frank$p.value), "NA", sprintf("%.4f", gof_frank$p.value))))
cat(sprintf("  Clayton:         p-value = %s\n",
            ifelse(is.na(gof_clayton$p.value), "NA", sprintf("%.4f", gof_clayton$p.value))))
cat(sprintf("  FGM:             p-value = %s\n",
            ifelse(is.na(gof_fgm$p.value), "NA", sprintf("%.4f", gof_fgm$p.value))))

# --- Plots comparativos ---
par(mfrow = c(2, 2))

# Gumbel
u_sim <- rCopula(1000, gumbelCopula(coef(fit_gumbel)))
plot(u_sim, pch = 16, cex = 0.3, col = "gray70",
     main = "Gumbel-Hougaard", xlab = "u", ylab = "v")
points(pseudo_obs, pch = 16, cex = 0.2, col = "red")

# Frank
u_sim <- rCopula(1000, frankCopula(coef(fit_frank)))
plot(u_sim, pch = 16, cex = 0.3, col = "gray70",
     main = "Frank", xlab = "u", ylab = "v")
points(pseudo_obs, pch = 16, cex = 0.2, col = "red")

# Clayton
u_sim <- rCopula(1000, claytonCopula(coef(fit_clayton)))
plot(u_sim, pch = 16, cex = 0.3, col = "gray70",
     main = "Clayton", xlab = "u", ylab = "v")
points(pseudo_obs, pch = 16, cex = 0.2, col = "red")

# FGM
u_sim <- rCopula(1000, fgmCopula(coef(fit_fgm)))
plot(u_sim, pch = 16, cex = 0.3, col = "gray70",
     main = "FGM", xlab = "u", ylab = "v")
points(pseudo_obs, pch = 16, cex = 0.2, col = "red")

par(mfrow = c(1, 1))

cat("\n=== Questão 1 finalizada ===\n")
