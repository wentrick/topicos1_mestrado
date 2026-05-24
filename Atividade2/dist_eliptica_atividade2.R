###############################################################################
# Geracao e Analise de Distribuicoes Elipticas
# Relatorio - Aula 8
###############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)

set.seed(42)

###############################################################################
# 1. GERADOR ESFERICO (base comum)
#    Gera n pontos uniformes na esfera S^{d-1} via rejeicao
###############################################################################

rspherical <- function(d, n) {
  resultado <- matrix(NA, nrow = n, ncol = d)
  for (i in 1:n) {
    norma2 <- 2
    while (norma2 > 1) {
      u <- runif(d, -1, 1)
      if (sum(u^2) <= 1) norma2 <- sum(u^2)
    }
    resultado[i, ] <- u / sqrt(norma2)
  }
  resultado
}

###############################################################################
# 2. GERADORES DE DISTRIBUICOES ELIPTICAS
###############################################################################

# 2.1 Normal Multivariada
# Raio: R^2 ~ chi2(d), ou seja, R ~ chi(d)
rmnorm <- function(n, mu, Sigma) {
  d <- length(mu)
  A <- t(chol(Sigma))          # Cholesky inferior: Sigma = A %*% t(A)
  U <- rspherical(d, n)        # n pontos na esfera
  R <- sqrt(rchisq(n, df = d)) # raio ~ chi(d)
  sweep(R * (U %*% t(A)), 2, mu, "+")
}

# 2.2 t de Student Multivariada
# Mistura de escala: Rt = Rn / sqrt(W/nu), W ~ chi2(nu)
rmvt <- function(n, mu, Sigma, nu) {
  d  <- length(mu)
  A  <- t(chol(Sigma))
  U  <- rspherical(d, n)
  Rn <- sqrt(rchisq(n, df = d))   # raio normal
  W  <- rchisq(n, df = nu)        # variavel auxiliar chi2(nu)
  R  <- Rn / sqrt(W / nu)         # raio da t
  sweep(R * (U %*% t(A)), 2, mu, "+")
}

# 2.3 Laplace Simetrica Multivariada
# Mistura de escala: RL = sqrt(W) * Rn, W ~ Exp(1)
rmlaplace <- function(n, mu, Sigma) {
  d  <- length(mu)
  A  <- t(chol(Sigma))
  U  <- rspherical(d, n)
  Rn <- sqrt(rchisq(n, df = d))
  W  <- rexp(n, rate = 1)
  R  <- sqrt(W) * Rn
  sweep(R * (U %*% t(A)), 2, mu, "+")
}

###############################################################################
# 3. PARAMETROS E GERACAO DAS AMOSTRAS
###############################################################################

n <- 1000

# --- Dimensao 2 ---
mu2    <- c(0, 0)
Sigma2 <- matrix(c(1, 0.6, 0.6, 1), 2, 2)

# --- Dimensao 3 ---
mu3    <- c(0, 0, 0)
Sigma3 <- matrix(c(1, 0.5, 0.3,
                   0.5, 1, 0.4,
                   0.3, 0.4, 1), 3, 3)

# Gerar amostras d=2
X_norm2 <- rmnorm(n, mu2, Sigma2)
X_t2    <- rmvt(n, mu2, Sigma2, nu = 3)
X_lap2  <- rmlaplace(n, mu2, Sigma2)

# Gerar amostras d=3
X_norm3 <- rmnorm(n, mu3, Sigma3)
X_t3    <- rmvt(n, mu3, Sigma3, nu = 3)
X_lap3  <- rmlaplace(n, mu3, Sigma3)

###############################################################################
# 4. ANALISE GRAFICA - d = 2
###############################################################################

df_2d <- bind_rows(
  tibble(X1 = X_norm2[,1], X2 = X_norm2[,2], Modelo = "Normal"),
  tibble(X1 = X_t2[,1],    X2 = X_t2[,2],    Modelo = "t de Student (nu = 3)"),
  tibble(X1 = X_lap2[,1],  X2 = X_lap2[,2],  Modelo = "Laplace")
) %>%
  mutate(Modelo = factor(Modelo, levels = c("Normal", "t de Student (nu = 3)", "Laplace")))

ggplot(df_2d, aes(x = X1, y = X2)) +
  geom_point(aes(color = Modelo), alpha = 0.4, size = 0.8) +
  facet_wrap(~ Modelo) +
  coord_fixed(xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_color_manual(values = c("Normal" = "#3366CC",
                                "t de Student (nu = 3)" = "#CC3333",
                                "Laplace" = "#33AA55")) +
  labs(x = expression(X[1]), y = expression(X[2])) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"))

###############################################################################
# 5. ANALISE GRAFICA - d = 3 (plotly)
###############################################################################

p1 <- plot_ly(x = X_norm3[,1], y = X_norm3[,2], z = X_norm3[,3],
              type = "scatter3d", mode = "markers",
              marker = list(size = 2, color = X_norm3[,3],
                            colorscale = "Blues", opacity = 0.6)) %>%
  layout(title = "Normal (d=3)",
         scene = list(xaxis = list(title = "X1"),
                      yaxis = list(title = "X2"),
                      zaxis = list(title = "X3")))

p2 <- plot_ly(x = X_t3[,1], y = X_t3[,2], z = X_t3[,3],
              type = "scatter3d", mode = "markers",
              marker = list(size = 2, color = X_t3[,3],
                            colorscale = "Reds", opacity = 0.6)) %>%
  layout(title = "t de Student (d=3, nu=3)",
         scene = list(xaxis = list(title = "X1"),
                      yaxis = list(title = "X2"),
                      zaxis = list(title = "X3")))

p3 <- plot_ly(x = X_lap3[,1], y = X_lap3[,2], z = X_lap3[,3],
              type = "scatter3d", mode = "markers",
              marker = list(size = 2, color = X_lap3[,3],
                            colorscale = "Greens", opacity = 0.6)) %>%
  layout(title = "Laplace (d=3)",
         scene = list(xaxis = list(title = "X1"),
                      yaxis = list(title = "X2"),
                      zaxis = list(title = "X3")))

p1
p2
p3

###############################################################################
# 6. SEMELHANCAS E DIFERENCAS - Distribuicao radial
###############################################################################

df_raio <- bind_rows(
  tibble(raio = sqrt(rowSums(X_norm2^2)), Modelo = "Normal"),
  tibble(raio = sqrt(rowSums(X_t2^2)),    Modelo = "t de Student (nu = 3)"),
  tibble(raio = sqrt(rowSums(X_lap2^2)),  Modelo = "Laplace")
) %>%
  mutate(Modelo = factor(Modelo, levels = c("Normal", "t de Student (nu = 3)", "Laplace")))

# Histogramas facetados
ggplot(df_raio, aes(x = raio, fill = Modelo)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 color = "white", alpha = 0.7) +
  facet_wrap(~ Modelo) +
  coord_cartesian(xlim = c(0, 10)) +
  scale_fill_manual(values = c("Normal" = "#3366CC",
                               "t de Student (nu = 3)" = "#CC3333",
                               "Laplace" = "#33AA55")) +
  labs(x = "||X||", y = "Densidade") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, face = "bold"))

# Densidades sobrepostas
ggplot(df_raio, aes(x = raio, color = Modelo, fill = Modelo)) +
  geom_density(alpha = 0.2, linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 8)) +
  scale_color_manual(values = c("Normal" = "#3366CC",
                                "t de Student (nu = 3)" = "#CC3333",
                                "Laplace" = "#33AA55")) +
  scale_fill_manual(values = c("Normal" = "#3366CC",
                               "t de Student (nu = 3)" = "#CC3333",
                               "Laplace" = "#33AA55")) +
  labs(x = "||X||", y = "Densidade",
       title = "Comparacao das distribuicoes radiais") +
  theme_minimal() +
  theme(legend.position = "bottom")

###############################################################################
# 7. ESTIMACAO DE PARAMETROS - Simulacao Monte Carlo
###############################################################################

M     <- 200    # numero de replicas
n_sim <- 500    # tamanho de cada amostra
mu0   <- c(0, 0)
Sig0  <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
nu0   <- 5

# Armazenar resultados em lista e empilhar com dplyr
resultados <- list()

for (m in 1:M) {
  # Normal
  X <- rmnorm(n_sim, mu0, Sig0)
  S <- cov(X)
  resultados[[length(resultados) + 1]] <- tibble(
    replica = m, Modelo = "Normal",
    mu1_hat = mean(X[,1]), mu2_hat = mean(X[,2]),
    sig11_hat = S[1,1], sig12_hat = S[1,2]
  )
  
  # t de Student
  X <- rmvt(n_sim, mu0, Sig0, nu = nu0)
  S <- cov(X) * (nu0 - 2) / nu0   # correcao: Var(t) = nu/(nu-2) * Sigma
  resultados[[length(resultados) + 1]] <- tibble(
    replica = m, Modelo = "t de Student (nu = 5)",
    mu1_hat = mean(X[,1]), mu2_hat = mean(X[,2]),
    sig11_hat = S[1,1], sig12_hat = S[1,2]
  )
  
  # Laplace
  X <- rmlaplace(n_sim, mu0, Sig0)
  S <- cov(X) / 2   # correcao: Var(Laplace) = 2*Sigma
  resultados[[length(resultados) + 1]] <- tibble(
    replica = m, Modelo = "Laplace",
    mu1_hat = mean(X[,1]), mu2_hat = mean(X[,2]),
    sig11_hat = S[1,1], sig12_hat = S[1,2]
  )
}

df_est <- bind_rows(resultados) %>%
  mutate(Modelo = factor(Modelo, levels = c("Normal", "t de Student (nu = 5)", "Laplace")))

###############################################################################
# 7.1 Grafico: estimacao de mu_1
###############################################################################

ggplot(df_est, aes(x = mu1_hat, fill = Modelo)) +
  geom_histogram(aes(y = after_stat(density)), bins = 25,
                 color = "white", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linewidth = 0.8) +
  facet_wrap(~ Modelo) +
  coord_cartesian(xlim = c(-0.3, 0.3)) +
  scale_fill_manual(values = c("Normal" = "#3366CC",
                               "t de Student (nu = 5)" = "#CC3333",
                               "Laplace" = "#33AA55")) +
  labs(x = expression(hat(mu)[1]), y = "Densidade",
       title = expression(paste("Distribuicao de ", hat(mu)[1], " ao longo de M = 200 replicas"))) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, face = "bold"))

###############################################################################
# 7.2 Grafico: estimacao de Sigma_11 e Sigma_12
###############################################################################

df_sig <- df_est %>%
  select(replica, Modelo, sig11_hat, sig12_hat) %>%
  pivot_longer(cols = c(sig11_hat, sig12_hat),
               names_to = "Parametro", values_to = "Estimativa") %>%
  mutate(Parametro = recode(Parametro,
                            "sig11_hat" = "Sigma[11]~(verdadeiro == 1)",
                            "sig12_hat" = "Sigma[12]~(verdadeiro == 0.6)"))

ref_lines <- tibble(
  Parametro = c("Sigma[11]~(verdadeiro == 1)", "Sigma[12]~(verdadeiro == 0.6)"),
  valor_ref = c(1, 0.6)
)

ggplot(df_sig, aes(x = Modelo, y = Estimativa, fill = Modelo)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  geom_hline(data = ref_lines, aes(yintercept = valor_ref),
             color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ Parametro, labeller = label_parsed, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#3366CC",
                               "t de Student (nu = 5)" = "#CC3333",
                               "Laplace" = "#33AA55")) +
  labs(y = "Estimativa", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1))

###############################################################################
# 7.3 Tabela resumo
###############################################################################

resumo <- df_est %>%
  group_by(Modelo) %>%
  summarise(
    Vies_mu1    = round(mean(mu1_hat), 4),
    DP_mu1      = round(sd(mu1_hat), 4),
    Media_Sig11 = round(mean(sig11_hat), 4),
    Media_Sig12 = round(mean(sig12_hat), 4),
    .groups = "drop"
  )

print(resumo)