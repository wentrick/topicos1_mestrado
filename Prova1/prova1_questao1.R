###############################################################################
# Prova 1 — Questão 1: Distribuição Bivariada com Marginais Fréchet
###############################################################################

library(tidyverse)
library(copula)
library(patchwork)
set.seed(42)

###############################################################################
# (a) Funções e plots da densidade
###############################################################################

pFrechet <- function(x, s, a) ifelse(x <= 0, 0, exp(-(x/s)^(-a)))

dFrechet <- function(x, s, a) ifelse(x <= 0, 0, a/s*(x/s)^(-a-1)*exp(-(x/s)^(-a)))

CDF_biv <- function(x, y, s1, s2, al, rho) {
  a <- (x/s1)^al; b <- (y/s2)^al
  exp(-1/a - 1/b + rho/(a + b))
}

PDF_biv <- function(x, y, s1, s2, al, rho) {
  a <- (x/s1)^al; b <- (y/s2)^al; S <- a + b
  Fv  <- CDF_biv(x, y, s1, s2, al, rho)
  gx  <- (al/x) * (1/a - rho*a/S^2)
  gy  <- (al/y) * (1/b - rho*b/S^2)
  gxy <- 2*rho*al^2*a*b / (x*y*S^3)
  Fv * (gx*gy + gxy)
}

grid_df <- crossing(x = seq(.05, 5, length.out = 60),
                    y = seq(.05, 5, length.out = 60),
                    rho = c(0, 0.5, 0.9)) %>%
  mutate(f = pmap_dbl(list(x, y, rho),
                      ~ PDF_biv(..1, ..2, s1=1, s2=1, al=2, rho=..3)))

ggplot(grid_df, aes(x, y, z = f)) +
  geom_contour_filled(bins = 12) +
  facet_wrap(~ paste0("rho == ", rho), labeller = label_parsed) +
  labs(title = "Contornos da densidade conjunta (alpha=2, sigma=1)",
       x = "x", y = "y") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none")

###############################################################################
# (b) Dados e MLE
###############################################################################

df_raw <- inner_join(
  read_csv("rfam08.csv",  show_col_types = FALSE) %>% select(nquest, income = Y),
  read_csv("risfam08.csv", show_col_types = FALSE) %>% select(nquest, consumption = C),
  by = "nquest"
) %>%
  mutate(income = as.numeric(income),
         consumption = as.numeric(consumption)) %>%
  filter(income > 0, consumption > 0,
         is.finite(income), is.finite(consumption))

x <- df_raw$income
y <- df_raw$consumption
n <- nrow(df_raw)

neg_ll <- function(par) {
  f <- PDF_biv(x, y, par[1], par[2], par[3], par[4])
  if (any(f <= 0 | !is.finite(f))) return(1e15)
  -sum(log(f))
}

t1 <- system.time(
  fit1 <- optim(c(median(x), median(y), 1, 0.5), neg_ll,
                method = "L-BFGS-B",
                lower = c(1e-4, 1e-4, 0.1, 0),
                upper = c(Inf, Inf, 10, 0.999))
)

neg_ll_repar <- function(phi) {
  par <- c(exp(phi[1]), exp(phi[2]), exp(phi[3]), plogis(phi[4]))
  neg_ll(par)
}
t2 <- system.time(
  fit2 <- optim(c(log(median(x)), log(median(y)), 0, 0), neg_ll_repar)
)
par2 <- c(exp(fit2$par[1:3]), plogis(fit2$par[4]))

t3 <- system.time(
  fit3 <- nlminb(c(median(x), median(y), 1, 0.5), neg_ll,
                 lower = c(1e-4, 1e-4, 0.1, 0),
                 upper = c(Inf, Inf, 10, 0.999))
)

resultados <- tibble(
  metodo   = c("L-BFGS-B", "Nelder-Mead", "nlminb"),
  neg_ll   = c(fit1$value, fit2$value, fit3$objective),
  tempo_s  = c(t1[3], t2[3], t3[3]),
  par_list = list(fit1$par, par2, fit3$par)
)

best_idx <- which.min(resultados$neg_ll)
par_full <- resultados$par_list[[best_idx]]
ll_full  <- -resultados$neg_ll[best_idx]

resultados %>%
  select(Método = metodo, NegLogLik = neg_ll, Tempo_s = tempo_s) %>%
  knitr::kable(digits = c(0, 2, 3))

cat("Parâmetros estimados:\n")
cat("  sigma1 =", par_full[1], "\n")
cat("  sigma2 =", par_full[2], "\n")
cat("  alpha  =", par_full[3], "\n")
cat("  rho    =", par_full[4], "\n")
cat("  LogLik =", ll_full, "\n")

###############################################################################
# (c) Qualidade do ajuste
###############################################################################

df_pit <- df_raw %>%
  mutate(u_pit = pFrechet(income,      par_full[1], par_full[3]),
         v_pit = pFrechet(consumption, par_full[2], par_full[3]))

ks_x <- ks.test(df_pit$u_pit, "punif")
ks_y <- ks.test(df_pit$v_pit, "punif")

p1 <- ggplot(df_pit, aes(sample = u_pit)) +
  stat_qq(distribution = qunif, size = .3, alpha = .3) +
  geom_abline(color = "red") +
  labs(title = "QQ Unif — Income") +
  theme_minimal(base_size = 9)

p2 <- ggplot(df_pit, aes(sample = v_pit)) +
  stat_qq(distribution = qunif, size = .3, alpha = .3) +
  geom_abline(color = "red") +
  labs(title = "QQ Unif — Consumption") +
  theme_minimal(base_size = 9)

grid_contorno <- crossing(
  x = seq(min(x), quantile(x, .99), length.out = 50),
  y = seq(min(y), quantile(y, .99), length.out = 50)
) %>%
  mutate(f = PDF_biv(x, y, par_full[1], par_full[2], par_full[3], par_full[4]))

p3 <- ggplot(df_raw, aes(income, consumption)) +
  geom_point(size = .3, alpha = .3) +
  geom_contour(data = grid_contorno, aes(x = x, y = y, z = f),
               color = "red", linewidth = .4) +
  labs(title = "Dados + contornos ajustados",
       x = "Income", y = "Consumption") +
  theme_minimal(base_size = 9)

p1 + p2 + p3

tibble(Marginal = c("Income (X)", "Consumption (Y)"),
       p_valor_KS = c(ks_x$p.value, ks_y$p.value)) %>%
  knitr::kable(digits = 4)

###############################################################################
# (d) Modelo reduzido sigma1 = sigma2
###############################################################################

neg_ll_red <- function(par) {
  f <- PDF_biv(x, y, par[1], par[1], par[2], par[3])
  if (any(f <= 0 | !is.finite(f))) return(1e15)
  -sum(log(f))
}

fit_red <- nlminb(c(mean(par_full[1:2]), par_full[3], par_full[4]), neg_ll_red,
                  lower = c(1e-4, 0.1, 0), upper = c(Inf, 10, 0.999))
par_red <- fit_red$par
ll_red  <- -fit_red$objective

D     <- 2 * (ll_full - ll_red)
p_lrt <- pchisq(D, df = 1, lower.tail = FALSE)

tab_modelos <- tibble(
  Modelo = c("Completo (4 par)", "Reduzido (3 par)"),
  k      = c(4, 3),
  LogLik = c(ll_full, ll_red)
) %>%
  mutate(AIC = -2*LogLik + 2*k,
         BIC = -2*LogLik + k*log(n))

knitr::kable(tab_modelos, digits = 2)
cat("LRT: D =", round(D, 3), ", p-valor =", round(p_lrt, 4), "\n")

###############################################################################
# (e) Cópulas concorrentes
###############################################################################

df_pseudo <- df_raw %>%
  mutate(u = pFrechet(income,      par_full[1], par_full[3]),
         v = pFrechet(consumption, par_full[2], par_full[3]),
         u = pmin(pmax(u, 1e-10), 1 - 1e-10),
         v = pmin(pmax(v, 1e-10), 1 - 1e-10))

pobs <- df_pseudo %>% select(u, v) %>% as.matrix()

fit_gu <- fitCopula(gumbelCopula(),  pobs, method = "mpl")
fit_fr <- fitCopula(frankCopula(),   pobs, method = "mpl")
fit_cl <- fitCopula(claytonCopula(), pobs, method = "mpl")
fit_fg <- fitCopula(fgmCopula(),     pobs, method = "itau")

ll_cop_orig <- ll_full -
  sum(log(dFrechet(x, par_full[1], par_full[3]))) -
  sum(log(dFrechet(y, par_full[2], par_full[3])))

tab_copulas <- tibble(
  Copula = c("Original (Fréchet biv.)", "Gumbel-Hougaard",
             "Frank", "Clayton", "FGM"),
  Param  = c(par_full[4], coef(fit_gu), coef(fit_fr),
             coef(fit_cl), coef(fit_fg)),
  LogLik = c(ll_cop_orig, as.numeric(logLik(fit_gu)),
             as.numeric(logLik(fit_fr)), as.numeric(logLik(fit_cl)),
             as.numeric(logLik(fit_fg)))
) %>%
  mutate(AIC = -2*LogLik + 2,
         BIC = -2*LogLik + log(n))

knitr::kable(tab_copulas, digits = 2)

copulas_ajustadas <- list(
  Gumbel  = gumbelCopula(coef(fit_gu)),
  Frank   = frankCopula(coef(fit_fr)),
  Clayton = claytonCopula(coef(fit_cl)),
  FGM     = fgmCopula(coef(fit_fg))
)

df_sims <- imap_dfr(copulas_ajustadas, function(cop, nome) {
  rCopula(800, cop) %>%
    as_tibble(.name_repair = ~ c("u", "v")) %>%
    mutate(copula = nome)
})

df_obs <- df_pseudo %>% select(u, v)

ggplot(df_sims, aes(u, v)) +
  geom_point(size = .2, color = "grey70") +
  geom_point(data = df_obs, size = .2, color = "red", alpha = .3) +
  facet_wrap(~ copula, nrow = 1) +
  labs(title = "Simulações (cinza) vs. dados (vermelho)",
       x = "u", y = "v") +
  theme_minimal(base_size = 8)
