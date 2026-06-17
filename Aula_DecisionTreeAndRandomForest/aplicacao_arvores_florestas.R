# =====================================================================
# Arvores de decisao e florestas aleatorias para classificacao
# Aplicacao com o conjunto Palmer Penguins
# Referencia: Izbicki & dos Santos, Secoes 8.3 e 8.4
# Metodos Multivariados e Machine Learning - PPGEST/UnB
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)      # dplyr, ggplot2, tidyr, purrr, tibble, readr
  library(palmerpenguins) # dados penguins
  library(rpart)          # arvores CART
  library(rpart.plot)     # diagrama da arvore
  library(ranger)         # florestas aleatorias (rapido)
  library(patchwork)      # composicao de graficos
})

set.seed(2026)
theme_set(theme_minimal(base_size = 12))
pal <- c(Adelie = "#E69F00", Chinstrap = "#0072B2", Gentoo = "#009E73")
dir.create("figs", showWarnings = FALSE)

# ---------------------------------------------------------------------
# 0. Funcoes auxiliares
# ---------------------------------------------------------------------

# Acuracia simples
acuracia <- function(obs, pred) mean(obs == pred)

# Matriz de confusao como tibble (linhas = observado, colunas = previsto)
matriz_confusao <- function(obs, pred) {
  table(observado = obs, previsto = pred)
}

# ---------------------------------------------------------------------
# 1. Dados: especie ~ medidas morfometricas (+ ilha, sexo)
# ---------------------------------------------------------------------

dados <- penguins |>
  drop_na() |>                      # remove ~11 linhas com NA
  mutate(species = factor(species))

glimpse(dados)
count(dados, species)

# Particao treino/teste estratificada por especie (70/30)
idx_treino <- dados |>
  mutate(linha = row_number()) |>
  group_by(species) |>
  slice_sample(prop = 0.70) |>
  pull(linha)

treino <- dados |> slice(idx_treino)
teste  <- dados |> slice(-idx_treino)

formula_mod <- species ~ bill_length_mm + bill_depth_mm +
  flipper_length_mm + body_mass_g + island + sex

# ---------------------------------------------------------------------
# 2. Secao 8.3 - Arvore de classificacao (CART)
# ---------------------------------------------------------------------

# Arvore grande, depois podada por validacao cruzada (custo-complexidade)
arvore_cheia <- rpart(
  formula_mod, data = treino, method = "class",
  control = rpart.control(cp = 0, minsplit = 5, xval = 10)
)

# cp que minimiza o erro de validacao cruzada (regra do menor xerror)
cp_otimo <- arvore_cheia$cptable |>
  as_tibble() |>
  slice_min(xerror, n = 1, with_ties = FALSE) |>
  pull(CP)

arvore <- prune(arvore_cheia, cp = cp_otimo)

# Avaliacao no conjunto de teste
pred_arvore <- predict(arvore, teste, type = "class")
acc_arvore  <- acuracia(teste$species, pred_arvore)
cat("\n[Arvore] acuracia (teste):", round(acc_arvore, 3), "\n")
print(matriz_confusao(teste$species, pred_arvore))

# Diagrama da arvore podada
pdf("figs/fig_arvore.pdf", width = 8, height = 5.2)
rpart.plot(
  arvore, type = 2, extra = 104, box.palette = "Blues",
  branch.lty = 3, nn = FALSE, fallen.leaves = TRUE,
  main = "Arvore de classificacao podada (Palmer Penguins)"
)
dev.off()

# Curva de custo-complexidade (erro VC vs cp / tamanho da arvore)
cp_tbl <- as_tibble(arvore_cheia$cptable)
p_cp <- ggplot(cp_tbl, aes(x = nsplit + 1, y = xerror)) +
  geom_line(color = "grey50") +
  geom_pointrange(aes(ymin = xerror - xstd, ymax = xerror + xstd),
                  color = "#0072B2") +
  geom_vline(xintercept = arvore$cptable[nrow(arvore$cptable), "nsplit"] + 1,
             linetype = 2, color = "#D55E00") +
  labs(x = "Numero de folhas", y = "Erro de validacao cruzada",
       title = "Poda por custo-complexidade") +
  theme_minimal(base_size = 13)
plot(p_cp)
ggsave("figs/fig_cp.pdf", p_cp, width = 6.5, height = 4)

# ---------------------------------------------------------------------
# 3. Secao 8.4 - Floresta aleatoria
# ---------------------------------------------------------------------

floresta <- ranger(
  formula_mod, data = treino,
  num.trees = 500, mtry = 2,            # mtry ~ sqrt(p) para classificacao
  importance = "impurity",              # importancia por reducao de Gini
  probability = FALSE, seed = 2026
)

pred_floresta <- predict(floresta, teste)$predictions
acc_floresta  <- acuracia(teste$species, pred_floresta)
cat("\n[Floresta] erro OOB:", round(floresta$prediction.error, 3), "\n")
cat("[Floresta] acuracia (teste):", round(acc_floresta, 3), "\n")
print(matriz_confusao(teste$species, pred_floresta))

# Importancia das variaveis
imp_tbl <- enframe(ranger::importance(floresta),
                   name = "variavel", value = "importancia") |>
  arrange(importancia)

p_imp <- imp_tbl |>
  mutate(variavel = fct_inorder(variavel)) |>
  ggplot(aes(x = importancia, y = variavel)) +
  geom_col(fill = "#0072B2", width = 0.7) +
  labs(x = "Reducao media de impureza (Gini)", y = NULL,
       title = "Importancia das variaveis (floresta aleatoria)") +
  theme_minimal(base_size = 13)
p_imp

ggsave("figs/fig_importancia.pdf", p_imp, width = 6.5, height = 4)

# Erro OOB em funcao do numero de arvores (usa randomForest p/ a curva)
suppressPackageStartupMessages(library(randomForest))
rf_curva <- randomForest(formula_mod, data = treino,
                         ntree = 500, mtry = 2)
oob_tbl <- tibble(arvores = seq_len(nrow(rf_curva$err.rate)),
                  oob = rf_curva$err.rate[, "OOB"])
p_oob <- ggplot(oob_tbl, aes(arvores, oob)) +
  geom_line(color = "#0072B2", linewidth = 0.8) +
  labs(x = "Numero de arvores", y = "Erro OOB",
       title = "Estabilizacao do erro out-of-bag") +
  theme_minimal(base_size = 13)

p_oob 
ggsave("figs/fig_oob.pdf", p_oob, width = 6.5, height = 4)

# ---------------------------------------------------------------------
# 4. Comparacao final
# ---------------------------------------------------------------------

comparacao <- tibble(
  modelo   = c("Arvore unica (CART podada)", "Floresta aleatoria"),
  acuracia = c(acc_arvore, acc_floresta)
)
cat("\n--- Comparacao (acuracia no teste) ---\n")
print(comparacao)

# ---------------------------------------------------------------------
# 5. Ilustracao: medidas de impureza (slide conceitual)
# ---------------------------------------------------------------------

# No de duas classes com proporcao p da classe 1
impureza <- tibble(p = seq(0, 1, length.out = 401)) |>
  mutate(
    Gini       = 2 * p * (1 - p),
    Entropia   = if_else(p %in% c(0, 1), 0,
                         -(p * log2(p) + (1 - p) * log2(1 - p))) / 2,
    Erro       = pmin(p, 1 - p)
  ) |>
  pivot_longer(-p, names_to = "medida", values_to = "valor")

p_imp_curva <- ggplot(impureza, aes(p, valor, color = medida)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(Gini = "#0072B2",
                                Entropia = "#D55E00",
                                Erro = "#009E73")) +
  labs(x = "Proporcao da classe 1 no no (p)", y = "Impureza",
       color = NULL,
       title = "Medidas de impureza (no com duas classes)",
       subtitle = "Entropia reescalada por 1/2 para comparacao") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
p_imp_curva

ggsave("figs/fig_impureza.pdf", p_imp_curva, width = 6.5, height = 4.2)

# ---------------------------------------------------------------------
# 6. Ilustracao: regioes de decisao em 2D (arvore vs floresta)
# ---------------------------------------------------------------------

# Modelo so com 2 variaveis para visualizar o particionamento do plano
f2 <- species ~ bill_length_mm + flipper_length_mm

arvore_2d <- rpart(f2, data = treino, method = "class",
                   control = rpart.control(cp = 0.01))
floresta_2d <- ranger(f2, data = treino, num.trees = 500,
                      mtry = 1, seed = 2026)

grade <- crossing(
  bill_length_mm    = seq(min(dados$bill_length_mm),
                          max(dados$bill_length_mm), length.out = 220),
  flipper_length_mm = seq(min(dados$flipper_length_mm),
                          max(dados$flipper_length_mm), length.out = 220)
)

grade <- grade |>
  mutate(
    arvore   = predict(arvore_2d, grade, type = "class"),
    floresta = predict(floresta_2d, grade)$predictions
  )

plot_regioes <- function(coluna, titulo) {
  ggplot() +
    geom_raster(data = grade,
                aes(bill_length_mm, flipper_length_mm,
                    fill = .data[[coluna]]), alpha = 0.35) +
    geom_point(data = dados,
               aes(bill_length_mm, flipper_length_mm, color = species),
               size = 1.1) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_color_manual(values = pal, name = "Especie") +
    labs(x = "Comprimento do bico (mm)",
         y = "Comprimento da nadadeira (mm)", title = titulo) +
    theme_minimal(base_size = 12)
}

p_fronteiras <- (plot_regioes("arvore", "Arvore unica") +
                 plot_regioes("floresta", "Floresta aleatoria")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("figs/fig_fronteiras.pdf", p_fronteiras, width = 9, height = 4.4)

# EDA: dispersao das duas variaveis usadas na fronteira
p_eda <- ggplot(dados, aes(bill_length_mm, flipper_length_mm,
                           color = species)) +
  geom_point(size = 1.6, alpha = 0.85) +
  scale_color_manual(values = pal, name = "Especie") +
  labs(x = "Comprimento do bico (mm)",
       y = "Comprimento da nadadeira (mm)",
       title = "Palmer Penguins: tres especies, medidas morfometricas") +
  theme_minimal(base_size = 13)
ggsave("figs/fig_eda.pdf", p_eda, width = 6.5, height = 4.2)

cat("\nFiguras salvas em figs/.\n")
