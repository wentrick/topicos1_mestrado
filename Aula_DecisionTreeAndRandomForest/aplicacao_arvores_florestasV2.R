# Arvores de decisao e florestas aleatorias para classificacao
# Palmer Penguins - Izbicki & dos Santos, Secoes 8.3 e 8.4

library(tidyverse)
library(palmerpenguins)
library(rpart)
library(rpart.plot)
library(ranger)

set.seed(2026)
theme_set(theme_minimal(base_size = 13))
cores <- c(Adelie = "#E69F00", Chinstrap = "#0072B2", Gentoo = "#009E73")
dir.create("figs", showWarnings = FALSE)

# Dados e particao treino/teste (70/30)
dados <- penguins |> drop_na() |> select(-year)
i_treino <- sample(nrow(dados), 0.7 * nrow(dados))
treino <- dados |> slice(i_treino)
teste  <- dados |> slice(-i_treino)

# Visao geral: duas variaveis mais informativas
ggplot(dados, aes(bill_length_mm, flipper_length_mm, color = species)) +
  geom_point() +
  scale_color_manual(values = cores, name = "Especie") +
  labs(x = "Comprimento do bico (mm)", y = "Comprimento da nadadeira (mm)")

ggsave("figs/fig_eda.pdf", width = 6.5, height = 4.2)

# ---- 8.3 Arvore de classificacao ----
arvore <- rpart(species ~ ., data = treino, method = "class")

rpart.plot(arvore, type = 2, extra = 104, box.palette = "Blues")  # mostra
pdf("figs/fig_arvore.pdf", width = 8, height = 5)                 # salva
rpart.plot(arvore, type = 2, extra = 104, box.palette = "Blues")
dev.off()

pred_arvore <- predict(arvore, teste, type = "class")
mean(pred_arvore == teste$species)
table(observado = teste$species, previsto = pred_arvore)

# ---- 8.4 Floresta aleatoria (mtry ~ sqrt(p) por padrao) ----
floresta <- ranger(species ~ ., data = treino, importance = "impurity")

floresta$prediction.error                       # erro out-of-bag
pred_floresta <- predict(floresta, teste)$predictions
mean(pred_floresta == teste$species)
table(observado = teste$species, previsto = pred_floresta)

# Importancia das variaveis
importancia <- enframe(importance(floresta), "variavel", "valor")
ggplot(importancia, aes(valor, fct_reorder(variavel, valor))) +
  geom_col(fill = "#0072B2") +
  labs(x = "Reducao de impureza (Gini)", y = NULL)
ggsave("figs/fig_importancia.pdf", width = 6.5, height = 4)

# ---- Fronteiras de decisao em 2D (arvore vs floresta) ----
arvore_2d   <- rpart(species ~ bill_length_mm + flipper_length_mm,
                     data = treino, method = "class")
floresta_2d <- ranger(species ~ bill_length_mm + flipper_length_mm, data = treino)

grade <- crossing(
  bill_length_mm    = seq(32, 60, length.out = 200),
  flipper_length_mm = seq(170, 232, length.out = 200)
)
grade <- grade |>
  mutate(`Arvore unica`       = predict(arvore_2d, grade, type = "class"),
         `Floresta aleatoria` = predict(floresta_2d, grade)$predictions) |>
  pivot_longer(c(`Arvore unica`, `Floresta aleatoria`),
               names_to = "modelo", values_to = "classe")

ggplot(grade, aes(bill_length_mm, flipper_length_mm)) +
  geom_raster(aes(fill = classe), alpha = 0.35) +
  geom_point(data = dados, aes(color = species), size = 1) +
  scale_fill_manual(values = cores, guide = "none") +
  scale_color_manual(values = cores, name = "Especie") +
  facet_wrap(~ modelo) +
  labs(x = "Comprimento do bico (mm)", y = "Comprimento da nadadeira (mm)")
ggsave("figs/fig_fronteiras.pdf", width = 9, height = 4.4)

# ---- Comparacao final ----
tibble(modelo   = c("Arvore", "Floresta"),
       acuracia = c(mean(pred_arvore == teste$species),
                    mean(pred_floresta == teste$species)))
