# =============================================================
# Arvores de decisao, bagging e florestas aleatorias
# Aplicacao: classificacao das especies de iris
# Referencia: Izbicki & dos Santos, Secoes 8.3 e 8.4
# =============================================================

library(tidyverse)
library(rpart)
library(rpart.plot)
library(ranger)

dir.create("figs_iris", showWarnings = FALSE)

# -------------------------------------------------------------
# 1. Dados e divisao treino/teste
# -------------------------------------------------------------
dados <- as_tibble(iris)

set.seed(2026)
id_treino <- sample(nrow(dados), 0.7 * nrow(dados))
treino <- dados %>% slice(id_treino)
teste  <- dados %>% slice(-id_treino)

# -------------------------------------------------------------
# 2. Analise exploratoria
# -------------------------------------------------------------
g_eda <- ggplot(dados, aes(x = Petal.Length, y = Petal.Width,
                           color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "Comprimento da petala (cm)", y = "Largura da petala (cm)",
       color = "Especie", title = "Iris: as petalas ja separam bem as especies") +
  theme_minimal(base_size = 13)

print(g_eda)
ggsave("figs_iris/eda.png", g_eda, width = 7, height = 4.5, dpi = 150)

# -------------------------------------------------------------
# 3. Exemplo 1: arvore de classificacao (CART)
# -------------------------------------------------------------
# cresce uma arvore grande e poda pelo erro de validacao cruzada
arvore_cheia <- rpart(Species ~ ., data = treino, method = "class", cp = 0)
cp_otimo <- arvore_cheia$cptable %>%
  as_tibble() %>%
  slice_min(xerror, n = 1, with_ties = FALSE) %>%
  pull(CP)
arvore <- prune(arvore_cheia, cp = cp_otimo)

rpart.plot(arvore, type = 2, extra = 104, box.palette = "Blues")
pdf("figs_iris/arvore.pdf", width = 8, height = 5)
rpart.plot(arvore, type = 2, extra = 104, box.palette = "Blues")
dev.off()

pred_arvore <- predict(arvore, teste, type = "class")
acc_arvore <- mean(pred_arvore == teste$Species)

# -------------------------------------------------------------
# 4. Exemplo 2: bagging (floresta com mtry = p)
# -------------------------------------------------------------
p <- ncol(treino) - 1   # numero de preditores

bagging <- ranger(Species ~ ., data = treino,
                  num.trees = 500, mtry = p, seed = 2026)
acc_bagging <- mean(predict(bagging, teste)$predictions == teste$Species)

# -------------------------------------------------------------
# 5. Exemplo 3: floresta aleatoria (mtry = sqrt(p), padrao)
# -------------------------------------------------------------
floresta <- ranger(Species ~ ., data = treino,
                   num.trees = 500, importance = "impurity", seed = 2026)
acc_floresta <- mean(predict(floresta, teste)$predictions == teste$Species)

# -------------------------------------------------------------
# 6. Comparacao dos tres modelos
# -------------------------------------------------------------
comparacao <- tibble(
  modelo    = c("Arvore", "Bagging (mtry = p)", "Floresta (mtry = sqrt(p))"),
  acc_teste = c(acc_arvore, acc_bagging, acc_floresta),
  erro_oob  = c(NA, bagging$prediction.error, floresta$prediction.error)
)
print(comparacao)

g_comp <- comparacao %>%
  mutate(modelo = fct_inorder(modelo)) %>%
  ggplot(aes(x = modelo, y = acc_teste)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_text(aes(label = scales::percent(acc_teste, accuracy = 0.1)),
            vjust = -0.4, size = 4.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  labs(x = NULL, y = "Acuracia no teste",
       title = "Arvore vs. bagging vs. floresta aleatoria") +
  theme_minimal(base_size = 13)

print(g_comp)
ggsave("figs_iris/comparacao.png", g_comp, width = 7, height = 4.5, dpi = 150)

# -------------------------------------------------------------
# 7. Importancia de variaveis (floresta)
# -------------------------------------------------------------
g_imp <- floresta$variable.importance %>%
  enframe(name = "variavel", value = "importancia") %>%
  ggplot(aes(x = importancia, y = fct_reorder(variavel, importancia))) +
  geom_col(fill = "steelblue") +
  labs(x = "Importancia (reducao de Gini)", y = NULL,
       title = "Quais variaveis mais ajudam a classificar?") +
  theme_minimal(base_size = 13)

print(g_imp)
ggsave("figs_iris/importancia.png", g_imp, width = 7, height = 4.5, dpi = 150)

# -------------------------------------------------------------
# 8. Fronteiras de decisao: arvore vs. floresta (2 variaveis)
# -------------------------------------------------------------
arvore_2d   <- rpart(Species ~ Petal.Length + Petal.Width,
                     data = treino, method = "class")
floresta_2d <- ranger(Species ~ Petal.Length + Petal.Width,
                      data = treino, num.trees = 500, seed = 2026)

grade <- expand_grid(
  Petal.Length = seq(min(dados$Petal.Length), max(dados$Petal.Length), length.out = 200),
  Petal.Width  = seq(min(dados$Petal.Width),  max(dados$Petal.Width),  length.out = 200)
)

fronteiras <- bind_rows(
  grade %>% mutate(modelo = "Arvore",
                  pred = predict(arvore_2d, grade, type = "class")),
  grade %>% mutate(modelo = "Floresta aleatoria",
                  pred = predict(floresta_2d, grade)$predictions)
)

g_front <- ggplot(fronteiras, aes(Petal.Length, Petal.Width)) +
  geom_raster(aes(fill = pred), alpha = 0.35) +
  geom_point(data = dados, aes(color = Species), size = 1.2) +
  facet_wrap(~ modelo) +
  labs(x = "Comprimento da petala (cm)", y = "Largura da petala (cm)",
       fill = "Previsao", color = "Especie",
       title = "Fronteiras de decisao: cortes retos vs. media de muitas arvores") +
  theme_minimal(base_size = 13)

print(g_front)
ggsave("figs_iris/fronteiras.png", g_front, width = 9, height = 4.5, dpi = 150)