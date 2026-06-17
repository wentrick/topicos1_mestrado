# =============================================================
# Arvores de decisao, bagging e florestas aleatorias
# Aplicacao: Sonar (mlbench) - classificar eco de sonar como
#            Rocha (R) ou Mina metalica (M)
# 60 preditores numericos, 208 observacoes: caso DIFICIL,
# onde uma arvore unica vai mal e a floresta da um salto.
# Referencia: Izbicki & dos Santos, Secoes 8.3 e 8.4
# =============================================================

library(tidyverse)
library(rpart)
library(rpart.plot)
library(ranger)

# mlbench traz o dataset Sonar; instala se necessario
if (!requireNamespace("mlbench", quietly = TRUE)) install.packages("mlbench")

dir.create("figs_sonar", showWarnings = FALSE)

# -------------------------------------------------------------
# 1. Dados e divisao treino/teste
# -------------------------------------------------------------
# package = "mlbench" garante que o Sonar e encontrado mesmo que
# o pacote ainda nao tenha sido carregado com library()
data("Sonar", package = "mlbench")
dados <- as_tibble(Sonar)   # colunas V1..V60 (numericas) e Class (R/M)

set.seed(2026)
id_treino <- sample(nrow(dados), 0.7 * nrow(dados))
treino <- dados %>% slice(id_treino)
teste  <- dados %>% slice(-id_treino)

p <- ncol(treino) - 1   # 60 preditores

# -------------------------------------------------------------
# 2. Exemplo 1: arvore de classificacao (CART) com poda por CV
# -------------------------------------------------------------
arvore_cheia <- rpart(Class ~ ., data = treino, method = "class", cp = 0)
cp_otimo <- arvore_cheia$cptable %>%
  as_tibble() %>%
  slice_min(xerror, n = 1, with_ties = FALSE) %>%
  pull(CP)
arvore <- prune(arvore_cheia, cp = cp_otimo)

rpart.plot(arvore, type = 2, extra = 104, box.palette = "Oranges")
pdf("figs_sonar/arvore.pdf", width = 8, height = 5)
rpart.plot(arvore, type = 2, extra = 104, box.palette = "Oranges")
dev.off()

acc_arvore <- mean(predict(arvore, teste, type = "class") == teste$Class)

# -------------------------------------------------------------
# 3. Exemplo 2: bagging (floresta com mtry = p)
# -------------------------------------------------------------
bagging <- ranger(Class ~ ., data = treino,
                  num.trees = 500, mtry = p, seed = 2026)
acc_bagging <- mean(predict(bagging, teste)$predictions == teste$Class)

# -------------------------------------------------------------
# 4. Exemplo 3: floresta aleatoria (mtry = sqrt(p), padrao)
# -------------------------------------------------------------
floresta <- ranger(Class ~ ., data = treino,
                   num.trees = 500, importance = "impurity", seed = 2026)
acc_floresta <- mean(predict(floresta, teste)$predictions == teste$Class)

# -------------------------------------------------------------
# 5. Comparacao dos tres modelos
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
  geom_col(fill = "darkorange", width = 0.6) +
  geom_text(aes(label = scales::percent(acc_teste, accuracy = 0.1)),
            vjust = -0.4, size = 4.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  labs(x = NULL, y = "Acuracia no teste",
       title = "Sonar: a floresta supera a arvore unica com folga") +
  theme_minimal(base_size = 13)

print(g_comp)
ggsave("figs_sonar/comparacao.png", g_comp, width = 7, height = 4.5, dpi = 150)

# -------------------------------------------------------------
# 6. Curva: erro OOB cai conforme adicionamos arvores
#    (mostra quando parar de adicionar arvores compensa)
# -------------------------------------------------------------
n_trees_grid <- c(1, 2, 5, 10, 25, 50, 100, 200, 350, 500)
curva_oob <- tibble(num_trees = n_trees_grid) %>%
  mutate(erro_oob = map_dbl(num_trees, \(B)
                            ranger(Class ~ ., data = treino, num.trees = B, seed = 2026)$prediction.error))

g_oob <- ggplot(curva_oob, aes(num_trees, erro_oob)) +
  geom_line(color = "darkorange", linewidth = 1) +
  geom_point(color = "darkorange", size = 2) +
  labs(x = "Numero de arvores (B)", y = "Erro OOB",
       title = "O erro OOB estabiliza: depois de ~100 arvores, pouco muda") +
  theme_minimal(base_size = 13)

print(g_oob)
ggsave("figs_sonar/curva_oob.png", g_oob, width = 7, height = 4.5, dpi = 150)

# -------------------------------------------------------------
# 7. Curva: efeito do mtry (liga bagging e floresta empiricamente)
#    mtry pequeno = floresta; mtry = p = bagging
# -------------------------------------------------------------
mtry_grid <- c(1, 2, 4, 8, 12, 20, 30, 45, 60)
curva_mtry <- tibble(mtry = mtry_grid) %>%
  mutate(erro_oob = map_dbl(mtry, \(m)
                            ranger(Class ~ ., data = treino, num.trees = 500, mtry = m,
                                   seed = 2026)$prediction.error))

g_mtry <- ggplot(curva_mtry, aes(mtry, erro_oob)) +
  geom_line(color = "darkorange", linewidth = 1) +
  geom_point(color = "darkorange", size = 2) +
  geom_vline(xintercept = floor(sqrt(p)), linetype = "dashed") +
  geom_vline(xintercept = p, linetype = "dotted") +
  annotate("text", x = floor(sqrt(p)), y = max(curva_mtry$erro_oob),
           label = "floresta (sqrt p)", hjust = -0.1, size = 3.5) +
  annotate("text", x = p, y = max(curva_mtry$erro_oob),
           label = "bagging (p)", hjust = 1.1, size = 3.5) +
  labs(x = "mtry (variaveis candidatas por corte)", y = "Erro OOB",
       title = "Por que a floresta vence: mtry intermediario minimiza o erro") +
  theme_minimal(base_size = 13)

print(g_mtry)
ggsave("figs_sonar/curva_mtry.png", g_mtry, width = 7, height = 4.5, dpi = 150)

# -------------------------------------------------------------
# 8. Importancia de variaveis (top 15 das 60)
# -------------------------------------------------------------
g_imp <- floresta$variable.importance %>%
  enframe(name = "variavel", value = "importancia") %>%
  slice_max(importancia, n = 15) %>%
  ggplot(aes(x = importancia, y = fct_reorder(variavel, importancia))) +
  geom_col(fill = "darkorange") +
  labs(x = "Importancia (reducao de Gini)", y = NULL,
       title = "As 15 variaveis mais informativas (de 60)") +
  theme_minimal(base_size = 13)

print(g_imp)
ggsave("figs_sonar/importancia.png", g_imp, width = 7, height = 4.5, dpi = 150)

# -------------------------------------------------------------
# 9. Fronteiras de decisao nas 2 variaveis mais importantes
# -------------------------------------------------------------
top2 <- floresta$variable.importance %>%
  enframe(name = "variavel", value = "importancia") %>%
  slice_max(importancia, n = 2) %>%
  pull(variavel)

f <- as.formula(paste("Class ~", paste(top2, collapse = " + ")))
arvore_2d   <- rpart(f, data = treino, method = "class")
floresta_2d <- ranger(f, data = treino, num.trees = 500, seed = 2026)

grade <- expand_grid(
  v1 = seq(min(dados[[top2[1]]]), max(dados[[top2[1]]]), length.out = 200),
  v2 = seq(min(dados[[top2[2]]]), max(dados[[top2[2]]]), length.out = 200)
) %>%
  set_names(top2)

fronteiras <- bind_rows(
  grade %>% mutate(modelo = "Arvore",
                  pred = predict(arvore_2d, grade, type = "class")),
  grade %>% mutate(modelo = "Floresta aleatoria",
                  pred = predict(floresta_2d, grade)$predictions)
)

g_front <- ggplot(fronteiras, aes(.data[[top2[1]]], .data[[top2[2]]])) +
  geom_raster(aes(fill = pred), alpha = 0.35) +
  geom_point(data = dados, aes(color = Class), size = 1.2) +
  facet_wrap(~ modelo) +
  labs(x = top2[1], y = top2[2], fill = "Previsao", color = "Classe",
       title = "Fronteiras nas 2 variaveis mais importantes") +
  theme_minimal(base_size = 13)

print(g_front)
ggsave("figs_sonar/fronteiras.png", g_front, width = 9, height = 4.5, dpi = 150)