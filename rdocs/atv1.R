consumo = factor(c("0","<1","1-2","3-5",">=6"),
                 levels = c("0","<1","1-2","3-5",">=6"),
                 ordered = TRUE)
ausente = c(17066,14464,788,126,37)
presente = c(48,38,5,1,1)
df1 = data.frame(consumo, ausente, presente)
pacman::p_load(tidyverse,epitools)

df = df1 %>%
  pivot_longer(cols = c(ausente, presente), names_to = "estado", values_to = "n") %>%
  uncount(weights = n)
df$estado = factor(df$estado)

outcome <- c('Ausente', 'Presente')
data <- matrix(c(17066, 48, 14464, 38,788,5,126,1,37,1), nrow=5, ncol=2, byrow=TRUE)
dimnames(data) <- list('Consumo de Álcool'=consumo, 'Presença de malformação'=outcome)
or = oddsratio(data)
or

chisq.test(data) # chi-squared contingency table tests and goodness-of-fit tests.

fisher.test(data) # testing the null of independence of rows and columns in a contingency table with fixed marginals.

# ---------------------------------------------------------------------------- #
# H_0) Não existe associação entre consumo de álcool e presença de malformação.
# H_1) c.c.

df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 0,
      consumo == "<1" ~ .5,
      consumo == "1-2" ~ 1.5,
      consumo == "3-5" ~ 4,
      consumo == ">=6" ~ 7),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value
# Com este p-valor, rejeitamos a hipótese nula, ou seja, aparenta existir associação entre consumo de álcool e presença de malformação.

# Para algum outro valor arbitrário de score:
# 1 a 5 ----
df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 1,
      consumo == "<1" ~ 2,
      consumo == "1-2" ~ 3,
      consumo == "3-5" ~ 4,
      consumo == ">=6" ~ 5),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
r
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value

# 0 a 4 ----
df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 0,
      consumo == "<1" ~ 1,
      consumo == "1-2" ~ 2,
      consumo == "3-5" ~ 3,
      consumo == ">=6" ~ 4),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
r
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value

# outro ----
df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 0,
      consumo == "<1" ~ 10,
      consumo == "1-2" ~ 20,
      consumo == "3-5" ~ 30,
      consumo == ">=6" ~ 40),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
r
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value
# Daqui, vemos que se o espaçamento entre os scores é igual, o p-valor é o mesmo, então não importa o valor do score, mas sim a diferença entre eles.

# Solicitados na atividade:

# 1,2,3,4 e 5 ----

df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 1,
      consumo == "<1" ~ 2,
      consumo == "1-2" ~ 3,
      consumo == "3-5" ~ 4,
      consumo == ">=6" ~ 5),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
r
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value

# 2,4,6,8 e 10 ----

df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 2,
      consumo == "<1" ~ 4,
      consumo == "1-2" ~ 6,
      consumo == "3-5" ~ 8,
      consumo == ">=6" ~ 10),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
r
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value

# 10,20,30,40 e 50 ----

df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 2,
      consumo == "<1" ~ 4,
      consumo == "1-2" ~ 6,
      consumo == "3-5" ~ 8,
      consumo == ">=6" ~ 10),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
r
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value

# 1,2,4,7 e 12 ----

df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 1,
      consumo == "<1" ~ 2,
      consumo == "1-2" ~ 4,
      consumo == "3-5" ~ 7,
      consumo == ">=6" ~ 12),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
r
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value

# 5,4,3,2 e 1 ----

df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 5,
      consumo == "<1" ~ 4,
      consumo == "1-2" ~ 3,
      consumo == "3-5" ~ 2,
      consumo == ">=6" ~ 1),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
r
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value

# exponencial do ponto medio ----
df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 0,
      consumo == "<1" ~ exp(.5),
      consumo == "1-2" ~ exp(1.5),
      consumo == "3-5" ~ exp(4),
      consumo == ">=6" ~ exp(7)),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
r
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value

# ln do ponto medio ----
df_score = df |>
  mutate(
    consumo_score = case_when(
      consumo == "0" ~ 0,
      consumo == "<1" ~ log(.5),
      consumo == "1-2" ~ log(1.5),
      consumo == "3-5" ~ log(4),
      consumo == ">=6" ~ log(7)),
    estado = ifelse(estado == "ausente", 0, 1)) |>
  select(consumo_score, estado)

r = cor(df_score$consumo_score, df_score$estado, method = "pearson")
r
M2 = (nrow(df_score)-1)*(r^2)
M2
p_value = 1 - pchisq(M2, 1)
p_value

# --------------------------------------------------------------------------- #