source("rdocs/source/source.R")

# -------------------------------- Gráficos --------------------------------- #

# 1. Barras/Colunas ----
# 1.1 Colunas com duas frequências ----
# 1.1.1 Univariado ----
classes <- mpg |>
  filter(!is.na(class)) |>
  count(class) |>
  mutate(
    freq = n |> percent(),
  ) |>
  mutate(
    freq = gsub("\\.", ",", freq) |> paste("%", sep = ""),
    label = formatC(n,big.mark = ".",decimal.mark = ",", format = "d"),
    label = str_c(label, " (", freq, ")") |> str_squish()
  )

ggplot(classes) +
  aes(
    x = fct_reorder(class, n, .desc = T),
    y = n,
    label = label
  ) +
  geom_bar(stat = "identity", fill = "#006eab", width = 0.7) +
  geom_text(
    position = position_dodge(width = .9),
    vjust = -0.5, # hjust = .5,
    size = 3
  ) +
  labs(x = "", y = "") +
  theme_minimal()

#ggsave("colunas-uni-freq.pdf", width = 158, height = 93, units = "mm")

# 1.1.2 Bivariado ----
trans_drv <- mpg |>
  mutate(trans = case_when(
    trans |> str_detect("auto") ~ "auto",
    trans |> str_detect("manual") ~ "manual"
  )) |>
  group_by(trans, drv) |>
  summarise(freq = n()) |>
  mutate(
    freq_relativa = freq |> percent()
  )

porcentagens <- str_c(trans_drv$freq_relativa, "%") |> str_replace("\\.", ",")

legendas <- str_squish(str_c(trans_drv$freq, " (", porcentagens, ")"))

ggplot(trans_drv) +
  aes(
    x = fct_reorder(trans, freq, .desc = T),
    y = freq,
    fill = drv,
    label = legendas
  ) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0)) +
  geom_text(
    position = position_dodge(width = .9),
    vjust = -0.5, hjust = 0.5,
    size = 3
  ) +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("colunas-bi-freq.pdf", width = 158, height = 93, units = "mm")

# 1.2 Barras com duas frequências ----
# 1.2.1 Univariado ----

ggplot(classes) +
  aes(
    x = fct_reorder(class, n, .desc = T),
    y = n,
    label = label
  ) +
  geom_bar(stat = "identity", fill = "#006eab", width = 0.7) +
  geom_text(
    position = position_dodge(width = .9),
    vjust = 0, hjust = -.1,
    size = 3
  ) +
  scale_y_continuous(breaks = seq(0, 70, 20), limits = c(0, 70)) +
  labs(x = "", y = "") +
  theme_minimal() +
  coord_flip()
#ggsave("barras-uni-freq.pdf", width = 158, height = 93, units = "mm")

# 1.2.2 Bivariado ----
class_drv <- mpg |>
  group_by(class, drv) |>
  summarise(freq = n()) |>
  mutate(freq_relativa = freq |> percent())

porcentagens <- str_c(class_drv$freq_relativa, "%") |> str_replace("\\.", ",")
legendas <- str_squish(str_c(class_drv$freq, " (", porcentagens, ")"))

ggplot(class_drv) +
  aes(
    x = fct_reorder(class, freq, .desc = T),
    y = freq,
    fill = drv,
    label = legendas
  ) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0)) +
  geom_text(
    position = position_dodge(width = .9),
    vjust = 0.2, hjust = -0.1,
    size = 3
  ) +
  labs(x = "", y = "") +
  theme_minimal() +
  scale_y_continuous(name = "Speed of cars", limits = c(0, 60)) +
  coord_flip()

#ggsave("barras-bi-freq.pdf", width = 158, height = 93, units = "mm")

# 1.3 Colunas ----
# 1.3.1 Univariado ----
ggplot(mpg) +
  aes(x = class) +
  geom_bar(fill = "#006eab") +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("colunas-uni-freq.pdf", width = 158, height = 93, units = "mm")

# 1.3.2 Univariado com porcentagem no gráfico e no eixo
ggplot(mpg) +
  aes(x = class) +
  geom_bar(aes(y = prop.table(after_stat(count)) * 100), fill = "#006eab") +
  geom_text(aes(
    y = prop.table(after_stat(count)) * 100 + 0.5,
    label = paste0(gsub("\\.", ",", round(prop.table(..count..) * 100, 2)), "%")
  ),
  stat = "count", vjust = 0, size = 4
  ) +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("colunas-uni-percent.pdf", width = 158, height = 93, units = "mm")

# 1.3.3 Univariado com porcentagem no gráfico e freq absoluta no eixo ----
ggplot(mpg$class |> vector_frequencies()) +
  aes(
    x = groups,
    y = absolute,
    label = relative
  ) +
  geom_bar(stat = "identity", fill = "#006eab") +
  geom_text(vjust = -0.5, size = 4) +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("colunas-uni-freq-percent.pdf", width = 158, height = 93, units = "mm")

# 1.3.4 Bivariado com dodge ----
class_trans <- as.data.frame(table(mpg$class, mpg$trans))
ggplot(class_trans) +
  aes(x = Var1, y = Freq, fill = Var2) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(name = "Transmissão") +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("colunas-bi-dodge.pdf", width = 158, height = 93, units = "mm")

# 1.3.5 Bivariado com stack ----
ggplot(class_trans, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(name = "Transmissão") +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("colunas-bi-stack.pdf", width = 158, height = 93, units = "mm")


# 1.3.6 Bivariado com fill ----
ggplot(class_trans, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(name = "Transmissão") +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("colunas-bi-fill.pdf", width = 158, height = 93, units = "mm")

# 1.3.7 Bivariado com porcentagem ----

trans_class <- table(mpg$trans, mpg$class) |>
  data.frame() |>
  mutate(Pct = Freq / sum(Freq))

orderclass <- c(
  "2seater", "compact", "midsize",
  "minivan", "pickup", "subcompact",
  "suv")

ggplot(trans_class) +
  aes(
    x = factor(Var2, level = orderclass),
    y = Pct,
    fill = Var1
  ) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(name = "Transmissão") +
  scale_y_continuous(
    limits = c(0, .15), expand = c(0, 0), breaks = seq(0, .3, .05),
    labels = paste0(seq(0, 30, 5), "%")
  ) +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("colunas-bivariado-percent.pdf", width = 158, height = 93, units = "mm")

# 1.4 Barras ----
# Basta adicionar coord_flip() nos códigos para Colunas

# 2. Setores ----
# 2.1 Com porcentagem ----
contagem <- mpg |>
  group_by(drv) |>
  summarise(Freq = n()) |>
  mutate(Prop = round(100 * (Freq / sum(Freq)), 2)) |>
  arrange(desc(drv)) |>
  mutate(posicao = cumsum(Prop) - 0.5 * Prop,
         ymax = cumsum(Prop),
         ymin = c(0, head(ymax, n=-1)),
         Freq = formatC(Freq,big.mark = ".",decimal.mark = ",", format = "d"),
         Prop_L = formatC(Prop,big.mark = ".",decimal.mark = ",")
         )

ggplot(contagem) +
  aes(
    x = factor(""),
    y = Prop,
    fill = factor(drv)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_x_discrete() +
  scale_fill_manual(name = "DRV") +
  geom_text(
    aes(x = 1.8, y = posicao, label = paste0(Freq," ","(", Prop, "%)")),
    color = "#000000") + 
  theme_minimal()
#ggsave("setor.pdf", width = 158, height = 93, units = "mm")

# Gráfico de roskinha kkk ----

ggplot(contagem, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=drv)) +
  geom_rect() +
  coord_polar(theta="y") + 
  xlim(c(2, 4.5)) + 
  theme_minimal() +
  theme_void() +
  geom_text(
    aes(x = 4.4, y = posicao, label = paste0(Freq," ","(", Prop, "%)")),
    color = "#000000") +
  scale_fill_manual(values = vetor_cores, name = 'Classe')

# 3. Boxplot ----
# 3.1 Univariado ----
ggplot(mpg) +
  aes(
    x = factor(""),
    y = cty
  ) +
  geom_boxplot(fill = c("#006eab"), width = 0.5) +
  stat_summary(
    fun = "mean", geom = "point", shape = 23, size = 3, fill = "#ffffff"
  ) +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("box_uni.pdf", width = 158, height = 93, units = "mm")

# 3.2 Bivariado ----
ggplot(mpg) +
  aes(
    x = trans,
    y = cty
  ) +
  geom_boxplot(fill = c("#006eab"), width = 0.5) +
  stat_summary(
    fun = "mean", geom = "point", shape = 23, size = 3, fill = "#ffffff"
  ) +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("box_bi.pdf", width = 158, height = 93, units = "mm")




# 4. Histograma ----
# 4.1 Univariado ----
ggplot(mpg) +
  aes(x = cty) +
  geom_histogram(colour = "#000000", fill = "#006eab", binwidth = 7) +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("hist_uni.pdf", width = 158, height = 93, units = "mm")

# 4.2 Univariado em porcentagem ----
ggplot(mpg) +
  aes(x = cty) +
  geom_histogram(
    aes(y = 100 * (..count..) / sum(..count..)),
    colour = "#000000",
    fill = "#006eab",
    binwidth = 7
  ) +
  labs(x = "", y = "") +
  theme_minimal()
#ggsave("hist_uni_porc.pdf", width = 158, height = 93, units = "mm")

# 4.3 Bivariado com facet grid ----
ggplot(mpg) +
  aes(x = cty) +
  geom_histogram(colour = "#000000", fill = "#006eab", binwidth = 7) +
  facet_grid(. ~ class) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    strip.background = element_rect(colour = "#000000", fill = "#ffffff")
  )
#ggsave("hist_grid.pdf", width = 200, height = 93, units = "mm")

# 5. Dispersão ----
# 5.1 Univariado com poucos pontos sobrepostos ----
ggplot(mpg, aes(x = cty, y = hwy)) +
  geom_point(colour = "#006eab", size = 3) +
  labs(x = "",y = "") +
  theme_minimal()

# 5.2 Univariado com muitos pontos sobrepostos ----
# 5.2.1 geom_jitter ----
ggplot(mpg, aes(x = cyl, y = cty)) +
  geom_jitter(colour = "#006eab", size = 3) +
  labs(x = "",y = "") +
  theme_minimal()

# 5.2.2 Alpha ----
ggplot(mpg, aes(x = cyl, y = cty)) +
  geom_point(
    colour = "#006eab",
    size = 3,
    alpha = 0.3
  ) +
  labs(x = "",y = "") +
  theme_minimal()

# 5.3 Bivariado ----
mpg$trans <- factor(substr(mpg$trans, 1, nchar(as.character(mpg$trans)) - 4))

ggplot(mpg, aes(x = cty, y = hwy)) +
  geom_point(aes(colour = trans)) +
  scale_colour_manual(
    name = "Transmissão", values = c("#006eab", "#00a8e7"),
    labels = c("Automático", "Manual")
  ) +
  labs(
    x = "",
    y = ""
  ) +
  theme_minimal()

# 6. Linhas ----
dados <- tibble(
  ano = c(
    "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014",
    "2015", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013",
    "2014", "2015"
  ),
  preco = c(
    2.5, 5.7, 3.4, 3.7, 4.5, 4.8, 4.1, 4.6, 4.8,
    5, 4.5, 6.5, 3.5, 4.6, 4.7, 4.9, 5, 5.5, 3.5,
    7.5
  ),
  produto = c(
    "a", "a", "a", "a", "a", "a", "a", "a", "b",
    "b", "b", "b", "b", "b", "b", "b", "b", "b",
    "b", "b"
  )
)

# 6.1 Univariado ----
ggplot(dados) +
  aes(x = ano, y = preco, group = 1) +
  geom_line(linewidth = 1, colour = "#00a8e7") +
  geom_point(colour = "#006eab", size = 2) +
  labs(x = "", y = "") +
  theme_minimal()

# 6.2 Bivariado ----
ggplot(dados) +
  aes(x = ano, y = preco, group = produto, colour = produto) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(name = "Produto", labels = c("A", "B")) +
  labs(x = "", y = "") +
  theme_minimal()

# 7.0 Diagrama de Sankey ----

prop <- mpg |>
  select(trans,class) |>
  count(trans, class) |>
  mutate(proptot = prop.table(n))

# Caso seja necessário ordenar os fatores, utilize a linha abaixo; além de remover
# os # do código ggplot abaixo para o filtro funcionar.

# fct_lvl <- c("Muito","Médio","Pouco","Nada")

library("ggalluvial")

ggplot(as.data.frame(prop),
       aes(y = proptot, axis1 = factor(trans
                                       #,level=fct_lvl
       ), axis2 = factor(class
                         #,level=fct_lvl
       ))) +
  geom_alluvium(aes(fill = factor(trans
                                  #,level=fct_lvl
  )), width = 1/12,alpha=.8,show.legend = FALSE) +
  geom_stratum(width = 1/12, fill = "#006eab", colour = "#000000",alpha=1) +
  geom_label(stat = "stratum") +
  aes(label = after_stat(stratum)) +
  scale_x_discrete(limits = c("Trans", "Class"),
                   expand = c(.05, .05),
                   labels = c("Trans", "Class")) +
  scale_fill_manual(values = vetor_cores) +
  scale_y_continuous(labels = NULL,
                     name = NULL,
                     breaks = NULL) +
  theme_minimal()
#ggsave("diagrama_de_sankey.pdf", width = 158, height = 93, units = "mm") 

# 9.0 Treemap ----

library("ggfittext")
library("treemapify")

cars <- mtcars
cars$carname <- rownames(cars)
cars <- mutate(cars, cyl = factor(cyl))

ggplot(cars, aes(area = disp, fill = cyl, label = carname)) +
  geom_treemap() +
  geom_treemap_text(
    fontface = "italic",
    colour = "#000000",
    place = "centre",
    grow = TRUE
  ) + 
  theme_minimal()

# 10.0 Clusters ----
# 10.1 K-means ----
# Gráficos necessários para visualização de clusterização k-means

data("mtcars")
df <- scale(mtcars)

library("factoextra")

# Graficar o número ótimo de clusters pelo método de elbow ----
fviz_nbclust(df, kmeans, 
             method = "wss", # soma de quadrados totais
             k.max = 10, # máximo de clusters
             nboot = 1000, # Máximo de bootstraps
             barfill = "#006eab",
             barcolor = "#006eab",
             linecolor = "#00a8e7") +
  geom_vline(xintercept = 2, linetype = 2) +
  theme_minimal()

set.seed(9803)
km.res <- kmeans(df, 2, nstart=25)

#aggregate(mtcars, by=list(cluster=km.res$cluster), mean)
mtcars2 <- cbind(mtcars, cluster=km.res$cluster)

# Visualizando os clusters
fviz_cluster(km.res, data=mtcars2,
             ellipse.type="euclid",
             star.plot=TRUE,
             repel=TRUE,
             ggtheme=theme_minimal())
