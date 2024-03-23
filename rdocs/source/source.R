# Pacotes ----

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse, data.table,
  readxl, readr, ggcorrplot, cowplot,
  RColorBrewer, scales, nortest, xlsx
  )

# Fonte ----

windowsFonts(Arial=windowsFont("sans"))

# Cores ----

# Definindo paleta de cores
vetor_cores <- c("#ffffff","#000000","#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF","#FF00FF","#C0C0C0","#808080","#800000","#808000","#008000","#800080","#008080","#000080")


# Funções ----

theme_proprio <- function(...) {
  theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(colour = "#000000", size = 10),
      axis.title.x = ggplot2::element_text(colour = "#000000", size = 10),
      axis.text = ggplot2::element_text(colour = "#000000", size = 10),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "#000000"),
      text = element_text(family = "sans", size = 12),
      legend.position = "top",
      ...
    )
  return(
    list(
      theme,
      scale_fill_manual(values = vetor_cores),
      scale_colour_manual(values = vetor_cores),
      scale_y_continuous(
        labels = scales::number_format(decimal.mark = ',',
#                                       accuracy = 0.01,
                                       big.mark = "."))
    )
  )
}

# Definindo função que retorna frequências relativas de um vetor
percent <- function(absolute, digits = 2) {
  return(round(100 * absolute / sum(absolute), digits))
}

# Definindo função que retorna banco de dados com frequências
# relativas e absolutas de uma variável categórica
vector_frequencies <- function(vector) {
  frequency <- vector %>%
    table() %>%
    as_tibble() %>%
    mutate(
      rel = n %>%
        percent() %>%
        paste("%", sep = "")
    )
  colnames(frequency) <- c("groups", "absolute", "relative")
  return(frequency)
}

'%notin%' = Negate('%in%')