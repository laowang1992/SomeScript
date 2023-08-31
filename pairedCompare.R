# paired comparison
# too simple, nothing interesting
library(tidyverse)
library(gghalves)
library(cowplot)
library(ggsignif)

pairedCompare <- function(df, id, comparison, color, xlab = waiver(), ylab = waiver()){
  dd <- df %>% select(ID = all_of(id), all_of(comparison)) %>% gather(x, y, -ID)
  dd$x <- factor(dd$x, levels = comparison)
  p <- ggplot() + 
    geom_half_violin(data = dd %>% filter(x == comparison[1]), mapping = aes(x = x, y = y), fill = color[1], side = "l") + 
    geom_half_violin(data = dd %>% filter(x == comparison[2]), mapping = aes(x = x, y = y), fill = color[2], side = "r") +
    geom_boxplot(data = dd, mapping = aes(x = x, y = y), width = 0.15) + 
    geom_line(data = dd, mapping = aes(x = x, y = y, group = ID), color = "gray") + 
    geom_signif(data = dd, mapping = aes(x = x, y = y), 
                comparisons = list(comparison), map_signif_level = F, 
                test = t.test, test.args = list(paired = TRUE)) +
    labs(x = xlab, y = ylab) +
    theme_half_open()
  return(p)
}

if (FALSE) {
  setwd("G:/其他/Ali/202205_Phenotype_plot")
  xlab <- NULL
  ylab <- "Ovlue Number"
  id = "ID"
  comparison <- c("DH-EZ16-2", "DH-EZ17")
  color = c("blue", "red")
  df <- read_tsv(file = "./ON_DH_data.txt") 
  p <- pairedCompare(df = df, id = id, comparison = comparison, color = color, xlab = xlab, ylab = ylab)
  p
  ggsave(p, filename = "comparision.pdf", width = 4, height = 4)
}

