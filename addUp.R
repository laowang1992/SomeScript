library(tidyverse)

if (FALSE) {
  setwd(dir = "E:/test/202209_QTLseq/03.Analysis/RY2P/")
  df <- read_tsv(file = "./RY2P.SlidingWindow.txt") %>% filter(str_detect(CHROM, "A"))
  len <- read_tsv(file = "../ref.len", col_names = c("CHROM", "Len")) %>% filter(str_detect(CHROM, "A"))
  group <- "CHROM"
  pos <- c("win_start", "win_end")
  band <- 0.01
}

addUp <- function(df, len = NULL, lenName = "Len", group, pos, band = 0.01){
  df_tmp <- df
  df <- df %>% ungroup() %>% dplyr::select(group = all_of(group), all_of(pos))
  # 取最大值
  if (is.null(len)) {
    len <- df %>% mutate(max = apply(df[,-1], 1, max)) %>% group_by(group) %>% summarise(Len = max(max))
  }else{
    len <- len %>% dplyr::select(group = all_of(group), Len = all_of(lenName))
  }
  accu <- rep(0, rep(nrow(len)))
  for (i in seq_along(accu)[-1]) {
    accu[i] <- sum(len$Len[1:i-1]) + sum(len$Len)*band*(i-1)
  }
  names(accu) <- len$group
  breaks <- accu + len$Len/2
  gaps <- accu[-1] - sum(len$Len)*band/2
  labels <- names(accu)
  
  for (p in pos) {
    dd <- df %>% select(group, p = all_of(p)) %>% mutate(p_addUp = 0)
    for (c in names(accu)) {
      dd <- dd %>% mutate(p_addUp = if_else(group == c, p+accu[c], p_addUp))
    }
    colnames(dd) <- c("group", p, paste(p, "addUp", sep = "_"))
    df <- df %>% left_join(dd, by = c("group", p))
  }
  
  colnames(df)[1] <- group
  df <- df %>% left_join(df_tmp, by = c(group, pos))
  outList <- list(df = df, breaks = breaks, labels = labels, gaps = gaps)
  return(outList)
}

# 测试函数
if (FALSE) {
  x <- addUp(df = df, len = len, group = group, pos = pos, band = 0)
  x$df$pos <- x$df$win_start_addUp/2+x$df$win_end_addUp/2
  ggplot(x$df %>% filter(nSNPs > 20), aes(x = pos, y = ED4, color = CHROM)) + 
    geom_line(linewidth = 1) +
    geom_vline(xintercept = x$gaps, linetype = "dashed", color = "gray") + 
    scale_x_continuous(breaks = x$breaks, labels = str_remove(x$labels, "scaffold"), expand = c(0, 0)) + 
    #scale_color_brewer(palette = "Set3") +
    ggsci::scale_color_npg() +
    labs(x = NULL) +
    cowplot::theme_half_open() + 
    theme(legend.position = "NA")
}
