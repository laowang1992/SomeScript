#!/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for background analysis")

## Add command line arguments
#
p <- add_argument(parser = p, arg = "--input", short = "-i", help = "input file name, GATK Table format", type = "character")

p <- add_argument(parser = p, arg = "--donor", short = "-d", help = "donor parent name", type = "character")
p <- add_argument(parser = p, arg = "--recurrent", short = "-r", help = "recurrent parent name", type = "character")
p <- add_argument(parser = p, arg = "--sample", short = "-s", help = "sample parent name", type = "character")

p <- add_argument(parser = p, arg = "--chromosome", short = "-c", help = "tab-separated file contain chromosome ID and label shown in figure. if not assigned, all chromosome in genome will be shown", type = "character")
p <- add_argument(parser = p, arg = "--length", short = "-l", help = "tab-separated file contain chromosome ID and length. if not assigned, the max position of SNP for each chromosome in GATK Table file will be regarded as chromosome length", type = "character")

p <- add_argument(parser = p, arg = "--minDdp", help = "Minimum depth for donor parent", type = "numeric", default = -Inf)
p <- add_argument(parser = p, arg = "--maxDdp", help = "Maxmum depth for donor parent", type = "numeric", default = Inf)
p <- add_argument(parser = p, arg = "--minRdp", help = "Minimum depth for recurrent parent", type = "numeric", default = -Inf)
p <- add_argument(parser = p, arg = "--maxRdp", help = "Maxmum depth for recurrent parent", type = "numeric", default = Inf)

p <- add_argument(parser = p, arg = "--path_to_snpbinner", short = "-p", help = "path to snpbinner", type = "character", default = "~/tools/SNPbinner-1.0.0-GondaEtAl2019/snpbinner")

# Parse the command line arguments
argv <- parse_args(p)

input <- argv$input
donor <- argv$donor
recurrent <- argv$recurrent
sample <- argv$sample
chromosome <- argv$chromosome
length <- argv$length
minDdp <- argv$minDdp
maxDdp <- argv$maxDdp
minRdp <- argv$minRdp
maxRdp <- argv$maxRdp
path_to_snpbinner <- argv$path_to_snpbinner

library(tidyverse)
library(cowplot)
library(ggh4x)

if (0) {
  input <- "SYX_bg.filter.SNPs.table.gz"
  donor <- "9172"
  recurrent <- "Y633"
  sample <- "WS100"
  chromosome <- "./chrom.txt"
  length <- "./ref.len"
  #winSize <- 15
  minDdp <- 3
  maxDdp <- 50
  minRdp <- 3
  maxRdp <- 50
  
  #fixSize <- 10
  
  path_to_snpbinner <- "~/tools/SNPbinner-1.0.0-GondaEtAl2019/snpbinner"
}

# 读取数据
#source("./Support_functions.R")
# 创建输出文件夹
dir.create(path = sample, showWarnings = F)

df <- read_tsv(file = input)
if (is.na(chromosome)) {
  chr <- df %>% distinct(CHROM) %>% mutate(LABEL = CHROM, y = rev(seq_along(rownames(.))*3))
} else {
  chr <- read_tsv("./chrom.txt", col_names = c("CHROM", "LABEL"), show_col_types = FALSE) %>% mutate(y = rev(seq_along(rownames(.))*3))
}
if (is.na(length)) {
  len <- df %>% group_by(CHROM) %>% summarise(Len = max(POS)) %>% right_join(chr, by = "CHROM")
} else {
  len <- read_tsv(file = "./ref.len", col_names = c("CHROM", "Len"), show_col_types = FALSE) %>% right_join(chr, by = "CHROM")
}
genome_size <- sum(len$Len)

dd1 <- df %>% select(CHROM, POS, REF, ALT, 
                     donor.GT = all_of(paste(donor, "GT", sep = ".")), 
                     donor.DP = all_of(paste(donor, "DP", sep = ".")), 
                     donor.GQ = all_of(paste(donor, "GQ", sep = ".")), 
                     recurrent.GT = all_of(paste(recurrent, "GT", sep = ".")), 
                     recurrent.DP = all_of(paste(recurrent, "DP", sep = ".")), 
                     recurrent.GQ = all_of(paste(recurrent, "GQ", sep = ".")), 
                     sample.GT = all_of(paste(sample, "GT", sep = ".")), 
                     sample.DP = all_of(paste(sample, "DP", sep = ".")), 
                     sample.GQ = all_of(paste(sample, "GQ", sep = "."))) %>% 
  # gatk4结果中有的GT是”1|1“这种形式
  mutate(donor.GT     = str_replace(donor.GT,     "\\|", "/"), 
         recurrent.GT = str_replace(recurrent.GT, "\\|", "/"), 
         sample.GT    = str_replace(sample.GT,    "\\|", "/")) %>% 
  # 有的结果中GT是”./.“，DP和GQ仍然有数值而不是NA，针对这种情况做过滤
  filter(donor.GT != "./." & recurrent.GT != "./." & sample.GT != "./.")

# 统计DP分布
dp <- dd1 %>% dplyr::select(donor = donor.DP, recurrent = recurrent.DP, Sample = sample.DP) %>%
  gather(key = "sample", value = "depth")
dp$sample <- factor(dp$sample, levels = c("donor", "recurrent", "Sample"), labels = c(donor, recurrent, sample))
## 根据DP的ave+2*sd设定x轴上线
upper <- dp %>% group_by(sample) %>% summarise(ave = mean(depth), sd = sd(depth)) %>% mutate(upper = ave+2*sd) %>% pull(upper) %>% max()
P_dp <- ggplot(dp, aes(x = depth)) + 
  geom_histogram(aes(y = after_stat(density), fill = sample), binwidth = 2) +
  geom_density() + 
  scale_x_continuous(limits = c(0, upper)) + 
  theme_half_open() +
  facet_wrap(~sample, nrow = 1)
ggsave(P_dp, filename = paste(sample, "/", sample, ".depth_desity.pdf", sep = ""), height = 3.5, width = 8)
ggsave(P_dp, filename = paste(sample, "/", sample, ".depth_desity.png", sep = ""), height = 3.5, width = 8, dpi = 500)
remove(P_dp)

# 过滤亲本DP、GQ
dd2 <- dd1 %>% filter(donor.DP > minDdp, donor.DP < maxDdp, recurrent.DP > minRdp, recurrent.DP < maxRdp) 
# 筛选亲本为aa x bb模式的位点，并将sample基因型recode（供体亲本为2，轮回亲本为0，杂合为1）
dd3 <- dd2 %>% filter((donor.GT == paste(REF, REF, sep = "/") & recurrent.GT == paste(ALT, ALT, sep = "/")) | 
                        (donor.GT == paste(ALT, ALT, sep = "/") & recurrent.GT == paste(REF, REF, sep = "/"))) %>% 
  mutate(code = if_else(sample.GT == donor.GT, 2, if_else(sample.GT == recurrent.GT, 0, 1)))
if (1) {
  # 导出SNPbinner输入文件格式
  snpbinner <- dd2 %>% filter((donor.GT == paste(REF, REF, sep = "/") & recurrent.GT == paste(ALT, ALT, sep = "/")) | 
                   (donor.GT == paste(ALT, ALT, sep = "/") & recurrent.GT == paste(REF, REF, sep = "/"))) %>% 
    mutate(genotype = if_else(sample.GT == donor.GT, "a", if_else(sample.GT == recurrent.GT, "b", if_else(sample.GT == str_c(REF, ALT, sep = "/"), "h", "-"))),
           marker = str_c(CHROM, POS, sep = "_")) %>% 
    select(CHROM, marker, POS, genotype)
  colnames(snpbinner) <- c("chr", "marker", "position(bp)", sample)
  bins <- tibble(chr = character(0), start = numeric(0), end = numeric(0), genotype = character(0))
  for (c in chr$CHROM) {
    snpbinner_chr <- snpbinner %>% filter(chr == c) %>% select(-chr)
    write_tsv(x = snpbinner_chr, file = paste(sample, "/", sample, ".", c, ".tsv", sep = ""))
    # python2 ~/tools/SNPbinner-1.0.0-GondaEtAl2019/snpbinner crosspoints --min-length 500000 --chrom-len $len --input ./$sample/$sample.$i.tsv --output ./$sample/$sample.$i.crosp.csv
    cmd1 <- paste("python2 ", path_to_snpbinner, " crosspoints --min-length 500000 --chrom-len ", len$Len[len$CHROM==c], " --input ./", sample, "/", sample, ".", c, ".tsv --output ./", sample, "/", sample, ".", c, ".crosp.csv", sep = "")
    # python2 ~/tools/SNPbinner-1.0.0-GondaEtAl2019/snpbinner bins --min-bin-size 10000 --input ./$sample/$sample.$i.crosp.csv --output ./$sample/$sample.$i.bins.csv
    cmd2 <- paste("python2 ", path_to_snpbinner, " bins --min-bin-size 10000 --input ./", sample, "/", sample, ".", c, ".crosp.csv --output ./", sample, "/", sample, ".", c, ".bins.csv", sep = "")
    # sed -i 's/^##//' ./$sample/$sample.$i.bins.csv
    cmd3 <- paste("sed -i 's/^##//' ./", sample, "/", sample, ".", c, ".bins.csv", sep = "")
    # csvtk transpose ./$sample/$sample.$i.bins.csv -o ./$sample/$i.bins.csv
    cmd4 <- paste("csvtk transpose ./", sample, "/", sample, ".", c, ".bins.csv -o ./", sample, "/", c, ".bins.csv", sep = "")
    system(cmd1, intern = FALSE)
    system(cmd2, intern = FALSE)
    system(cmd3, intern = FALSE)
    system(cmd4, intern = FALSE)
    bins_chr <- read_csv(file = paste("./", sample, "/", c, ".bins.csv", sep = ""), show_col_types = FALSE) %>% 
      mutate(chr = c) %>% 
      select(chr, start = `bin start`, end = `bin end`, genotype = all_of(sample))
    bins <- rbind(bins, bins_chr)
    system(command = paste("rm ./", sample, "/", c, ".bins.csv", sep = ""), intern = FALSE)
    system(command = paste("rm ./", sample, "/", sample, ".", c, ".bins.csv", sep = ""), intern = FALSE)
    system(command = paste("rm ./", sample, "/", sample, ".", c, ".crosp.csv", sep = ""), intern = FALSE)
    system(command = paste("rm ./", sample, "/", sample, ".", c, ".tsv", sep = ""), intern = FALSE)
  }
  # export table
  outbins <- bins %>% mutate(genotype = if_else(genotype == "b", "Recurrent", if_else(genotype == "a", "Donor", if_else(genotype == "h", "Heterozygous", "-"))))
  write_csv(x = outbins, file = paste("./", sample, "/", sample, ".bg.csv", sep = ""))
  bg_stat <- outbins %>% mutate(length = end - start) %>% group_by(genotype) %>% summarize(`length(bp)` = sum(length), rate = sum(length)/genome_size)
  write_csv(x = bg_stat, file = paste("./", sample, "/", sample, ".bg_stat.csv", sep = ""))
  # 
  genoforplot <- snpbinner %>% 
    select(CHROM = chr, x = `position(bp)`, genotype = all_of(sample)) %>% 
    right_join(chr, by = "CHROM")
  binsforplot <- bins %>% 
    select(CHROM = chr, START = start, END = end, genotype) %>% 
    filter(genotype == "a" | genotype == "h") %>% 
    left_join(chr, by = "CHROM")
  p <- ggplot(len) + 
    geom_rect(mapping = aes(xmin = 0, xmax = Len, ymin = y-1, ymax = y+1), fill = "gray60", color = "black") + 
    #geom_segment(data = genoforplot, mapping = aes(x = x, xend = x, y = y-1, yend = y, color = genotype), linewidth = 0.1) + 
    geom_rect(data = binsforplot, mapping = aes(xmin = START, xmax = END, ymin = y, ymax = y+1, fill = genotype)) + 
    scale_x_continuous(breaks = seq(0, max(len$Len) + 1000000, 10000000), 
                       labels = paste(seq(0, (max(len$Len) %/% 1000000) + 1, 10), "Mb", sep = " "),
                       limits = c(0, (max(len$Len+1000000) %/% 1000000) * 1000000), 
                       expand = c(0.025, 0)) + 
    scale_y_continuous(breaks = chr$y, labels = chr$LABEL, expand = c(0.025, 0)) + 
    scale_fill_manual(breaks = c("b", "h", "a"), values = c("blue", "green", "red")) + 
    scale_color_manual(breaks = c("b", "h", "a"), values = c("blue", "green", "red"), labels = c(recurrent, "Het", donor)) + 
    guides(x = guide_axis_truncated(), y = guide_axis_truncated(), fill = "none", 
           color = guide_legend(override.aes = list(linewidth = 3))) + 
    labs(x = NULL, y = NULL, color = "Genotype", fill = NULL) + 
    theme_half_open()
  if (require(ggrastr, quietly = TRUE)) {
    p <- p + rasterise(geom_segment(data = genoforplot, mapping = aes(x = x, xend = x, y = y-1, yend = y, color = genotype), linewidth = 0.1), dpi = 500)
  } else {
    p <- p + geom_segment(data = genoforplot, mapping = aes(x = x, xend = x, y = y-1, yend = y, color = genotype), linewidth = 0.1)
  }
  ggsave(p, filename = paste("./", sample, "/", sample, ".bg.png", sep = ""), width = 8, height = 7, dpi = 500)
  ggsave(p, filename = paste("./", sample, "/", sample, ".bg.pdf", sep = ""), width = 8, height = 7)
}

if (0) {
  dd3 <- dd2 %>% filter((donor.GT == paste(REF, REF, sep = "/") & recurrent.GT == paste(ALT, ALT, sep = "/")) | 
                          (donor.GT == paste(ALT, ALT, sep = "/") & recurrent.GT == paste(REF, REF, sep = "/"))) %>% 
    mutate(code = if_else(sample.GT == donor.GT, 2, if_else(sample.GT == recurrent.GT, 0, 1)))
  genotype <- tibble()
  for (c in chr$CHROM) {
    wind_sum <- dd3 %>% filter(CHROM == c) %>% arrange(POS) %>% 
      mutate(group = as.numeric(rownames(.)) %/% 15 + 1) %>%
      group_by(group, CHROM) %>%
      summarise(start = min(POS), end=max(POS), code_sum = sum(code)) %>%
      ungroup() %>% select(-group)
    wind_geno <- mutate(wind_sum, code = ifelse(code_sum < 6, 0, ifelse(code_sum > 24, 2, 1)))
    wind_geno$fix <- fixgeno.func(wind_geno$code, fix.size = fixSize)
    genotype <- rbind(genotype, wind_geno)
    #bin_map <- gather(wind_geno, key='method', value='value', code_sum:fix)
  }
  
  introgressionRegion <- getIntrogressionRegion(data = genotype)
  
  # draw figure
  genoforplot <- genotype %>% left_join(chr, by = "CHROM") %>% select(-LABEL) %>% mutate(x = (start + end) / 2)
  intrforplot <- introgressionRegion %>% left_join(chr, by = "CHROM") %>% select(-LABEL)
  p_bg <- ggplot(len) + 
    geom_rect(mapping = aes(xmin = 0, xmax = Len, ymin = y-1, ymax = y+1), fill = "gray60", color = "black") + 
    geom_segment(data = genoforplot, mapping = aes(x = x, xend = x, y = y-1, yend = y, color = code_sum)) + 
    geom_rect(data = intrforplot, mapping = aes(xmin = START, xmax = END, ymin = y, ymax = y+1, fill = as.factor(GENO)), linewidth = 3) + 
    scale_x_continuous(breaks = seq(0, max(len$Len) + 1000000, 10000000), 
                       labels = paste(seq(0, (max(len$Len) %/% 1000000) + 1, 10), "Mb", sep = " "),
                       limits = c(0, (max(len$Len+1000000) %/% 1000000) * 1000000), 
                       expand = c(0.025, 0)) + 
    scale_y_continuous(breaks = chr$y, labels = chr$LABEL, expand = c(0.025, 0)) + 
    scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 15) + 
    scale_fill_manual(breaks = c(1, 2), values = c("green", "red"), labels = c("Het", "Homo")) + 
    guides(x = guide_axis_truncated(), y = guide_axis_truncated()) + 
    labs(x = NULL, y = NULL, color = NULL, fill = NULL) + 
    theme_half_open()
  ggsave(p_bg, filename = paste(sample, "background.png", sep = "."), width = 8, height = 7, dpi = 500)
  ggsave(p_bg, filename = paste(sample, "background.pdf", sep = "."), width = 8, height = 7)
  
  bin_map <- genotype %>% gather(key = "method", value = "value", code_sum:fix) %>% right_join(chr, by = "CHROM")
  p <- ggplot(bin_map, aes(x=(start+end)/2, y=value)) + 
    geom_point() + 
    labs(x = NULL, y = NULL) + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing.x = unit(1, "mm")) + 
    facet_grid(method~LABEL , scales = "free", space = "free_x") 
  ggsave(p, filename = paste(sample, "pdf", sep = "."), width = 25, height = 5)
  ggsave(p, filename = paste(sample, "png", sep = "."), width = 25, height = 5, dpi = 500)
}

# self-defining function
fixgeno.func <- function(w.geno, fix.size = NULL){
  wind.geno.rle <- rle(w.geno)
  error.id <- which(wind.geno.rle$lengths < fix.size)
  for(i in error.id){
    left.id <- sum(wind.geno.rle$lengths[1:i]) - wind.geno.rle$lengths[i]
    right.id <- sum(wind.geno.rle$lengths[1:i])
    if(i==1){ 
      w.geno[(left.id+1):right.id] <-  w.geno[right.id+1]
    }else{
      w.geno[(left.id+1):right.id] <- w.geno[left.id]
    }
  }
  return(w.geno)
}

getIntrogressionRegion <- function(data) {
  introgressionRegion <- tibble(CHROM = character(0), START = numeric(0), END = numeric(0), GENO = integer(0))
  for (c in unique(data$CHROM)) {
    df <- data %>% filter(CHROM == c) %>% arrange(start)
    start <- Inf
    for (i in seq_along(rownames(df))) {
      if (i > 1) {
        if (df$fix[i] != df$fix[i-1] & df$fix[i-1]) {
          end <- mean(df$end[i:(i+1)])
          introgressionRegion <- rbind(introgressionRegion, tibble(CHROM = c, START = start, END = end, GENO = df$fix[i-1]))
        }
      }
      if (i == nrow(df) & df$fix[i]) {
        end <- df$end[i]
        introgressionRegion <- rbind(introgressionRegion, tibble(CHROM = c, START = start, END = end, GENO = df$fix[i]))
      }
      if (df$fix[i]) {
        if (i == 1) {
          start <- df$start[i]
        } else if (i > 1 & df$fix[i] != df$fix[i-1]) {
          start <- mean(df$start[(i-1):i])
        }
      }
    }
  }
  introgressionRegion
}
