#!/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for background analysis")

## Add command line arguments
#
p <- add_argument(parser = p, arg = "--createParameter", short = "-C", help = "create a Parameter csv file", flag = TRUE)
p <- add_argument(parser = p, arg = "--parameter", help = "Parameter csv file containing parameters, if the same parameter is listed in CMD line, this CMD line parameter will be omitted", type = "character")

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

p <- add_argument(p, "--width", help = "Plot width", type = "numeric", default = 8)
p <- add_argument(p, "--height", help = "Plot height", type = "numeric", default = 7)

p <- add_argument(parser = p, arg = "--path_to_snpbinner", short = "-p", help = "path to snpbinner", type = "character", default = "~/tools/snpbinner_python3/snpbinner")
p <- add_argument(parser = p, arg = "--path_to_csvtk", short = "-P", help = "path to csvtk, default uses system PATH", type = "character", default = "csvtk")

# Parse the command line arguments
argv <- parse_args(p)

# 创建一个参数模板文件
if (argv$createParameter) {
  cat("You will get a template file of parameter.\n")
  parameter <- data.frame(sample = argv$sample, donor = argv$donor, recurrent = argv$recurrent, 
                          minDdp = argv$minDdp, maxDdp = argv$maxDdp, 
                          minRdp = argv$minRdp, maxRdp = argv$maxRdp, 
                          width = argv$width, height = argv$height)
  write.csv(x = parameter, file = "Parameter.csv", na = "", row.names = F)
  quit(save = "no")
}

library(tidyverse)
library(cowplot)
#library(ggh4x)

input <- argv$input
chromosome <- argv$chromosome
length <- argv$length
path_to_snpbinner <- argv$path_to_snpbinner
path_to_csvtk <- argv$path_to_csvtk

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
  
  width <- 8
  height <- 7
  
  #fixSize <- 10
  
  path_to_snpbinner <- "~/tools/SNPbinner-1.0.0-GondaEtAl2019/snpbinner"
  csvtk <- "csvtk"
}

df <- read_tsv(file = input)

if (is.na(chromosome)) {
  chr <- df %>% distinct(CHROM) %>% mutate(LABEL = CHROM, y = rev(seq_along(rownames(.))*3))
} else {
  chr <- read_tsv(chromosome, col_names = c("CHROM", "LABEL"), show_col_types = FALSE) %>% mutate(y = rev(seq_along(rownames(.))*3))
}
if (is.na(length)) {
  len <- df %>% group_by(CHROM) %>% summarise(Len = max(POS)) %>% right_join(chr, by = "CHROM")
} else {
  len <- read_tsv(file = length, col_names = c("CHROM", "Len"), show_col_types = FALSE) %>% right_join(chr, by = "CHROM")
}
genome_size <- sum(len$Len)

# 
if (is.na(argv$parameter)) {
  cat("Get parameters from CMD line.\n")
  parameter <- data.frame(sample = argv$sample, donor = argv$donor, recurrent = argv$recurrent, 
                          minDdp = argv$minDdp, maxDdp = argv$maxDdp, 
                          minRdp = argv$minRdp, maxRdp = argv$maxRdp, 
                          width = argv$width, height = argv$height)
} else {
  cat("Get parameters from parameter csv file.\n")
  parameter <- read.csv(argv$parameter, header = T)
}

for (i in seq_along(parameter[,1])) {
  #outPrefix <- parameter$outPrefix[i]
  sample <- parameter$sample[i]
  donor <- parameter$donor[i]
  recurrent <- parameter$recurrent[i]
  minDdp <- parameter$minDdp[i]
  maxDdp <- parameter$maxDdp[i]
  minRdp <- parameter$minRdp[i]
  maxRdp <- parameter$maxRdp[i]
  width <- parameter$width[i]
  height <- parameter$height[i]
  
  # 创建输出文件夹
  dir.create(path = sample, showWarnings = F)
  
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
    #coord_cartesian(xlim = c(0, upper)) + 
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
    write_csv(snpbinner %>% select(-marker), file = paste("./", sample, "/", sample, ".genotype.csv", sep = ""))
    
    bins <- tibble(chr = character(0), start = numeric(0), end = numeric(0), genotype = character(0))
    for (c in chr$CHROM) {
      snpbinner_chr <- snpbinner %>% filter(chr == c) %>% select(-chr)
      write_tsv(x = snpbinner_chr, file = paste(sample, "/", sample, ".", c, ".tsv", sep = ""))
      # python ~/tools/SNPbinner-1.0.0-GondaEtAl2019/snpbinner crosspoints --min-length 500000 --chrom-len $len --input ./$sample/$sample.$i.tsv --output ./$sample/$sample.$i.crosp.csv
      cmd1 <- paste("python ", path_to_snpbinner, " crosspoints --min-length 500000 --chrom-len ", len$Len[len$CHROM==c], " --input ./", sample, "/", sample, ".", c, ".tsv --output ./", sample, "/", sample, ".", c, ".crosp.csv", sep = "")
      # python ~/tools/SNPbinner-1.0.0-GondaEtAl2019/snpbinner bins --min-bin-size 10000 --input ./$sample/$sample.$i.crosp.csv --output ./$sample/$sample.$i.bins.csv
      cmd2 <- paste("python ", path_to_snpbinner, " bins --min-bin-size 10000 --input ./", sample, "/", sample, ".", c, ".crosp.csv --output ./", sample, "/", sample, ".", c, ".bins.csv", sep = "")
      # sed -i 's/^##//' ./$sample/$sample.$i.bins.csv
      cmd3 <- paste("sed -i 's/^##//' ./", sample, "/", sample, ".", c, ".bins.csv", sep = "")
      # csvtk transpose ./$sample/$sample.$i.bins.csv -o ./$sample/$i.bins.csv
      cmd4 <- paste(path_to_csvtk, " transpose ./", sample, "/", sample, ".", c, ".bins.csv -o ./", sample, "/", c, ".bins.csv", sep = "")
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
    write_csv(x = outbins %>% filter(genotype != "Recurrent"), file = paste("./", sample, "/", sample, ".bg.csv", sep = ""))
    bg_stat <- outbins %>% mutate(length = end - start) %>% group_by(genotype) %>% summarize(`length(bp)` = sum(length), rate = sum(length)/genome_size, .groups = "drop") %>%
      complete(genotype = c("Donor", "Heterozygous", "Recurrent"), fill = list(`length(bp)` = 0, rate = 0))
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
      guides(x = guide_axis(cap = TRUE), y = guide_axis(cap = TRUE), fill = "none", 
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
}
