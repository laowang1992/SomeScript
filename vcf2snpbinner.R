#!/usr/bin/env Rscript
# Author: Wang Pengfei <wangpf0608@126.com>

# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for converting vcf to table of snpbinner. genotype same as parent_1 is designated 'a', genotype same as parent_2 is designated 'b', heterozygous genotype is designated 'h', missing genotype is designated '-'")

## Add command line arguments
##
p <- add_argument(p, "--input", help = "vcf or vcf.gz file containing two parents and progeny lines", type = "character")
p <- add_argument(p, "--out", help = "output file prefix", type = "character")
#
p <- add_argument(p, "--parent1", help = "name of parent_1", type = "character")
p <- add_argument(p, "--parent2", help = "name of parent_2", type = "character")
#
p <- add_argument(p, "--minDP_p1", help = "Minimum depth of parent_1", type = "numeric", default = 5)
p <- add_argument(p, "--minDP_p2", help = "Minimum depth of parent_2", type = "numeric", default = 5)

## Parse the command line arguments
argv <- parse_args(p)

library(tidyverse)
library(vcfR)

vcfFile  <-  argv$input
out_prefix <-  argv$out
parent1  <-  argv$parent1
parent2  <-  argv$parent2
minDP_p1 <-  argv$minDP_p1
minDP_p2 <-  argv$minDP_p2

test <- FALSE
if (test) {
  vcfFile <- "250_Chr11.missing_maf.vcf.gz"
  out_prefix <- "250_Chr11"
  parent1 <- "875"
  parent2 <- "876"
  minDP_p1 <- 5
  minDP_p2 <- 5
}

vcf <- read.vcfR(file = vcfFile)
samples <- colnames(extract.gt(vcf,element="GT"))
samples <- samples[samples != parent1 & samples != parent2]
gt <- extract.gt(vcf,element="GT") %>% gsub(pattern = "\\|", replacement = "/") %>% as.data.frame() %>% rownames_to_column(var = "ID") %>% as_tibble()
dp <- extract.gt(vcf,element="DP", as.numeric = T) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>% as_tibble()

gt_parent <- gt %>% select(ID, p1_gt = all_of(parent1), p2_gt = all_of(parent2))
dp_parent <- dp %>% select(ID, p1_dp = all_of(parent1), p2_dp = all_of(parent2))
parent <- gt_parent %>% left_join(dp_parent, by = "ID") %>% 
  filter(p1_dp >= minDP_p1, p2_dp >= minDP_p2) %>% 
  filter((p1_gt=="1/1" | p1_gt=="0/0") & (p2_gt=="1/1" | p2_gt=="0/0") & p1_gt != p2_gt)

RIL_gt <- gt %>% select(ID, all_of(samples))

gt <- parent %>% select(-p1_dp, -p2_dp) %>%
  left_join(RIL_gt, by = "ID")

#reform_gt <- parent %>% select(ID) %>% mutate(ID2 = ID) %>% separate(col = ID2, into = c("chr", "position(bp)")) %>% select(-chr)
reform_gt <- parent %>% select(ID) %>% mutate(`position(bp)` = str_match(string = ID, pattern = "\\d+$")[,1])
stat <- tibble()

for (i in samples) {
  df <- gt %>% select(ID, p1_gt, p2_gt, sample = all_of(i)) %>% 
    mutate(sample = if_else(is.na(sample), "-", if_else(sample == p1_gt, "a", if_else(sample == p2_gt, "b", "h")))) %>%
    select(ID, sample)
  
  total <- nrow(df)
  n_miss <- df %>% filter(sample == "-") %>% nrow()
  miss_rate <- n_miss / total
  n_a <- df %>% filter(sample == "a") %>% nrow()
  a_rate <- n_a / (total - n_miss)
  n_b <- df %>% filter(sample == "b") %>% nrow()
  b_rate <- n_b / (total - n_miss)
  n_h <- df %>% filter(sample == "h") %>% nrow()
  het_rate <- n_h / (total - n_miss)
  stat_tmp <- tibble(sample = i, a_number = n_a, a_percent = a_rate, b_number = n_b, b_percent = b_rate, het_number = n_h, heterozygosity_rate = het_rate, missing_number = n_miss, missing_rate = miss_rate)
  stat <- rbind(stat, stat_tmp)
  
  colnames(df) <- c("ID", i)
  reform_gt <- reform_gt %>% left_join(df, by = "ID")
}

pdf(file = paste(out_prefix, "miss_het.stat.pdf", sep = "."), width = 5, height = 4)
hist(stat$heterozygosity_rate, breaks = 50, xlab = 'Individual heterozygosity rate', border = 'orange', main = NULL, col = 'orange')
hist(stat$missing_rate, breaks = 50, xlab = 'Individual missing rate', border = 'orange', main = NULL, col = 'orange')
dev.off()

stat <- stat %>% mutate(a_percent = scales::percent(a_percent, accuracy = 0.01),
                        b_percent = scales::percent(b_percent, accuracy = 0.01),
                        heterozygosity_rate = scales::percent(heterozygosity_rate, accuracy = 0.01),
                        missing_rate = scales::percent(missing_rate, accuracy = 0.01))

reform_gt <- reform_gt %>% rename(marker = ID)
write_tsv(x = reform_gt, file = paste(out_prefix, "recode.txt", sep = "."))
write_tsv(x = stat, file = paste(out_prefix, "stat.txt", sep = "."))
