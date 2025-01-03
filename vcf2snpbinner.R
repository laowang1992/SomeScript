#!/usr/bin/env Rscript
# Author: Wang Pengfei <wangpf0608@126.com>

# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for converting vcf to table of snpbinner. genotype same as parent_1 is designated 'a', genotype same as parent_2 is designated 'b', heterozygous genotype is designated 'h', missing genotype is designated '-'")

## Add command line arguments
##
p <- add_argument(parser = p, arg = "--input", short = "-i", help = "vcf or vcf.gz file containing two parents and progeny lines", type = "character")
p <- add_argument(parser = p, arg = "--out", short = "-o", help = "output file prefix", type = "character")
#
p <- add_argument(parser = p, arg = "--parent1", short = "-p", help = "name of parent_1", type = "character")
p <- add_argument(parser = p, arg = "--parent2", short = "-P", help = "name of parent_2", type = "character")
#
p <- add_argument(parser = p, arg = "--minDP_p1", short = "-d", help = "Minimum depth of parent_1", type = "numeric", default = 5)
p <- add_argument(parser = p, arg = "--minDP_p2", short = "-D", help = "Minimum depth of parent_2", type = "numeric", default = 5)
#
p <- add_argument(parser = p, arg = "--max_missing", short = "-m", help = "Maximum missing rate of SNP", type = "numeric", default = 0.3)

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
max_missing <- argv$max_missing

test <- FALSE
if (test) {
  vcfFile <- "../00.data/Genotype/filter.snps.vcf.gz"
  out_prefix <- "MelonF2"
  parent1 <- "DW"
  parent2 <- "HPM"
  minDP_p1 <- 5
  minDP_p2 <- 5
  max_missing <- 0.3
}

vcf <- read.vcfR(file = vcfFile)
gt <- extract.gt(vcf,element="GT") %>% gsub(pattern = "\\|", replacement = "/") %>% as.data.frame() %>% rownames_to_column(var = "ID") %>% as_tibble()
dp <- extract.gt(vcf,element="DP", as.numeric = T) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>% as_tibble()
samples <- colnames(gt)
samples <- samples[!(samples %in% c(parent1, parent2, "ID"))]

gt_parent <- gt %>% select(ID, p1_gt = all_of(parent1), p2_gt = all_of(parent2))
dp_parent <- dp %>% select(ID, p1_dp = all_of(parent1), p2_dp = all_of(parent2))
parent <- gt_parent %>% left_join(dp_parent, by = "ID") %>% 
  filter(p1_dp >= minDP_p1, p2_dp >= minDP_p2) %>% 
  filter((p1_gt=="1/1" | p1_gt=="0/0") & (p2_gt=="1/1" | p2_gt=="0/0") & p1_gt != p2_gt)

RIL_gt <- gt %>% select(ID, all_of(samples)) %>% 
  # filter missing rate
  mutate(na_count = rowSums(is.na(.)),                # 计算每行的NA数量
         na_percent = na_count / length(samples)) %>% 
  filter(na_percent < max_missing) %>% 
  select(-na_count, -na_percent)

# inner_join
gt <- parent %>% select(-p1_dp, -p2_dp) %>%
  inner_join(RIL_gt, by = "ID")
# recode
df <- gt %>% 
  mutate(across(c(-ID, -p1_gt, -p2_gt), ~ case_when(
    . == p1_gt ~ "a",          # 与p1_gt相同的元素替换为"a"
    . == p2_gt ~ "b",          # 与p2_gt相同的元素替换为"b"
    . == "0/1" ~ "h",          # 等于"0/1"的元素替换为"h"
    is.na(.) ~ "-",            # NA替换为"-"
    TRUE ~ "-"                    # 其他值保持不变
  ))) %>% select(-p1_gt, -p2_gt)

reform_gt <- df %>% select(ID) %>% 
  mutate(`position(bp)` = str_match(string = ID, pattern = "\\d+$")[,1]) %>% 
  left_join(df, by = "ID") %>% rename(marker = ID)
write_tsv(x = reform_gt, file = paste(out_prefix, "recode.txt", sep = "."))

# stat
##
stat <- df %>%
  pivot_longer(cols = -ID, names_to = "sampleID", values_to = "geno") %>%
  group_by(sampleID, geno) %>%
  summarize(count = n(), .groups = "drop") %>%
  # 完成所有可能的value（a, b, h, -）
  complete(sampleID, geno = c("a", "b", "h", "-"), fill = list(count = 0)) %>%
  pivot_wider(names_from = geno, values_from = count, values_fill = list(count = 0)) %>% 
  rename(a_count = a, b_count = b, het_count = h, missing_count = `-`) %>% 
  mutate(a_percent = a_count/(a_count+b_count+het_count), 
         b_percent = b_count/(a_count+b_count+het_count), 
         het_percent = het_count/(a_count+b_count+het_count), 
         missing_percent = missing_count/(a_count+b_count+het_count+missing_count)) %>% 
  select(sampleID, a_count, a_percent, b_count, b_percent, het_count, het_percent, missing_count, missing_percent)
## export
write_tsv(x = stat, file = paste(out_prefix, "stat.txt", sep = "."))
## figure
p_het <- ggplot(stat, aes(x = het_percent)) + 
  geom_histogram(bins = 50, fill = "orange") + 
  scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Individual heterozygosity rate", y = "Count") + 
  cowplot::theme_half_open()
ggsave(p_het, filename = paste(out_prefix, "heterozygosityRate.pdf", sep = "."), width = 5, height = 3.5)
ggsave(p_het, filename = paste(out_prefix, "heterozygosityRate.png", sep = "."), width = 5, height = 3.5, units = "in", dpi = 500)
p_miss <- ggplot(stat, aes(x = missing_percent)) + 
  geom_histogram(bins = 50, fill = "orange") + 
  scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Individual missing rate", y = "Count") + 
  cowplot::theme_half_open()
ggsave(p_miss, filename = paste(out_prefix, "missingRate.pdf", sep = "."), width = 5, height = 3.5)
ggsave(p_miss, filename = paste(out_prefix, "missingRate.png", sep = "."), width = 5, height = 3.5, units = "in", dpi = 500)
