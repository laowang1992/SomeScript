# A R script for drawing gene structure and the variation of this gene in 
# a population. A gtf file containing target gene and a vcf file containing 
# variation of this gene is needed. Hierarchical clustering algorithm was 
# adopted to distinguish different haplotype, the number of haplotype can 
# be designated according clusting result. Some polymorphism may exsit 
# within samples belonging to the same haplotype, you can divided them 
# into different haplotypes by setting more haplotypes.
# Author: Wang Pengfei <wangpf0608@126.com>
library(vcfR)
library(tidyverse)
library(ggnewscale)
library(cowplot)
library(ggh4x)
library(ggbio)
library(GenomicRanges)

# 
vcf_file <- "./example_data/Haplotype/input/myb28.C2.vcf.gz"
gtf_file <- "./example_data/Haplotype/input/BnaC02_MYB28.gtf"
pheno_file <- "./example_data/Haplotype/input/WH_GLU_blup.txt"
pheno_y_lab <- "Glucosinolate content (μmol/g)"

# 导入数据
## 基因结构
{
  gene <- read_tsv(file = gtf_file, col_names = F)
  x_upper <- max(gene$X4, gene$X5)
  x_lower <- min(gene$X4, gene$X5)
  # 构造画图用到的数据集
  gene <- GRanges("chr",IRanges(gene$X4,end=gene$X5,group=gene$X3))
  gene
}
## 变异信息
{
  geno <- read.vcfR(file = vcf_file)

  gt <- extract.gt(geno,element="GT", mask = FALSE, as.numeric = FALSE) %>% 
    as.data.frame() %>% rownames_to_column("ID") %>% as.data.frame()
  #extract.haps(geno, unphased_as_NA = F)
  
  pos <- getFIX(geno)[, c(1, 2, 3, 4, 5)] %>% as.data.frame()
  pos$POS <- as.numeric(pos$POS)
  
  df <- pos %>% left_join(gt, by = "ID") %>% 
    gather(key = sample, value = geno, -CHROM, -POS, -ID, -REF, -ALT) %>% 
    mutate(geno = if_else((geno == "0/0" | geno == "0|0"), 0, 
                          if_else((geno == "0/1" | geno == "0|1"), 1, 
                                  if_else((geno == "1/1" | geno == "1|1"), 2, NA)))) %>%
    filter(!is.na(geno))
  df_tmp <- df %>% spread(sample, geno)
  rownames(df_tmp) <- df_tmp$ID
  df_tmp <- df_tmp %>% select(-CHROM, -POS, -ID, -REF, -ALT) %>% t()
}
# 定义一些位置信息
{
  # x limit
  xlim_lower <- x_lower - (x_upper-x_lower)*0.33
  xlim_upper <- x_upper + (x_upper-x_lower)*0.1
  # 表示单倍型的位置
  hap_x_pos <- x_lower-(x_upper-x_lower)*0.05
  # haplotype位置
  hap_text_x_pos <- x_lower-(x_upper-x_lower)*0.20
  # haplotype width
  hap_width <- (x_upper-x_lower)*0.015
  # variation width
  #vari_width <- (x_upper-x_lower)/nrow(pos)*1.1
  # gene名称位置
  gene_name_x_pos <- x_lower-(x_upper-x_lower)*0.02
  # gene名称
  gene_name <- "italic(BnaC07.MYB28)"
}


# 层次聚类
d <- dist(df_tmp, method = "euclidean") # 计算各个样本点之间的欧氏距离
fit <- hclust(d) #进行Ward层次聚类
df$sample <- factor(df$sample, levels = fit$labels[fit$order])
ggplot(df, aes(x = ID, y = sample)) + 
  geom_tile(aes(fill = factor(geno))) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "NA",
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.ticks.y = element_blank())
group <- cutree(fit, k = 2) # 设定聚类个数
## 给聚成的类别加上红色边框
plot(fit, hang = -1) # 绘制树状图展示聚类结果
rect.hclust(fit, k=2, border="red")
dev.off()

hap <- tibble(sample = names(group)[fit$order], Haplotype = paste("Hap", group[fit$order], sep = ""))

Hap_out <- hap %>% left_join(df_tmp %>% as.data.frame() %>% rownames_to_column("sample"), by = "sample")
Hap_out
write_tsv(x = Hap_out, file = "./example_data/Haplotype/output/haplotype_with_geno.txt")


df <- Hap_out %>% rownames_to_column("y_pos")
df$y_pos <- as.numeric(df$y_pos)
df <- df %>% mutate(y_pos = as.numeric(y_pos), y_pos = (y_pos-nrow(df))/nrow(df)*20-2)
hap_text <- df %>% group_by(Haplotype) %>% summarise(y_min = min(y_pos), y_max = max(y_pos), n = n()) %>% 
  mutate(y_pos = (y_min+y_max)/2, Haplotype = paste(Haplotype, "\n(n = ", n, ")" , sep = "")) 
df <- df %>% gather(key = ID, value = geno, -y_pos, -sample, -Haplotype)
#lower <- min(pos$POS)
#upper <- max(pos$POS)
dd1 <- pos %>% mutate(x_pos = seq(from = x_lower, to = x_upper, length.out = nrow(pos))) %>%
  left_join(df, by = "ID") %>% as_tibble() %>% 
  mutate(geno = if_else(geno == 0, "0/0", if_else(geno == 1, "0/1", if_else(geno == 2, "1/1", NA))))
dd1 <- pos %>% mutate(x_pos = seq(from = x_lower + (x_upper-x_lower)/nrow(pos)/2, to = x_upper - (x_upper-x_lower)/nrow(pos)/2, length.out = nrow(pos))) %>%
  left_join(df, by = "ID") %>% as_tibble() %>% 
  mutate(geno = if_else(is.na(geno), "./.", if_else(geno == 0, "0/0", if_else(geno == 1, "0/1", if_else(geno == 2, "1/1", "./.")))))
#dd1$ID <- factor(dd1$ID, levels = pos$ID)
#dd1$sample <- factor(dd1$sample, levels = names(group)[fit$order])

P_geno <- ggplot() +
  geom_tile(data = dd1, mapping = aes(x = x_pos, y = y_pos, fill = geno)) +
  scale_fill_manual(breaks = c("0/0", "0/1", "1/1", "./."), values = c("#99CCCC", "#FFFFCC", "#FFCC99", "gray")) +
  new_scale_fill() +
  geom_tile(data = dd1, mapping = aes(x = hap_x_pos, y = y_pos, fill = factor(Haplotype)), width = hap_width) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  geom_text(data = hap_text, mapping = aes(x = hap_text_x_pos, y = y_pos, label = Haplotype), size = 4, hjust = 0.5) +
  scale_x_continuous(
    limits = c(xlim_lower, x_upper), 
    expand = c(0, 0)
    ) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave(P_geno, filename = "./example_data/Haplotype/output/Genotype.pdf", width = 5, height = 4, units = "in")
ggsave(P_geno, filename = "./example_data/Haplotype/output/Genotype.png", width = 5, height = 4, units = "in", dpi = 500)

pheno <- read_tsv(file = pheno_file, col_names = F) %>% dplyr::rename(sample = X1, pheno = X2)
pheno <- pheno %>% right_join(hap, by = "sample")
P_pheno <- ggplot(pheno, aes(x = Haplotype, y = pheno)) + 
  geom_boxplot(aes(fill = Haplotype), width = 0.6) +
  scale_fill_brewer(palette = "Set1") + 
  labs(x = NULL, y = pheno_y_lab) +
  cowplot::theme_half_open() +
  theme(legend.position = "none")
ggsave(P_pheno, filename = "./example_data/Haplotype/output/Phenotype.pdf", width = 3, height = 4, units = "in")
ggsave(P_pheno, filename = "./example_data/Haplotype/output/Phenotype.png", width = 3, height = 4, units = "in", dpi = 500)

p <- plot_grid(P_geno, P_pheno, rel_widths = c(2.5, 1))
ggsave(p, filename = "./example_data/Haplotype/output/GenotypeWithPhenotype.pdf", width = 8, height = 4.5, units = "in")
ggsave(p, filename = "./example_data/Haplotype/output/GenotypeWithPhenotype.png", width = 8, height = 4.5, units = "in", dpi = 500)

# 生成变异位置信息
snp_info <- pos %>% mutate(x_pos = seq(from = x_lower + (x_upper-x_lower)/nrow(pos)/2, to = x_upper - (x_upper-x_lower)/nrow(pos)/2, length.out = nrow(pos)))

P_gene_vari <- ggbio::autoplot(gene, aes(fill = group),geom="alignment") +
  theme_bw() +
  scale_x_continuous(limits = c(xlim_lower, x_upper),
                     expand = c(0, 0),
                     #breaks = c(seq(62971175,62975175,by=1000)),
                     position = "top"
                     )+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12, color = "black")) +
  guides(x = guide_axis_truncated()) +
  scale_fill_manual(values = c("#92d050", "#a6a6a6", "#a6a6a6")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate(geom = "text", x = gene_name_x_pos, y = 1, hjust = 1,
           label = gene_name, size = 4, parse = T) +
  coord_cartesian(clip = "off") +
  new_scale_fill() +
  geom_tile(data = dd1, mapping = aes(x = x_pos, y = y_pos, fill = geno)) +
  scale_fill_manual(breaks = c("0/0", "0/1", "1/1", "./."), values = c("#99CCCC", "#FFFFCC", "#FFCC99", "gray")) +
  new_scale_fill() +
  geom_tile(data = dd1, mapping = aes(x = hap_x_pos, y = y_pos, fill = factor(Haplotype)), width = hap_width) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  geom_text(data = hap_text, mapping = aes(x = hap_text_x_pos, y = y_pos, label = Haplotype), size = 4) +
  geom_segment(data = snp_info, mapping = aes(x = POS, xend = x_pos, y = 0.5, yend = -1.8), lty = "dashed", color = "grey") +
  geom_segment(data = snp_info, mapping = aes(x = POS, xend = POS, y = 1.4, yend = 0.6), color = "red")
pdf(file = "./example_data/Haplotype/output/GeneStructWithVari.pdf", width = 6, height = 5)
P_gene_vari
dev.off()
png(filename = "./example_data/Haplotype/output/GeneStructWithVari.png", width = 6, height = 5, units = "in", res = 500)
P_gene_vari
dev.off()
