#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("Generates plots of alignment data produced by show-coords.")

# Add command line arguments
p <- add_argument(parser = p, arg = "--input", short = "-i", help = "coords file from mummer program 'show-coords'", type="character")
p <- add_argument(parser = p, arg = "--out", short = "-o", help = "outfile prefix", type = "character")
p <- add_argument(parser = p, arg = "--ref", short = "-r", help = "reference name shown in plot", type = "character", default = "Reference")
p <- add_argument(parser = p, arg = "--query", short = "-q", help = "query name shown in plot", type = "character", default = "Query")
p <- add_argument(parser = p, arg = "--refIDs", short = "-R", 
                  help = "a tab-separated file containing two columns, first column is IDs in reference genome, second column is names to shown in plot, 
                  or a string like 'scaffoldA01:A01,scaffoldA03:A03'. If there is no file named that, we will threat it as a string", 
                  type = "character", default = NULL)
p <- add_argument(parser = p, arg = "--queryIDs", short = "-Q", 
                  help = "a tab-separated file containing two columns, first column is IDs in query genome, second column is names to shown in plot, 
                  or a string like 'scaffoldA01:A01,scaffoldA03:A03'. If there is no file named that, we will threat it as a string", 
                  type = "character", default = NULL)
p <- add_argument(parser = p, arg = "--refLen", short = "-l", 
                  help = "a tab-separated file containing two columns, first column is IDs in reference genome, second column is chr length", 
                  default = NULL)
p <- add_argument(parser = p, arg = "--queryLen", short = "-L", 
                  help = "a tab-separated file containing two columns, first column is IDs in query genome, second column is chr length", 
                  default = NULL)
p <- add_argument(parser = p, arg = "--min-query-length", short = "-M", help = "filter queries with total alignments less than cutoff X bp", 
                  type = "integer", default = 400000)
p <- add_argument(parser = p, arg = "--min-alignment-length", short = "-m", help = "filter alignments less than cutoff X bp", type = "integer", 
                  default = 2000)
p <- add_argument(parser = p, arg = "--min-identity", short = "-s", help = "filter alignments with identity less than X %", type = "numeric", default = 90)
p <- add_argument(parser = p, arg = "--color-by", short = "-c", help = "turn on color alignments by 'direction' or 'identity', no color if not assign", 
                  type = "character", default = NULL)
p <- add_argument(parser = p, arg = "--size", short = "-S", help = "line width of alignments in figure", type = "numeric", default = NULL)
p <- add_argument(parser = p, arg = "--width", short = "-W", help = "plot width (inches)", type = "numeric", default = 10)
p <- add_argument(parser = p, arg = "--height", short = "-H", help = "plot height (inches)", type = "numeric", default = 10)

# Parse the command line arguments
argv <- parse_args(p)

input_filename <- argv$input
output_filename <- argv$out
ref <- argv$ref
query <- argv$query
refIDs_filename <- argv$refIDs
queryIDs_filename <- argv$queryIDs
refLen_filename <- argv$refLen
queryLen_filename <- argv$queryLen
min_query_aln <- argv$min_query_length
min_align <- argv$min_alignment_length
min_idy <- argv$min_identity
color_by <- argv$color_by # "direction", "identity" or NULL
linewidth <- argv$size
plot_width <- argv$width
plot_height <- argv$height

# test
if (FALSE) {
  input_filename <- "./test.coords"
  output_filename <- "outprefix"
  ref <- "ZS11"
  query <- "Y317"
  refIDs_filename <- "refIDs.txt"
  queryIDs_filename <- "queryIDs.txt"
  refLen_filename <- "refLen.txt"
  queryLen_filename <- "queryLen.txt"
  min_query_aln <- 400000
  min_align <- 10000
  min_idy <- 90
  color_by <- "direction" # "direction", "identity"
  plot_width <- 10
  plot_height <- 10
}

library(tidyverse)

if (is.na(ref)) {
  ref <- "Reference"
}
if (is.na(query)) {
  query <- "Query"
}

# read in alignments
alignments = read.table(input_filename, stringsAsFactors = F, skip = 5)
alignments = alignments[,-c(3,6,9,11)]

# set column names
colnames(alignments) = c("refStart", "refEnd", "queryStart", "queryEnd", "lenAlnRef", "lenAlnQuery", "pctIDY", "refID", "queryID")
alignments <- as_tibble(alignments)

# replace ID to Name
## filter refs and querys contained in refIDs and queryIDs
if (is.na(refIDs_filename)) {
  refIDs <- alignments %>% distinct(refID) %>% mutate(refName = refID)
} else {
  if (file.exists(refIDs_filename)) {
    refIDs <- read_tsv(file = refIDs_filename, col_names = c("refID", "refName"), col_types = "cc")
  } else {
    refIDs <- tibble(
      pair = str_split_1(refIDs_filename, ",")
    ) %>%
      separate(pair, into = c("refID", "refName"), sep = ":")
  }
  if (nrow(refIDs) != length(unique(refIDs$refID)) || nrow(refIDs) != length(unique(refIDs$refName))) {
    stop("refID or refName is not unique!")
  }
}
if (is.na(queryIDs_filename)) {
  queryIDs <- alignments %>% distinct(queryID) %>% mutate(queryName = queryID)
} else {
  if (file.exists(queryIDs_filename)) {
    queryIDs <- read_tsv(file = queryIDs_filename, col_names = c("queryID", "queryName"), col_types = "cc")
  } else {
    queryIDs <- tibble(
      pair = str_split_1(queryIDs_filename, ",")
    ) %>%
      separate(pair, into = c("queryID", "queryName"), sep = ":")
  }
  if (nrow(queryIDs) != length(unique(queryIDs$queryID)) || nrow(queryIDs) != length(unique(queryIDs$queryName))) {
    stop("queryID or queryName is not unique!")
  }
}

# refID and queryID rename
alignments <- alignments %>% 
  filter(refID %in% refIDs$refID & queryID %in% queryIDs$queryID) %>% 
  left_join(refIDs, by = "refID") %>% 
  left_join(queryIDs, by = "queryID") %>% 
  select(-refID, -queryID) %>% 
  rename(refID = refName, queryID = queryName)

# filter
## filter queries by alignment length, for now include overlapping intervals

###########################################################################
## if `queryIDs_filename` is not NA, this step will only keep query ID   ##
## if `queryAlnLenTotal >= min_query_aln` if we want to keep all query   ##
## ID in the `queryIDs_filename` regardless of whether `queryAlnLenTotal ##
## >= min_query_aln`, we should modifiy these code below                 ##
###########################################################################
queryLenAgg <- alignments %>% group_by(queryID) %>% 
  summarise(queryAlnLenTotal = sum(lenAlnQuery)) %>% 
  filter(queryAlnLenTotal >= min_query_aln) %>% 
  pull(queryID)
alignments <- alignments %>% filter(queryID %in% queryLenAgg)
## filter alignment by length and identity
alignments <- alignments %>% filter(lenAlnQuery >= min_align, pctIDY >= min_idy)
## update `queryLenAgg`, if any query in `queryLenAgg` have no alignment longer than `min_align` and more than `min_idy` identity, then remove it
queryLenAgg <- alignments %>% distinct(queryID) %>% pull(queryID)

# sort queryID
## queryID in coords is same in query genome
## filter queryIDs according `queryLenAgg`
queryIDsToKeepOrdered <- queryIDs %>% filter(queryName %in% queryLenAgg) %>% pull(queryName)
alignments$queryID <- factor(alignments$queryID, levels = queryIDsToKeepOrdered)

# sort refID
correspond <- alignments %>% 
  group_by(refID, queryID) %>% 
  summarise(refAlnLenTotal = sum(lenAlnRef)) %>% 
  group_by(queryID) %>% 
  slice_max(refAlnLenTotal, with_ties = FALSE) %>% 
  select(-refAlnLenTotal) %>% ungroup()
if (is.na(refIDs_filename)) {
  ## sometimes, there are multiple reference IDs corresponding to one query ID, so `unique()`
  refIDsToKeepOrdered <- correspond %>% pull(refID) %>% unique()
  ## 同时更新refIDs中的顺序
  refIDs <- refIDs %>% filter(refName %in% refIDsToKeepOrdered) %>% 
    mutate(refName = factor(refName, levels = refIDsToKeepOrdered)) %>%
    arrange(refName) %>% 
    mutate(refName = as.character(refName))
} else {
  refIDsToKeepOrdered <- refIDs %>% pull(refName)
}

# filter alignments to remove refID and queryID not in refIDsToKeepOrdered and queryIDsToKeepOrdered
alignments <- alignments %>% filter(refID %in% refIDsToKeepOrdered, queryID %in% queryIDsToKeepOrdered)

alignments$refID <- factor(alignments$refID, levels = refIDsToKeepOrdered)

# ref length and query length
if (is.na(refLen_filename)) {
  refLen <- alignments %>% group_by(refID) %>% summarise(refLen = max(refStart, refEnd))
} else {
  refLen <- read_tsv(file = refLen_filename, col_names = c("refID", "refLen"), col_types = "ci")
  refLen <- refIDs %>% left_join(refLen) %>% select(refID = refName, refLen)
  refLen$refID <- factor(refLen$refID, levels = refIDsToKeepOrdered)
}
if (is.na(queryLen_filename)) {
  queryLen <- alignments %>% group_by(queryID) %>% summarise(queryLen = max(queryStart, queryEnd))
} else {
  queryLen <- read_tsv(file = queryLen_filename, col_names = c("queryID", "queryLen"), col_types = "ci")
  queryLen <- queryIDs %>% left_join(queryLen) %>% select(queryID = queryName, queryLen)
  queryLen$queryID <- factor(queryLen$queryID, levels = queryIDsToKeepOrdered)
}
refLen <- refLen %>% 
  mutate(refCumSumLen = cumsum(refLen), 
         refBreaks = refCumSumLen - refLen/2)
queryLen <- queryLen %>% 
  mutate(queryCumSumLen = cumsum(queryLen),
         queryBreaks = queryCumSumLen - queryLen/2)

# calc weighted mean percent identity
weightedAveIDY <- alignments %>% 
  group_by(refID, queryID) %>% 
  summarise(weighted_avg_identity = weighted.mean(pctIDY, (lenAlnRef+lenAlnQuery)/2)) %>% 
  ungroup()
alignments <- correspond %>% left_join(weightedAveIDY, by = c("refID", "queryID")) %>% right_join(alignments, by = c("refID", "queryID"))

# calc new start and end, calc strand
#rL <- refLen$refLen[1]
#qL <- queryLen$queryLen[1]
alignments <- alignments %>% left_join(refLen, by = "refID") %>% left_join(queryLen, by = "queryID")
alignments <- alignments %>% 
  mutate(newRefStart = refStart + refCumSumLen - refLen, 
         newRefEnd   = refEnd   + refCumSumLen - refLen,
         newQueryStart = queryStart + queryCumSumLen - queryLen, 
         newQueryEnd   = queryEnd   + queryCumSumLen - queryLen, 
         strand = if_else(is.na(weighted_avg_identity), NA_character_, # 和weighted_avg_identity保持一致，只在主要ref-query上加颜色，使用 NA_character_ 表示字符类型的 NA
                          if_else(sign(refStart-refEnd)*sign(queryStart-queryEnd) > 0, "+", "-"))) %>% # (refStart-refEnd)*(queryStart-queryEnd)数值有可能很大，超出integer的限制，因此用sign
  select(-refLen, -refCumSumLen, -refBreaks, -queryLen, -queryCumSumLen, -queryBreaks)

# plot
auto_label <- function(x) {
  if (max(x) < 1e6) {
    # 使用 kb 单位
    scales::label_number(scale = 1e-3, suffix = "kb")(x)
  } else {
    # 使用 Mb 单位
    scales::label_number(scale = 1e-6, suffix = "Mb")(x)
  }
}
p <- ggplot(alignments) + 
  geom_vline(data = refLen, mapping = aes(xintercept = refCumSumLen), linetype = "dashed", color = "dimgrey") + 
  geom_hline(data = queryLen, mapping = aes(yintercept = queryCumSumLen), linetype = "dashed", color = "dimgrey") + 
  #geom_segment(aes(x = newRefStart, xend = newRefEnd, y = newQueryStart, yend = newQueryEnd, color = weighted_avg_identity), linewidth = 1.5) + 
  scale_x_continuous(expand = c(0, 0), breaks = refLen$refBreaks, labels = refLen$refID, name = ref, 
                     sec.axis = sec_axis(~., labels = auto_label), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), breaks = queryLen$queryBreaks, labels = queryLen$queryID, name = query, 
                     sec.axis = sec_axis(~., labels = auto_label), limits = c(0, NA)) + 
  #scale_color_gradient(low = "blue", high = "red", name = "Identity") + 
  coord_fixed() + 
  theme(panel.border = element_rect(color = "black", fill = NA), 
        panel.background = element_rect(fill = NA), 
        axis.text.y.left = element_text(angle = 90, hjust = 0.5), 
        axis.text.y.right = element_text(angle = 90, hjust = 0.5))

if (is.na(color_by)) {
  color_by <- "noColor"
} else if (color_by != "identity" & color_by != "direction") {
  warning("'--color-by' is not 'identity', 'direction' or default!!! There will be no color assigning to alignment")
}
if (is.na(linewidth)) {
  linewidth <- NULL
}
p_color <- switch(
  color_by,
  direction = {
    ## colored by strand
    p + geom_segment(aes(x = newRefStart, xend = newRefEnd, y = newQueryStart, yend = newQueryEnd, color = strand), linewidth = linewidth) +
      scale_color_manual(breaks = c("+", "-"), values = c("red", "blue"), labels = c("Forward", "Reverse"), name = "Orientation", na.value = "black")
  },
  identity = {
    ## colored by identity
    p + geom_segment(aes(x = newRefStart, xend = newRefEnd, y = newQueryStart, yend = newQueryEnd, color = weighted_avg_identity), linewidth = linewidth) +
      scale_color_gradient(low = "blue", high = "red", name = "Identity", na.value = "black")
  },
  {
    ## no color
    p + geom_segment(aes(x = newRefStart, xend = newRefEnd, y = newQueryStart, yend = newQueryEnd), color = "black", linewidth = linewidth)
  }
)
ggsave(p_color, filename = paste(output_filename, "png", sep = "."), width = plot_width, height = plot_height, units = "in", dpi = 500)
ggsave(p_color, filename = paste(output_filename, "pdf", sep = "."), width = plot_width, height = plot_height, units = "in")
