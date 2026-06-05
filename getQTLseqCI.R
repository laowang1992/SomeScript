#!/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for calculating QTLseq confidence interval")

## Add command line arguments
#
p <- add_argument(p, "--popType", help = "Population type, 'F2' or 'RIL'", type = "character")

p <- add_argument(p, "--bulkSizeH", help = "Bulk size with high phenotype", type = "numeric")
p <- add_argument(p, "--bulkSizeL", help = "Bulk size with low phenotype", type = "numeric")

p <- add_argument(p, "--minDepth", help = "Minimum depth for simulation", type = "numeric", default = 5)
p <- add_argument(p, "--maxDepth", help = "Maxmum depth for simulation", type = "numeric", default = 150)

p <- add_argument(p, "--repN", help = "how many times for simulation", type = "numeric", default = 10000)

# Parse the command line arguments
argv <- parse_args(p)

popType <- argv$popType
bulkSizeH <- argv$bulkSizeH
bulkSizeL <- argv$bulkSizeL

minDepth <- argv$minDepth
maxDepth <- argv$maxDepth

repN <- argv$repN


library(tidyverse)
library(cowplot)

if (FALSE) {
  # 第一个版本
  getQTLseqCI <- function(popType, bulkSizeH, bulkSizeL, minDepth = 5, maxDepth = 150, repN = 10000){
    #set.seed(123)
    dltIndex_CI <- data.frame(HB.DP = integer(0), LB.DP = integer(0), CI95upper = numeric(0), CI95lower = numeric(0), CI99upper = numeric(0), CI99lower = numeric(0))

    cat(date(), ", program start runing ...\n", sep = "")
    ##
    {
      n <- 1
      for (depthH in minDepth:maxDepth) {
        for (depthL in minDepth:maxDepth) {
          indexH <- vector(length = repN)
          indexL <- vector(length = repN)
          dltIndex <- vector(length = repN)
          # 以循环方式速度太慢
          #for (i in 1:repN) {
          #  if (popType == "RIL") {
          #    PH <- mean(sample(c(0,1), bulkSizeH, replace = TRUE))
          #    PL <- mean(sample(c(0,1), bulkSizeL, replace = TRUE))
          #  }else if (popType == "F2") {
          #    PH <- mean(sample(c(0, 0.5, 1), bulkSizeH, replace = TRUE, prob = c(0.25, 0.5, 0.25)))
          #    PL <- mean(sample(c(0, 0.5, 1), bulkSizeL, replace = TRUE, prob = c(0.25, 0.5, 0.25)))
          #  }
          #  indexH[i] <- mean(sample(c(0, 1), depthH, replace = TRUE, prob = c(1-PH, PH)))
          #  indexL[i] <- mean(sample(c(0, 1), depthL, replace = TRUE, prob = c(1-PL, PL)))
          #  dltIndex[i] <- indexL[i] - indexH[i]
          #}
          # 改为取随机数方式
          if (popType == "RIL") {
            PH <- rbinom(repN, bulkSizeH, 0.5) / bulkSizeH
            PL <- rbinom(repN, bulkSizeL, 0.5) / bulkSizeL
            #PH <- apply(rmultinom(repN, bulkSizeH, c(1, 1)) * c(1, 0) / bulkSizeH, 2, sum)
            #PL <- apply(rmultinom(repN, bulkSizeL, c(1, 1)) * c(1, 0) / bulkSizeL, 2, sum)
          } else if (popType == "F2") {
            PH <- apply(rmultinom(repN, bulkSizeH, c(1, 2, 1)) * c(1, 0.5, 0) / bulkSizeH, 2, sum)
            PL <- apply(rmultinom(repN, bulkSizeL, c(1, 2, 1)) * c(1, 0.5, 0) / bulkSizeL, 2, sum)
          }
          indexH <- rbinom(repN, depthH, PH) / depthH
          indexL <- rbinom(repN, depthL, PL) / depthL
          dltIndex <- indexH - indexL
          dltIndex_CI[n, ] <- c(depthH, depthL,
                                quantile(dltIndex, 0.975), quantile(dltIndex, 0.025),
                                quantile(dltIndex, 0.995), quantile(dltIndex, 0.005))
          n = n + 1
        }
        cat(date(), ", ", paste(depthH, maxDepth, sep = "/"), " have been done ...\n", sep = "")
      }
    }
    cat(date(), ", finish, export data ...\n", sep = "")


    df <- gather(dltIndex_CI, type, CI, -HB.DP, -LB.DP) %>%
      mutate(level = if_else(str_detect(type, "95"), "95 CI", "99 CI"),
             direction = if_else(str_detect(type, "upper"), "Upper CI", "Lower CI"))
    df$direction <- factor(df$direction, levels = c("Upper CI", "Lower CI"))
    P_CI <- ggplot(df, aes(x = HB.DP, y = LB.DP, fill = CI)) +
      geom_tile() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_distiller(palette = "RdBu") +
      labs(title = paste("CI",
                         paste(paste("indH", bulkSizeH, sep = ""),
                               paste("indL", bulkSizeL, sep = ""),
                               popType, sep = "_"),
                         paste("Depth", minDepth, maxDepth, sep = "_"),
                         paste("Rep", repN, sep = "_"), sep = "."),
           fill = NULL) +
      facet_grid(direction ~ level) +
      theme_half_open() +
      theme(strip.background = element_rect(fill = "#90EE90"),
            plot.title = element_text(hjust = 0.5))
    ggsave(P_CI, filename = paste("QTLseqCI",
                                  paste(paste("indH", bulkSizeH, sep = ""),
                                        paste("indL", bulkSizeL, sep = ""),
                                        popType, sep = "_"),
                                  paste("Depth", minDepth, maxDepth, sep = "_"),
                                  paste("Rep", repN, sep = "_"),
                                  "pdf", sep = "."),
           width = 6.5, height = 6)
    ggsave(P_CI, filename = paste("QTLseqCI",
                                  paste(paste("indH", bulkSizeH, sep = ""),
                                        paste("indL", bulkSizeL, sep = ""),
                                        popType, sep = "_"),
                                  paste("Depth", minDepth, maxDepth, sep = "_"),
                                  paste("Rep", repN, sep = "_"),
                                  "png", sep = "."),
           width = 6.5, height = 6, units = "in", dpi = 500)


    save(dltIndex_CI,
         file = paste("QTLseqCI",
                      paste(paste("indH", bulkSizeH, sep = ""),
                            paste("indL", bulkSizeL, sep = ""),
                            popType, sep = "_"),
                      paste("Depth", minDepth, maxDepth, sep = "_"),
                      paste("Rep", repN, sep = "_"),
                      "RData", sep = "."))
  }
  # 以上版本在“RIL 30 30 5 200 10000” “执行时间：278.360000000015秒”
}

# 第二个版本
getQTLseqCI <- function(popType, bulkSizeH, bulkSizeL, minDepth = 5, maxDepth = 150, repN = 10000){
  dltIndex_CI <- tibble(HB.DP = rep(minDepth:maxDepth, times = maxDepth-minDepth+1),
                        LB.DP = rep(minDepth:maxDepth, each = maxDepth-minDepth+1),
                        CI95upper = NA, CI95lower = NA, CI99upper = NA, CI99lower = NA)
  cat(date(), ", program start runing ...\n", sep = "")
  # 定义一个进度条
  width <- options()$width
  pb <- progress::progress_bar$new(
    format = 'Progress [:bar] :percent eta: :eta',
    total = nrow(dltIndex_CI), clear = FALSE, width = width
  )
  for (i in 1:nrow(dltIndex_CI)) {
    depthH <- dltIndex_CI[i, 1][[1]]
    depthL <- dltIndex_CI[i, 2][[1]]

    if (popType == "RIL") {
      PH <- rbinom(repN, bulkSizeH, 0.5) / bulkSizeH
      PL <- rbinom(repN, bulkSizeL, 0.5) / bulkSizeL
      #PH <- apply(rmultinom(repN, bulkSizeH, c(1, 1)) * c(1, 0) / bulkSizeH, 2, sum)
      #PL <- apply(rmultinom(repN, bulkSizeL, c(1, 1)) * c(1, 0) / bulkSizeL, 2, sum)
    } else if (popType == "F2") {
      PH <- apply(rmultinom(repN, bulkSizeH, c(1, 2, 1)) * c(1, 0.5, 0) / bulkSizeH, 2, sum)
      PL <- apply(rmultinom(repN, bulkSizeL, c(1, 2, 1)) * c(1, 0.5, 0) / bulkSizeL, 2, sum)
    }

    indexH <- rbinom(repN, depthH, PH) / depthH
    indexL <- rbinom(repN, depthL, PL) / depthL
    dltIndex <- indexH - indexL

    dltIndex_CI[i, c("CI95upper", "CI95lower", "CI99upper", "CI99lower")] <- t(
      c(quantile(dltIndex, 0.975), quantile(dltIndex, 0.025),
        quantile(dltIndex, 0.995), quantile(dltIndex, 0.005))
    )

    # 打印进度条
    pb$tick()
  }
  cat(date(), ", finish, export data ...\n", sep = "")


  df <- gather(dltIndex_CI, type, CI, -HB.DP, -LB.DP) %>%
    mutate(level = if_else(str_detect(type, "95"), "95 CI", "99 CI"),
           direction = if_else(str_detect(type, "upper"), "Upper CI", "Lower CI"))
  df$direction <- factor(df$direction, levels = c("Upper CI", "Lower CI"))
  P_CI <- ggplot(df, aes(x = HB.DP, y = LB.DP, fill = CI)) +
    geom_tile() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_distiller(palette = "RdBu") +
    #scale_fill_gradientn(colors = c("#4682B4", "white", "white", "#FF4500"),
    #                     values = c(0, (max(dltIndex_CI$CI95lower, dltIndex_CI$CI99lower) - min(dltIndex_CI$CI95lower, dltIndex_CI$CI99lower)) / (max(dltIndex_CI$CI95upper, dltIndex_CI$CI99upper) - min(dltIndex_CI$CI95lower, dltIndex_CI$CI99lower)),
    #                                (min(dltIndex_CI$CI95upper, dltIndex_CI$CI99upper) - min(dltIndex_CI$CI95lower, dltIndex_CI$CI99lower)) / (max(dltIndex_CI$CI95upper, dltIndex_CI$CI99upper) - min(dltIndex_CI$CI95lower, dltIndex_CI$CI99lower)), 1)) +
    labs(title = paste("CI",
                       paste(paste("indH", bulkSizeH, sep = ""),
                             paste("indL", bulkSizeL, sep = ""),
                             popType, sep = "_"),
                       paste("Depth", minDepth, maxDepth, sep = "_"),
                       paste("Rep", repN, sep = "_"), sep = "."),
         fill = NULL) +
    facet_grid(direction ~ level) +
    theme_half_open() +
    theme(strip.background = element_rect(fill = "#90EE90"),
          plot.title = element_text(hjust = 0.5))
  ggsave(P_CI, filename = paste("QTLseqCI",
                                paste(paste("indH", bulkSizeH, sep = ""),
                                      paste("indL", bulkSizeL, sep = ""),
                                      popType, sep = "_"),
                                paste("Depth", minDepth, maxDepth, sep = "_"),
                                paste("Rep", repN, sep = "_"),
                                "pdf", sep = "."),
         width = 6.5, height = 6)
  ggsave(P_CI, filename = paste("QTLseqCI",
                                paste(paste("indH", bulkSizeH, sep = ""),
                                      paste("indL", bulkSizeL, sep = ""),
                                      popType, sep = "_"),
                                paste("Depth", minDepth, maxDepth, sep = "_"),
                                paste("Rep", repN, sep = "_"),
                                "png", sep = "."),
         width = 6.5, height = 6, units = "in", dpi = 500)


  save(dltIndex_CI,
       file = paste("QTLseqCI",
                    paste(paste("indH", bulkSizeH, sep = ""),
                          paste("indL", bulkSizeL, sep = ""),
                          popType, sep = "_"),
                    paste("Depth", minDepth, maxDepth, sep = "_"),
                    paste("Rep", repN, sep = "_"),
                    "RData", sep = "."))
}
# 以上版本在“RIL 30 30 5 200 10000” “执行时间：264.410000000003秒”

# 一个加速版本
getQTLseqCI_fast <- function(popType, bulkSizeH, bulkSizeL, minDepth = 5, maxDepth = 150, repN = 10000){

  # 1. 计算所有独立的深度值
  depths <- minDepth:maxDepth
  n_depths <- length(depths)
  total_comb <- n_depths * n_depths

  cat(date(), " | Step 1: calculating allele frequency (PH & PL)...\n", sep = "")

  # 【优化一】将 PH 和 PL 移出循环，只计算一次！
  # 同时用高速的矩阵算术代替 F2 中的 apply(..., 2, sum)
  if (popType == "RIL") {
    PH <- rbinom(repN, bulkSizeH, 0.5) / bulkSizeH
    PL <- rbinom(repN, bulkSizeL, 0.5) / bulkSizeL
  } else if (popType == "F2") {
    mH <- rmultinom(repN, bulkSizeH, c(1, 2, 1))
    PH <- (mH[1, ] + 0.5 * mH[2, ]) / bulkSizeH

    mL <- rmultinom(repN, bulkSizeL, c(1, 2, 1))
    PL <- (mL[1, ] + 0.5 * mL[2, ]) / bulkSizeL
  } else {
    stop("Unknown population type, select 'RIL' or 'F2'")
  }

  cat(date(), " | Step 2: Precalculating SNP-index matrices across different depths...\n", sep = "")

  # 【优化二】预先生成所有独立深度的测序抽样矩阵 (行=repN, 列=独立深度数)
  # 这样 rbinom 的调用次数从 38416 次骤降到 196 次
  indexH_mat <- matrix(NA, nrow = repN, ncol = n_depths)
  indexL_mat <- matrix(NA, nrow = repN, ncol = n_depths)

  for(j in 1:n_depths) {
    indexH_mat[, j] <- rbinom(repN, depths[j], PH) / depths[j]
    indexL_mat[, j] <- rbinom(repN, depths[j], PL) / depths[j]
  }

  cat(date(), " | Step 3: calculating delta index and confidence interval...\n", sep = "")

  # 【优化三】利用矩阵索引和预分配内存，彻底告别 Tibble 按行写入
  h_indices <- rep(1:n_depths, times = n_depths)
  l_indices <- rep(1:n_depths, each = n_depths)

  # 创建原代码对应的基础网格
  dltIndex_CI <- tibble(
    HB.DP = rep(depths, times = n_depths),
    LB.DP = rep(depths, each = n_depths)
  )

  # 预分配存储置信区间的矩阵
  res_cl <- matrix(NA, nrow = total_comb, ncol = 4)

  pb <- progress::progress_bar$new(
    format = 'Progress [:bar] :percent eta: :eta',
    total = total_comb, clear = FALSE, width = options()$width
  )

  # 核心循环：现在内部只有极快的矩阵列减法和分位数计算
  for (i in 1:total_comb) {
    dltIndex <- indexH_mat[, h_indices[i]] - indexL_mat[, l_indices[i]]

    # names = FALSE 可以略微加速 quantile 函数
    res_cl[i, ] <- quantile(dltIndex, probs = c(0.975, 0.025, 0.995, 0.005), names = FALSE)
    pb$tick()
  }

  # 将计算结果快速合并回 tibble
  dltIndex_CI <- dltIndex_CI %>%
    mutate(
      CI95upper = res_cl[, 1],
      CI95lower = res_cl[, 2],
      CI99upper = res_cl[, 3],
      CI99lower = res_cl[, 4]
    )

  cat(date(), " | Step 4: finish calculating, export figures and data...\n", sep = "")

  # ========================================
  df <- gather(dltIndex_CI, type, CI, -HB.DP, -LB.DP) %>%
    mutate(level = if_else(str_detect(type, "95"), "95 CI", "99 CI"),
           direction = if_else(str_detect(type, "upper"), "Upper CI", "Lower CI"))
  df$direction <- factor(df$direction, levels = c("Upper CI", "Lower CI"))

  P_CI <- ggplot(df, aes(x = HB.DP, y = LB.DP, fill = CI)) +
    geom_tile() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_distiller(palette = "RdBu") +
    labs(title = paste("CI",
                       paste(paste("indH", bulkSizeH, sep = ""),
                             paste("indL", bulkSizeL, sep = ""),
                             popType, sep = "_"),
                       paste("Depth", minDepth, maxDepth, sep = "_"),
                       paste("Rep", repN, sep = "_"), sep = "."),
         fill = NULL) +
    facet_grid(direction ~ level) +
    theme_half_open() +
    theme(strip.background = element_rect(fill = "#90EE90"),
          plot.title = element_text(hjust = 0.5))

  file_base <- paste("QTLseqCI", paste(paste("indH", bulkSizeH, sep = ""), paste("indL", bulkSizeL, sep = ""), popType, sep = "_"), paste("Depth", minDepth, maxDepth, sep = "_"), paste("Rep", repN, sep = "_"), sep = ".")
  ggsave(P_CI, filename = paste(file_base, "pdf", sep = "."), width = 6.5, height = 6)
  ggsave(P_CI, filename = paste(file_base, "png", sep = "."), width = 6.5, height = 6, units = "in", dpi = 500)
  save(dltIndex_CI, file = paste(file_base, "RData", sep = "."))

  cat(date(), " | finish！\n", sep = "")
  return(dltIndex_CI)
}

getQTLseqCI_fast(popType, bulkSizeH, bulkSizeL, minDepth, maxDepth, repN)

## test
if (FALSE) {
  popType <- "RIL"
  bulkSizeH <- 30
  bulkSizeL <- 30
  minDepth <- 5
  maxDepth <- 150
  repN <- 1000
  system.time(getQTLseqCI(popType, bulkSizeH, bulkSizeL, minDepth, maxDepth, repN))
  system.time(getQTLseqCI_fast(popType, bulkSizeH, bulkSizeL, minDepth, maxDepth, repN))
}
