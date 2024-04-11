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

if (FALSE) {
  popType <- "RIL"
  bulkSizeH <- 30
  bulkSizeL <- 30
  minDepth <- 5
  maxDepth <- 150
  repN <- 10000
}
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
getQTLseqCI(popType, bulkSizeH, bulkSizeL, minDepth, maxDepth, repN)
