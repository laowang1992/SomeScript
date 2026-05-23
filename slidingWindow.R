# A sliding window function in R.
# Author: Wang Pengfei <wangpf0608@126.com>

library(tidyverse)

slidingWindow <- function(df, winSize, winStep, groups, position, values, fun){
  # winStep应小于等于winSize
  if (winStep>winSize) {
    stop("winStep should not be bigger than winSize!")
  }
  df <- df %>% dplyr::select(groups = all_of(groups), position = all_of(position), all_of(values))
  chrLen <- df %>% group_by(groups) %>% summarise(Len = max(position))
  # 生成windows
  wid <- tibble()
  for (i in 1:nrow(chrLen)) {
    if (chrLen$Len[i]<=winSize) {
      start <- 1
      end <- chrLen$Len[i]
      #end <- winSize
    }else{
      end <- seq(from = winSize, to = chrLen$Len[i], by = winStep)
      start <- end - winSize + 1
      if (end[length(end)] < chrLen$Len[i]) {
        start <- c(start, start[length(start)]+winStep)
        end <- c(end, chrLen$Len[i])
      }
      #start = seq(from = 1, to = chrLen$Len[i]-winSize+winStep+1, by = winStep)
      #end = seq(from = winSize, to = chrLen$Len[i]+winStep, by = winStep)
    }
    wid_tmp <- tibble(groups = chrLen$groups[i], Start = start, End = end)
    wid <- rbind(wid, wid_tmp)
  }
  # 添加统计列（N，sum）
  # 分窗口统计
  cat(date(), ", Sliding window statistics start ...\n", sep = "")
  # 定义一个进度条
  width <- options()$width
  pb <- progress::progress_bar$new(
    format = 'Progress [:bar] :percent eta: :eta',
    total = nrow(wid), clear = FALSE, width = width
  )
  for (i in 1:nrow(wid)) {
    # 删掉na.omit()，这会导致 N 的数值小于实际的个数，后面apply添加na.rm = TRUE参数
    df_tmp <- df %>% filter(groups == wid$groups[i], position >= wid$Start[i], position <= wid$End[i])
    wid[i, "N"] <- nrow(df_tmp)
    # 添加na.rm = TRUE参数，处理含有NA的行
    x <- apply(df_tmp[, values], 2, fun, na.rm = TRUE)
    wid[i, paste(values, fun, sep = "_")] <- t(x[values])
    
    # 打印一个进度条
    ## 使用progress扩展包打印进度条
    pb$tick()
    ## 以下是一个手动打印进度条的代码
    #########################################################
    #cat('[', paste0(rep('#', i/nrow(wid)*width), collapse=''),
    #    paste0(rep('-', width - i/nrow(wid)*width), collapse=''),
    #    ']',
    #    round(i/nrow(wid)*100),'%')
    #if(i==nrow(wid))cat('\n', date(), ' DONE!\n', sep = "")
    #else cat('\r')
    ##########################################################
  }
  colnames(wid)[1:3] <- c(groups, "win_start", "win_end")
  return(wid)
}


slidingWindow_fast <- function(df, winSize, winStep, groups, position, values, fun) {
  # 1. 核心修复：安全转换函数对象，完美兼容 "mean" 或 mean
  actual_fun <- match.fun(fun)
  
  # 2. 智能提取函数文本名称，防止列名变成尴尬的 "_fun"
  fun_char <- deparse(substitute(fun))
  if (is.character(fun)) {
    fun_char <- fun
  } else if (fun_char == "fun" || grepl("^function", fun_char)) {
    fun_char <- "stat"
  }
  
  if (winStep > winSize) {
    stop("winStep should not be bigger than winSize!")
  }
  
  # 规范化列名
  df_clean <- df %>% 
    dplyr::select(groups = all_of(groups), position = all_of(position), all_of(values))
  
  cat(date(), " | 正在高效生成窗口坐标...\n", sep = "")
  
  # 生成 windows
  wid <- df_clean %>%
    group_by(groups) %>%
    summarise(Len = max(position, na.rm = TRUE), .groups = "drop") %>%
    group_by(groups) %>%
    reframe({
      if (Len <= winSize) {
        tibble(win_start = 1, win_end = Len)
      } else {
        ends <- seq(from = winSize, to = Len, by = winStep)
        starts <- ends - winSize + 1
        if (ends[length(ends)] < Len) {
          starts <- c(starts, starts[length(starts)] + winStep)
          ends <- c(ends, Len)
        }
        tibble(win_start = starts, win_end = ends)
      }
    }) %>%
    ungroup()
  
  cat(date(), " | 正在进行并行化滑窗统计 (Non-equi join)...\n", sep = "")
  
  # 3. 核心修复：这里调用 actual_fun 确保万无一失
  result <- wid %>%
    left_join(
      df_clean, 
      by = join_by(groups, win_start <= position, win_end >= position)
    ) %>%
    group_by(groups, win_start, win_end) %>%
    summarise(
      N = sum(!is.na(position)), 
      across(
        all_of(values), 
        ~ if(N[1] == 0) {NA} else {actual_fun(.x, na.rm = TRUE)}, 
        .names = "{.col}_{fun_char}"
      ),
      .groups = "drop"
    )
  
  colnames(result)[1] <- groups
  cat(date(), " | 统计完成！\n", sep = "")
  
  return(result)
}


# 测试这个函数
if (FALSE) {
  df <- read_tsv(file = "../../test/202209_QTLseq/03.Analysis/RY2P/RY2P.Depth_information.txt")
  
  winSize <- 2000000
  winStep <- 200000
  groups <- "CHROM"
  position <- "POS"
  values <- c("R.R3.depth", "R.qY.depth")
  fun <- "mean"
  
  system.time(x <- slidingWindow(df = df, winSize = winSize, winStep = winStep, groups = groups, position = position, values = values, fun = fun))
  system.time(x_fast <- slidingWindow_fast(df = df, winSize = winSize, winStep = winStep, groups = groups, position = position, values = values, fun = fun))
}


# 一些测试代码，以后可能更新到函数中
if(FALSE){
  fuc <- function(wid, values, df){
    df_tmp <- df %>% filter(groups == wid["groups"], position > wid["Start"], position < wid["End"])
    apply(df_tmp[, values], 2, sum)
  }
  x <- apply(wid, 1, fuc, values = values, df = df)
  
}
