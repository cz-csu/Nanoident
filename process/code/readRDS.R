
library(foreach)
library(plyr)
library(dplyr)
library(tidyr)
# library(rgl) # Not needed
library(purrr)
library(stringr)
# 假设df是你的数据框
df <- data.frame(
  contig = rep("CP003732.1", 7),
  position = c(1,1,2,2,3,3,4),
  dir = c("fwd", "wd", "fwd", "wd", "fwd", "wd", "fwd"),
  strand = rep("t", 7),
  N_wga = c(34, 230, 49, 34, 53, 91, 43),
  N_nat = c(16, 45, 30, 8, 31, 17, 20),
  mean_diff = c(0.78726120296, 0.07191890693, 0.43346003977, -0.41532554330, 0.45829297534, -0.90710000478, 0.03529801069),
  t_test_pval = c(8.392534e-02, 9.554498e-01, 1.574673e-01, 4.387530e-01, 2.175693e-01, 5.123250e-01, 9.3557)
)
print(df)
# 创建一个包含每个position周围5个位点的mean_diff值的列表列
df <- df %>%
  #arrange(contig, dir,strand) %>%
  group_by(contig, dir,strand) %>%
  mutate(mean_diff_around = map(position, ~ mean_diff[pmax(.x - 5, 1):pmin(.x + 5, n())]))

res<- readRDS("analysis/NO_subset_with_ks_difference.RDS")
print(head(res,30))
res$ks_pval<- NULL
res <- res %>%
  mutate(position = as.integer(position)) %>%
  group_by(contig, dir, strand) %>%
  mutate(mean_nat_around = map(row_number(), ~ mean_nat[pmax(.x - 2, 1):pmin(.x + 2, n())]))%>%
  mutate(mean_wga_around = map(row_number(), ~ mean_wga[pmax(.x - 2, 1):pmin(.x + 2, n())])) %>%
  mutate(ks_pval = map2(mean_wga_around, mean_nat_around, ~ try(wilcox.test(unlist(.x), unlist(.y), correct=FALSE)$p.value, silent=TRUE)))
  #ungroup()
#res$mean_diff_around <- NULL
res$mean_nat_around <- NULL
res$mean_wga_around <- NULL
res$ks_pval=unlist(res$ks_pval)
res$ks_pval=as.numeric(res$ks_pval)
#sink("mid_data/TP_diff_ks.txt")
print(head(res$ks_pval,30))
saveRDS(res, "analysis/NO_subset_with_ks_difference2.RDS")
#analysis/TP_subset_with_ks_difference.RDS

  #mutate(mean_diff_around = map(row_number(), ~ mean_diff[pmax(.x - 2, 1):pmin(.x + 2, n())]))%>%
  #mutate(mean_wga_around = map(row_number(), ~ mean_wga[pmax(.x - 6, 1):pmin(.x + 6, n())])) %>%
  #mutate(ks_pval = sum(unlist(mean_diff_around))) %>%