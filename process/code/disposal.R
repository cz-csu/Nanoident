
library(XML)

# 定义一个函数来合并相似的模体
merge_similar_motifs <- function(df, threshold) {
  i <- 1
  while (i <= nrow(df)) {
    # 找到与当前模体相似的模体
    similar_indices <- agrep(df$motif[i], df$motif, max.distance = threshold)
    
    # 如果找到了相似的模体
    if (length(similar_indices) > 1) {
      # 合并相似的模体
      df$motif[i] <- paste(df$motif[similar_indices], collapse = "/")
      df$log_p_value[i] <- min(df$log_p_value[similar_indices])
      
      # 去除相似的模体的行，除了当前行
      similar_indices <- similar_indices[similar_indices != i]
      df <- df[-similar_indices, ]
    }
    
    i <- i + 1
  }
  
  return(df)
}

doc <- xmlParse("analysis/NO/motif_detection/meme_NO/dreme_1/dreme.xml")

motif_nodes <- getNodeSet(doc, "//motif")
motifs <- sapply(motif_nodes, function(node) {xmlGetAttr(node, "seq")})
p_value <- sapply(motif_nodes, function(node) {xmlGetAttr(node, "pvalue")})

dreme_result <- data.frame(
  motif = motifs,
  log_p_value = -log10(as.numeric(p_value))
)

#print(nrow(dreme_result))
#print(dreme_result)

# 读取文件
data <- read.csv("analysis/NO/motif_detection/meme_NO/homer_1/knownResults.txt", header = TRUE, sep = "\t",fileEncoding="utf-8")
motifs <- data$Consensus
log_p_value <-data$Log.P.value
homer_old_result <- data.frame(
  motif = motifs,
  log_p_value = as.numeric(log_p_value)
)
homer_old_result <- homer_old_result[homer_old_result$log_p_value < -30, ]
#print(homer_old_result)

# 读取文件
lines <- readLines("analysis/NO/motif_detection/meme_NO/homer_1/homerMotifs.all.motifs")
motifs <- c() 
log_p_values <- c()
# 遍历每一行
for (line in lines) {
  # 如果这一行是一个模体的描述信息
  if (startsWith(line, ">")) {
    # 提取模体的名称
    name <- strsplit(line, "\t")[[1]][1]
    motif <- strsplit(name, ">")[[1]][2]
    log_p_value<- strsplit(line, "\t")[[1]][4]
    motifs <- c(motifs, motif)
    log_p_values <- c(log_p_values, log_p_value)
  } 
}
homer_new_result <- data.frame(
  motif = motifs,
  log_p_value = as.numeric(log_p_values)
)
homer_new_result <- homer_new_result[homer_new_result$log_p_value < -300, ]
#print(homer_new_result)


final_result <- rbind(dreme_result, homer_old_result, homer_new_result)
final_result$motif=as.character(final_result$motif)
print(final_result)
df <- merge_similar_motifs(final_result, threshold = 0.2)
print(df)
#print(class(final_result$motif))
print(dreme_result)