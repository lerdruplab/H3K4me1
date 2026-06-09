library('tidyverse')
library('ggplot2')
library('pcaMethods')
library('ggrepel')
library('corrplot')

matrix_1 <- read.table("R projects/Mixture modeling/Fig3_chip_big_matrix_10kb_263389_rawRC.txt",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = TRUE)

matrix_2 <- read.table("R projects/Mixture modeling/downsampledP10_10kbpmatrix_rawReadCounts.txt",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = TRUE)

matrix_3 <- read.table("R projects/Mixture modeling/2026_Matrix_encode_k4me1_rawReadCounts.txt",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = TRUE)

collected_df <- merge(matrix_1, matrix_2, by = "Bins")
collected_df <- merge(collected_df, matrix_3, by = "Bins")

rownames(collected_df) <- collected_df$Bins
collected_df$Bins <- NULL

names(collected_df)

#loading in the data from the right matrices to be used for analysis
df_pca <- collected_df %>% select(CD1_p10_rep1=p10A_H3K4me1_rep1,
                                  CD1_p10_rep2=p10B_H3K4me1_rep2,
                                  BL6_p10_rep1=P10_k4me1_rep1,
                                  BL6_p10_rep2=P10_k4me1_rep2,
                                  ES_k4me1,
                                  Kidney_k4me1,
                                  Liver_k4me1,        
                                  SmallIntestine_k4me1,
                                  Testis_k4me1)

normalize_data <- function(data) {
  index <- 0.05 * nrow(data)
  
  normalized_data <- apply(data, 2, function(col) {
    sorted_col <- sort(col, decreasing = TRUE)
    normalization_value <- sorted_col[index]
    col / normalization_value
  })
  
  return(as.data.frame(normalized_data))
}

df_pca_normalized <- normalize_data(df_pca)

pcoutput <- pca(t(df_pca_normalized))
dat <- as.data.frame(scores(pcoutput)[,1:2])
dat$celltype <- rownames(dat)

#for the percentages
pcoutput@R2

#percentages taken from results above
p <- ggplot(dat, aes(PC1, PC2, colour=celltype, label=celltype)) + 
  geom_point(size=3) + 
  geom_text_repel() + 
  theme_minimal() +
  labs(x = "PC1 66%", y = "PC2 12%")

p

corrplot(cor(df_pca_normalized, df_pca_normalized, method = "spearman"), method = "color", addCoef.col = "black", tl.col = "black", tl.srt = 45)  

