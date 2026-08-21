#ADD mutant analysis
library(tidyverse)
library(ggplot2)

df <- read.table("data_for_R/WTandMUT_10kbp_wk4me1_signal.txt", sep = '\t', header = TRUE)

df_correlation <- df %>%
  rename(
    WTmet = WT,
    MUTmet = MUT) %>%
  filter(WTmet > 1 | MUTmet > 1)  %>%
  mutate(log2fold = log2(MUTmet/WTmet)) %>%
  filter(log2fold != Inf & log2fold != -Inf & NSN.k4me1 != 0 & NSN.k4me1 <= 1) %>%
  na.omit()

# histogram of the log2fold change with added separator
hist(df_correlation$log2fold, breaks = 300)
abline(v=0, col = 'red')

# defining the different log2fold "groups"
df_correlation <- df_correlation %>% mutate(foldlvl = if_else(log2fold > 0, "plus", "minus"))
df_correlation <- df_correlation %>% mutate(foldlvlbinary = if_else(log2fold > 0, 1, 0))

#plotting the log2fold groups on the original plot and as boxplots based on NSN k4me1 level
ggplot(df_correlation, aes(x = WTmet, y = MUTmet, color = factor(foldlvl))) +
  geom_point(size = 0.4) +
  scale_color_manual(values = c("lightsalmon2", "lightskyblue2")) +
  theme_classic()

ggplot(df_correlation, aes(x = foldlvl, y = NSN.k4me1, fill = foldlvl)) +
  geom_violin() +
  geom_boxplot(width=0.15, color="black", alpha=0.1, linewidth = 0.6) +
  labs(title = "Methylation change from WT to ADD and its equivalent k4me1 lvl in NSN", x = "log2fold change mutant/wt", y = "NSN k4me1 level") +
  scale_fill_manual(values = c("lightsalmon2", "lightskyblue2")) +
  theme_minimal()

# n from the different categories
sum(df_correlation$foldlvl == "plus", na.rm = TRUE)
sum(df_correlation$foldlvl == "minus", na.rm = TRUE)
