library(lme4)
library(ggplot2)

# 假设pheno_df包含了列'id'（基因型标识）和'y.Trait1'（表型值）
# 构建模型，其中基因型作为随机效应
model <- lmer(y.Trait1 ~ 1 + (1 | id), data = pheno_df)

blups <- ranef(model)$id


# generate data frame for blup
blup_df <- data.frame(id = rownames(blups), blup = unlist(blups[,"(Intercept)"]))

# plot blup
# 确保id为数值类型并正确排序
blup_df$id <- factor(blup_df$id, levels = sort(unique(as.numeric(as.character(blup_df$id)))))

# 绘制散点图
ggplot(blup_df, aes(x = id, y = blup)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = c("1", "50", "100", "150", "200")) +  # 只显示特定的几个标签
  labs(title = "BLUP Values for Each Genotype",
       x = "Genotype ID",
       y = "BLUP Value")

