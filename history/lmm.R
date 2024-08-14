library(lme4)
# 使用混合模型分析GxE互作
analyze_interaction <- function(data) {
  # 模型包括随机效应
  lmm_result <- lmer(Accuracy ~ nGenos * nEnvs + (1|ID), data = data)
  return(summary(lmm_result))
}

# 进行GxE互作分析
interaction_results <- analyze_interaction(results_df)
print(interaction_results)
