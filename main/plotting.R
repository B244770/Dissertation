library(ggplot2)
library(reshape2)

# turn matrix to melt
design_melt <- melt(design_matrix)
replicate_melt <- melt(rep_matrix)

# Design Matrix
ggplot(design_melt, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("white", "blue"), 
                    name = "Presence",
                    labels = c("Absent", "Planted")) +
  theme_minimal() +
  labs(title = "Design Matrix",
       x = "Environment",
       y = "Genotype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Replication Matrix
ggplot(replicate_melt, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Replication Matrix",
       x = "Environment",
       y = "Genotype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
