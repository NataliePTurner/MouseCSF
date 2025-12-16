library(tidyverse)

# Create the data
data <- tibble(
  Sample = c("3", "4", "2", "1", "3", "4", "5", "1", "2"),
  Sex = c("F", "F", "M", "M", "F", "F", "F", "M", "M"),
  Condition = c("Dep", "Dep", "Dep", "Dep", "Non-dep", "Non-dep", "Non-dep", "Non-dep", "Non-dep"),
  Precursors = c(3836, 2617, 2214, 2001, 1813, 2035, 2127, 2255, 848),
  Proteins = c(572, 459, 423, 411, 357, 420, 313, 403, 159)
)

# Pivot to long format for faceting
data_long <- data %>%
  pivot_longer(cols = c(Precursors, Proteins),
               names_to = "Metric",
               values_to = "Count") %>%
  mutate(Metric = factor(Metric, levels = c("Precursors", "Proteins")))

# Calculate means per condition and metric
means <- data_long %>%
  group_by(Condition, Metric) %>%
  summarise(Mean = mean(Count), .groups = "drop")

library(ggrepel)

# Plot
ggplot(data_long, aes(x = Condition, y = Count, color = Condition, shape = Sex)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_text_repel(aes(label = Sample), position = position_dodge(width = 0.3),
                  size = 6, color = "black", max.overlaps = 20) +
  stat_summary(aes(group = 1), fun = mean, geom = "crossbar", 
               width = 0.5, color = "black", linewidth = 0.5) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_color_manual(values = c("Dep" = "#388ECC", "Non-dep" = "#F68B33")) +
  scale_shape_manual(values = c("F" = 16, "M" = 17),
                     labels = c("Female", "Male")) +
  labs(title = "Peptide Precursors and Proteins Identified in Mouse CSF",
       x = "Condition",
       y = "Count") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# With p value

library(ggpubr)

# Calculate p-values
precursor_p <- t.test(Precursors ~ Condition, data = data, var.equal = TRUE)$p.value
protein_p <- t.test(Proteins ~ Condition, data = data, var.equal = TRUE)$p.value

# Create annotation dataframe
p_labels <- data.frame(
  Metric = factor(c("Precursors", "Proteins"), levels = c("Precursors", "Proteins")),
  label = c(paste0("p = ", round(precursor_p, 3)), 
            paste0("p = ", round(protein_p, 3))),
  x = c(1.5, 1.5),
  y = c(max(data$Precursors) * 1.05, max(data$Proteins) * 1.05)
)

# Plot
ggplot(data_long, aes(x = Condition, y = Count, color = Condition, shape = Sex)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_text_repel(aes(label = Sample), position = position_dodge(width = 0.3),
                  size = 6, color = "black", max.overlaps = 20) +
  stat_summary(aes(group = 1), fun = mean, geom = "crossbar", 
               width = 0.5, color = "black", linewidth = 0.5) +
  geom_text(data = p_labels, aes(x = x, y = y, label = label), 
            inherit.aes = FALSE, size = 5) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_color_manual(values = c("Dep" = "#388ECC", "Non-dep" = "#F68B33")) +
  scale_shape_manual(values = c("F" = 16, "M" = 17),
                     labels = c("Female", "Male")) +
  labs(title = "Peptide Precursors and Proteins Identified in Mouse CSF",
       x = "Condition",
       y = "Count") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))