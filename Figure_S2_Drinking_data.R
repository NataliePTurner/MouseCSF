library(tidyverse)

# Create raw data (excluding M15)
drinking_raw <- tibble(
  Mouse = c(3, 5, 8, 11, 14, 16, 18, 19, 20),
  Sex = c("Male", "Male", "Male", "Male", "Female", "Female", "Female", "Female", "Female"),
  Group = c("Non-dep", "Non-dep", "Dep", "Dep", "Non-dep", "Non-dep", "Non-dep", "Dep", "Dep"),
  Baseline = c(0.75, 1.40, 0.73, 0.52, 4.24, 2.30, 3.38, 3.94, 2.38),
  CIE1 = c(1.67, 2.15, 2.45, 0.96, 3.65, 2.57, 1.76, 2.87, 4.43),
  CIE2 = c(2.93, 2.35, 2.04, 1.53, 3.65, 2.84, 4.41, 5.17, 4.57),
  CIE3 = c(2.21, 2.63, 1.91, 2.00, 2.86, 2.63, 4.01, 5.50, 3.24),
  CIE4 = c(2.63, 1.45, 4.17, 2.57, 3.60, 3.10, 3.28, 5.38, 3.38),
  CIE5 = c(2.78, 2.73, 4.21, 2.59, 3.25, 3.33, 3.17, 5.72, 4.05),
  CIE6 = c(2.39, 2.07, 3.00, 1.10, 3.25, 3.51, 1.66, 5.59, 1.86)
)

# Map mouse numbers to D/ND labels
mouse_to_label <- c("3" = "ND1", "5" = "ND2", "8" = "D1", "11" = "D2", 
                    "14" = "ND3", "16" = "ND4", "18" = "ND5", "19" = "D3", "20" = "D4")

# Pivot to long format
drinking_long <- drinking_raw %>%
  pivot_longer(cols = Baseline:CIE6,
               names_to = "Timepoint",
               values_to = "Intake") %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("Baseline", "CIE1", "CIE2", "CIE3", "CIE4", "CIE5", "CIE6")),
    Group_Sex = factor(paste(Group, Sex), 
                       levels = c("Non-dep Male", "Dep Male", "Non-dep Female", "Dep Female")),
    Animal_ID = factor(mouse_to_label[as.character(Mouse)], 
                       levels = c("ND1", "ND2", "D1", "D2", "ND3", "ND4", "ND5", "D3", "D4"))
  )

# Calculate summary statistics
drinking_summary <- drinking_long %>%
  group_by(Timepoint, Group_Sex) %>%
  summarise(
    Mean = mean(Intake),
    SEM = sd(Intake) / sqrt(n()),
    .groups = "drop"
  )

# Plot with bars and individual points
ggplot() +
  geom_bar(data = drinking_summary, 
           aes(x = Timepoint, y = Mean, fill = Group_Sex),
           stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(data = drinking_summary,
                aes(x = Timepoint, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM, group = Group_Sex),
                position = position_dodge(width = 0.8), width = 0.25) +
  geom_point(data = drinking_long,
             aes(x = Timepoint, y = Intake, group = Group_Sex, shape = Animal_ID),
             position = position_dodge(width = 0.8), size = 2.5, alpha = 0.7) +
  scale_fill_manual(values = c("Non-dep Male" = "#4472C4", 
                               "Dep Male" = "#ED7D31",
                               "Non-dep Female" = "#A5A5A5", 
                               "Dep Female" = "#FFC000")) +
  scale_shape_manual(values = c("ND1" = 0, "ND2" = 1,           # Non-dep Male
                                "D1" = 2, "D2" = 6,             # Dep Male
                                "ND3" = 15, "ND4" = 16, "ND5" = 17,  # Non-dep Female
                                "D3" = 18, "D4" = 8)) +         # Dep Female
  labs(x = "CIE-2BC Week Number",
       y = "Total EtOH Intake (g/kg)",
       fill = "",
       shape = "Animal ID") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        panel.grid.major.x = element_blank()) +
  guides(fill = guide_legend(order = 1, nrow = 1),
         shape = guide_legend(order = 2, nrow = 2))