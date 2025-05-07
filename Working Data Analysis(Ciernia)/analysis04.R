#Microglia - neuron stat testing (not used)

library(MicrogliaMorphologyR)
library(factoextra)
library(ppclust)
library(dplyr)
library(ggsignif)
library(ggrepel)
library(performance)  # for model checks
library(lme4)
library(lmerTest) # For p-values
library(MASS)
library(writexl)
calculate_cluster_percentage <- function(data, cluster_var, first_var = NULL, second_var = NULL, third_var=NULL, fourth_var=NULL) {
  
  # Dynamically determine grouping variables based on user input
  grouping_vars <- c()
  if (!is.null(first_var)) grouping_vars <- c(grouping_vars, first_var)
  if (!is.null(second_var)) grouping_vars <- c(grouping_vars, second_var)
  if (!is.null(third_var)) grouping_vars <- c(grouping_vars, third_var)
  if (!is.null(fourth_var)) grouping_vars <- c(grouping_vars, fourth_var)
  
  if (length(grouping_vars) == 0) {
    stop("You must specify at least one grouping variable (brain_region or mouse_id).")
  }
  
  # Group by the selected variables
  data %>%
    group_by(across(all_of(grouping_vars))) %>%  # Dynamic grouping
    mutate(total_within_group = sum(n, na.rm = TRUE)) %>%  # Sum only within each group
    group_by(across(all_of(grouping_vars)), !!sym(cluster_var)) %>%
    summarise(
      cluster_count = sum(n, na.rm = TRUE),  # Sum counts for the cluster
      total_within_group = first(total_within_group),  # Keep the total count per group
      output = (cluster_count / total_within_group),  
      .groups = "drop"
    ) %>%
    arrange(across(all_of(grouping_vars)), desc(output))
}


interactions_m <- read_excel("/Users/alexlawson/Desktop/Masters-Work/microglia-neurons/males/males-stats.xlsx")
interactions_f <- read_excel("/Users/alexlawson/Desktop/Masters-Work/microglia-neurons/females/females-stats.xlsx")
interactions_combined <- rbind(interactions_m, interactions_f)

cp_treatment_touchdegree_mouseid <- calculate_cluster_percentage(interactions_combined, "Cluster", "Treatment", "TouchDegree", "MouseID", "sex")
cp_treatment_mouseid <- calculate_cluster_percentage(interactions_combined, "Cluster", "Treatment", "MouseID", "sex")
cp_treatment_nocluster <- calculate_cluster_percentage(interactions_combined, "Treatment", "MouseID", "sex")
cp_treatment_nocluster_touchdegree <- calculate_cluster_percentage(interactions_combined, "Treatment", "MouseID", "TouchDegree", "sex")

View(cp_treatment_nocluster_touchdegree)
#----MALES----
#Overall-Counts
t_test_all_m <- t.test(cluster_count ~ Treatment, 
                     data = cp_treatment_nocluster %>% filter(sex=="M"),
                     var.equal = FALSE)  # Welch's t-test (recommended)
t_test_all_m
#By TD
#One
t_test_td1_m <- t.test(cluster_count ~ Treatment, 
                     data = cp_treatment_nocluster_touchdegree   %>% filter(TouchDegree=="1"&sex=="M"),
                     var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td1_m)
#Two
t_test_td2_m <- t.test(cluster_count ~ Treatment, 
                     data = cp_treatment_nocluster_touchdegree   %>% filter(TouchDegree=="2"&sex=="M"),
                     var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td2_m)

#ANOVA by cluster
View(cp_treatment_mouseid)
cp_treatment_mouseid$Treatment <- as.factor(cp_treatment_mouseid$Treatment)
cp_treatment_mouseid$Cluster <- as.factor(cp_treatment_mouseid$Cluster)
anova_model_m <- aov(cluster_count ~ Treatment * Cluster, data = cp_treatment_mouseid %>% filter(sex=="M"))
summary(anova_model_m)
tukey_results_m <- TukeyHSD(anova_model_m)
tukey_results_m
emm_m <- emmeans(anova_model_m, ~ Treatment | Cluster)
contrast(emm_m, method = "pairwise", adjust = "tukey")

# Compute estimated marginal means (least-square means)
em_results <- emmeans(anova_model_m, pairwise ~ Treatment * Cluster, adjust = "sidak")
# Print results
em_results


#ANOVA by cluster + td
cp_treatment_touchdegree_mouseid$Treatment <- as.factor(cp_treatment_touchdegree_mouseid$Treatment)
cp_treatment_touchdegree_mouseid$Cluster <- as.factor(cp_treatment_touchdegree_mouseid$Cluster)
cp_treatment_touchdegree_mouseid$TouchDegree <- as.factor(cp_treatment_touchdegree_mouseid$TouchDegree)
anova_model_td_m <- aov(cluster_count ~ Treatment * Cluster * TouchDegree, data = cp_treatment_touchdegree_mouseid %>% filter(sex=="M"))
summary(anova_model_td_m)
tukey_results_td_m <- TukeyHSD(anova_model_td_m)
print(tukey_results_td_m)
emm_td_m <- emmeans(anova_model_td_m, ~ Treatment | Cluster * TouchDegree)
contrast(emm_td_m, method = "pairwise", adjust = "tukey")


#----FEMALES----
#Overall-Counts
t_test_all_f <- t.test(cluster_count ~ Treatment, 
                       data = cp_treatment_nocluster %>% filter(sex=="F"),
                       var.equal = FALSE)  # Welch's t-test (recommended)
t_test_all_f
#By TD
#One
t_test_td1_f <- t.test(cluster_count ~ Treatment, 
                       data = cp_treatment_nocluster_touchdegree   %>% filter(TouchDegree=="1"&sex=="F"),
                       var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td1_f)
#Two
t_test_td2_f <- t.test(cluster_count ~ Treatment, 
                       data = cp_treatment_nocluster_touchdegree   %>% filter(TouchDegree=="2"&sex=="F"),
                       var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td2_f)

#ANOVA by cluster
cp_treatment_mouseid$Treatment <- as.factor(cp_treatment_mouseid$Treatment)
cp_treatment_mouseid$Cluster <- as.factor(cp_treatment_mouseid$Cluster)
anova_model_f <- aov(cluster_count ~ Treatment * Cluster, data = cp_treatment_mouseid %>% filter(sex=="F"))
summary(anova_model_f)
emm_f <- emmeans(anova_model_f, ~ Treatment | Cluster)
contrast(emm_f, method = "pairwise", adjust = "sidak")
tukey_results_f <- TukeyHSD(anova_model_f)
print(tukey_results_f)


#ANOVA by cluster + td
cp_treatment_touchdegree_mouseid$Treatment <- as.factor(cp_treatment_touchdegree_mouseid$Treatment)
cp_treatment_touchdegree_mouseid$Cluster <- as.factor(cp_treatment_touchdegree_mouseid$Cluster)
cp_treatment_touchdegree_mouseid$TouchDegree <- as.factor(cp_treatment_touchdegree_mouseid$TouchDegree)
anova_model_td_f <- aov(cluster_count ~ Treatment * Cluster * TouchDegree, data = cp_treatment_touchdegree_mouseid %>% filter(sex=="F"))
summary(anova_model_td_f)
emm_td_f <- emmeans(anova_model_td_f, ~ Treatment | Cluster * TouchDegree)
contrast(emm_td_f, method = "pairwise", adjust = "sidak")
tukey_results_td_f <- TukeyHSD(anova_model_td_f)
print(tukey_results_td_f)





#----PLOTS----

# totals M
total_touch_number_plot_m <- ggplot(cp_treatment_nocluster %>% filter(sex=="M"), aes(x = Treatment, y = cluster_count, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Transparent boxplot with no outliers shown
  geom_jitter(aes(color = Treatment), width = 0.15, size = 2, alpha = 0.8) +  # Overlay jittered data points
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Boxplot fill colors
  scale_color_manual(values = c("COLD" = "lightblue", "CON" = "lightblue")) +  # Jitter point colors
  ggtitle("Comparison of Total Number of Interactions in Cold and Control Groups in the PVN of Males") +
  labs(x = "Treatment", y = "Count", fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 14) 
total_touch_number_plot_m

# totals M, degree faceted
total_touch_number_plot_m_facet <- ggplot(cp_treatment_nocluster_touchdegree %>% filter(sex=="M"), aes(x = Treatment, y = cluster_count, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Transparent boxplot with no outliers shown
  geom_jitter(aes(color = Treatment), width = 0.15, size = 2, alpha = 0.8) +  # Overlay jittered data points
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Boxplot fill colors
  scale_color_manual(values = c("COLD" = "lightblue", "CON" = "lightblue")) +  # Jitter point colors
  ggtitle("Comparison of Total Number of Interactions in Cold and Control Groups in the PVN of Males") +
  labs(x = "Treatment", y = "Count", fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 14) +
  facet_wrap(~TouchDegree)  # Separate by TouchDegree
total_touch_number_plot_m_facet

# Display the plot
touch_number_plot
#count with degree colored
touch_numbers_count_m <- cp_treatment_touchdegree_mouseid %>% filter(sex=="M") %>%
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = as.factor(TouchDegree)), 
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), 
              size = 2, alpha = 0.6) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("1" = "lightblue", "2" = "darkblue")) +  # Custom dot colors
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN of Males: Colored by Degree") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Touch Degree") +
  theme_bw(base_size = 14)
touch_numbers_count_m

#count with degree faceted 
touch_numbers_faceted_count_m <- cp_treatment_touchdegree_mouseid %>% filter(sex=="M") %>%
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment, TouchDegree))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = as.factor(TouchDegree)), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 2, alpha = 0.6) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("1" = "lightblue", "2" = "darkblue")) +  # Custom dot colors
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN of Males: Separated by Degree") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Touch Degree") +
  theme_bw(base_size = 14) +
  facet_wrap(~TouchDegree)  # Separate by TouchDegree
touch_numbers_faceted_count_m

#count not separated by degree
touch_number_nodegree_count_m <- cp_treatment_mouseid %>% filter(sex=="M") %>%
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = Treatment),  # Remove TouchDegree, color by Treatment instead
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
              size = 2, alpha = 0.8) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("COLD" = "lightblue", "CON" = "lightblue")) +  # Jitter colors match Treatment
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN of Males") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 14)
touch_number_nodegree_count_m

# totals F
total_touch_number_plot_f <- ggplot(cp_treatment_nocluster %>% filter(sex=="F"), aes(x = Treatment, y = cluster_count, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Transparent boxplot with no outliers shown
  geom_jitter(aes(color = Treatment), width = 0.15, size = 2, alpha = 0.8) +  # Overlay jittered data points
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Boxplot fill colors
  scale_color_manual(values = c("COLD" = "lightpink", "CON" = "lightpink")) +  # Jitter point colors
  ggtitle("Comparison of Total Number of Interactions in Cold and Control Groups in the PVN of Females") +
  labs(x = "Treatment", y = "Count", fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 14) 
total_touch_number_plot_f

# totals F, degree faceted
total_touch_number_plot_f_facet <- ggplot(cp_treatment_nocluster_touchdegree %>% filter(sex=="F"), aes(x = Treatment, y = cluster_count, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Transparent boxplot with no outliers shown
  geom_jitter(aes(color = Treatment), width = 0.15, size = 2, alpha = 0.8) +  # Overlay jittered data points
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Boxplot fill colors
  scale_color_manual(values = c("COLD" = "lightpink", "CON" = "lightpink")) +  # Jitter point colors
  ggtitle("Comparison of Total Number of Interactions in Cold and Control Groups in the PVN of Females") +
  labs(x = "Treatment", y = "Count", fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 14) +
  facet_wrap(~TouchDegree)  # Separate by TouchDegree
total_touch_number_plot_f_facet

touch_numbers_count_f <- cp_treatment_touchdegree_mouseid %>% filter(sex=="F") %>%
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = as.factor(TouchDegree)), 
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), 
              size = 2, alpha = 0.6) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("1" = "#FFB6C1", "2" = "#FF1493")) +  # Custom dot colors
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN of Females: Colored by Degree") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Touch Degree") +
  theme_bw(base_size = 14)
touch_numbers_count_f

#count with degree faceted 
touch_numbers_faceted_count_f <- cp_treatment_touchdegree_mouseid %>% filter(sex=="F") %>%
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment, TouchDegree))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = as.factor(TouchDegree)), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 2, alpha = 0.6) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("1" = "#FFB6C1", "2" = "#FF1493")) +  # Custom dot colors
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN of Females: Separated by Degree") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Touch Degree") +
  theme_bw(base_size = 14) +
  facet_wrap(~TouchDegree)  # Separate by TouchDegree
touch_numbers_faceted_count_f

#count not separated by degree
touch_number_nodegree_count_f <- cp_treatment_mouseid %>% filter(sex=="F") %>%
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = Treatment),  # Remove TouchDegree, color by Treatment instead
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
              size = 2, alpha = 0.8) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("COLD" = "#FFB6C1", "CON" = "#FFB6C1")) +  # Jitter colors match Treatment
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN of Females") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 14)
touch_number_nodegree_count_f

#everything COLDvCON no facets 
touch_numbers_combined <- cp_treatment_touchdegree_mouseid %>%
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment))) +
  geom_boxplot(aes(fill = Treatment), 
               outlier.shape = NA) +
  geom_jitter(aes(color = interaction(sex, TouchDegree)), 
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), 
              size = 2, alpha = 0.6) +  # Align jitter points with boxes
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Box colors
  scale_color_manual(values = c("F.1" = "lightpink", 
                                "F.2" = "#FF1493",  # Dark pink
                                "M.1" = "lightblue", 
                                "M.2" = "darkblue")) +  # Jitter dot colors
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Sex + Touch Degree") +
  theme_bw(base_size = 14)

touch_numbers_combined

#everything COLDvCON no facets 
touch_numbers_combined <- cp_treatment_mouseid %>%
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment))) +
  geom_boxplot(aes(fill = Treatment), 
               outlier.shape = NA) +
  geom_jitter(aes(color = sex), 
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
              size = 2, alpha = 0.6) +  # Align jitter points with boxes
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Box colors
  scale_color_manual(values = c("F" = "#FF1493",  # Dark pink
                                "M" = "darkblue")) +  # Jitter dot colors
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Sex + Touch Degree") +
  theme_bw(base_size = 14)

touch_numbers_combined

#everything faceted
touch_numbers_faceted_count_all <- cp_treatment_touchdegree_mouseid %>%
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment, TouchDegree, sex))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = interaction(sex, TouchDegree)), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 2, alpha = 0.6) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("F.1" = "lightpink", 
                                "F.2" = "#FF1493", 
                                "M.1" = "lightblue", 
                                "M.2" = "darkblue")) +  # Custom dot colors
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN: Separated by Degree and Sex") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Sex + Touch Degree") +
  theme_bw(base_size = 14) +
  facet_grid(sex ~ TouchDegree)  # Facet by both sex and TouchDegree

touch_numbers_faceted_count_all


