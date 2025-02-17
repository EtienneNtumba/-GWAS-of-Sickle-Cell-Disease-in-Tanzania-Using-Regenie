# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read the data
gwas_data <- read.table("/home/p0129674/Documents/Emile_Analysis/Data/resultats/qc_results/step2_P.txt", 
                        header = TRUE, 
                        stringsAsFactors = FALSE)

# Convert p-values correctly
gwas_data <- gwas_data %>% 
  mutate(P = 10^(-LOG10P)) %>% 
  filter(!is.infinite(-log10(P)))  # Remove infinite values

# Create expected -log10(P) values under null hypothesis
gwas_data <- gwas_data %>%
  arrange(P) %>%
  mutate(Expected = -log10(ppoints(n())))

# Generate QQ plot
qq_plot <- ggplot(gwas_data, aes(x = Expected, y = -log10(P))) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(
    title = "QQ Plot of GWAS p-values",
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)"
  ) +
  theme_minimal()

# Save the QQ plot as an image
ggsave("/home/p0129674/Documents/Emile_Analysis/Data/resultats/qc_results/qq_plot.png", 
       plot = qq_plot, width = 8, height = 6, dpi = 300)

# Display the plot
print(qq_plot)
