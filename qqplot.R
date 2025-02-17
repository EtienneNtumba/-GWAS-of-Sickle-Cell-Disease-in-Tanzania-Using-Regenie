# Load necessary libraries
library(ggplot2)
library(dplyr)

# Lire les données
gwas_data <- read.table("/home/p0129674/Documents/Emile_Analysis/Data/resultats/qc_results/step2_P.txt", 
                        header = TRUE, 
                        stringsAsFactors = FALSE)

# Convertir correctement les p-values
gwas_data <- gwas_data %>% 
  mutate(P = 10^(-LOG10P)) %>% 
  filter(!is.infinite(-log10(P)))  # Supprimer les valeurs infinies

# Créer les valeurs attendues de -log10(P) sous l'hypothèse nulle
gwas_data <- gwas_data %>%
  arrange(P) %>%
  mutate(Expected = -log10(ppoints(n())))

# Générer le QQ plot avec fond blanc
qq_plot <- ggplot(gwas_data, aes(x = Expected, y = -log10(P))) +
  geom_point(alpha = 0.6, color = "blue") +  # Points en bleu
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # Ligne de référence
  labs(
    title = "QQ Plot of GWAS p-values",
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)"
  ) +
  theme_classic() +  # Fond blanc
  theme(
    panel.grid.major = element_blank(),  # Supprimer les grilles
    panel.grid.minor = element_blank()
  )

# Sauvegarder le QQ plot avec fond blanc
ggsave("/home/p0129674/Documents/Emile_Analysis/Data/resultats/qc_results/qq_plot.png", 
       plot = qq_plot, width = 8, height = 6, dpi = 300)

# Afficher le plot
print(qq_plot)
