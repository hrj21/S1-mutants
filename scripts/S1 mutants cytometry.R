# Load packages -----------------------------------------------------------
library(tidyverse)
library(rethinking)
library(tidybayes)
library(tidybayes.rethinking)
library(modelr)
library(patchwork)

# Read data ---------------------------------------------------------------
flow <- read_csv("data/Cytometry_data.csv")

# Summarise data ----------------------------------------------------------
flow
summary(flow)
glimpse(flow)

# Convert Protein variable into factor ------------------------------------
flow$Protein <- factor(
  flow$Protein, 
  levels = c("WT", "WT+BLV", "MUTANT", "MUTANT+BLV")
)

# Extract sample identifiers ----------------------------------------------
flow$Sample <- purrr::map_chr(flow$Sample, ~{
  str_split(., "BV ") %>% `[[`(1) %>% `[`(2)
})

# Plot data ---------------------------------------------------------------
ggplot(flow, aes(x = Protein, y = MFI, group = Sample)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  theme_bw()

# Numerically encode categorical variables --------------------------------
model_data <- flow %>%
  mutate(Protein = as.integer(Protein),
         Sample = as.integer(as.factor(Sample)))

# Model data --------------------------------------------------------------
model_2 <- ulam(
  alist(
    # likelihood function
    MFI ~ normal(mu, sigma),
    # linear model
    mu <- a[Sample] + b[Protein],
    # priors
    b[Protein] ~ normal(0, 400),
    sigma ~ exponential(0.1),
    a[Sample] ~ normal(mu_sample, tau_sample),
    # hyperpriors
    mu_sample ~ normal(2500, 200),
    tau_sample ~ exponential(0.01)
  ), 
  data = model_data, 
  chains = 4, 
  cores = 4, 
  iter = 10000,
  log_lik = FALSE
)

plot(model_2, depth = 2)
precis(extract.samples(model_2), depth = 2)

# Comparing models --------------------------------------------------------
PSIS(model_1)
PSIS(model_2) # superior
PSIS(model_3)

# Visualize posterior ------------------------------------------------------
model_data %>%
  add_fitted_draws(model_2) %>%
  mutate(MFI = .value,
         Protein = fct_recode(as.factor(Protein),
                              "WT" = "1", "WT+BLV" = "2", 
                              "MUTANT" = "3", "MUTANT+BLV" = "4")) %>%
  ggplot(aes(x = Protein, y = MFI)) +
  stat_lineribbon(aes(y = .value), .width = seq(0.2, 0.9, 0.1)) +
  geom_point(stat = "summary", fun = "median", col = "red", aes(group = Protein)) +
  scale_fill_brewer(palette = "Blues") +
  theme_bw()

# Parameter estimates -----------------------------------------------------
get_variables(model_2)

model_2 %>%
  recover_types(flow) %>% # converts index variables to original variables
  gather_draws(b[Protein], mu_sample, tau_sample) %>%
  median_hdi(.width = 0.95)

# Plotting parameter uncertainty ------------------------------------------
model_2 %>%
  recover_types(flow) %>% # converts index variables to original variables
  gather_draws(
    b[Protein],
    mu_sample,
    tau_sample
  ) %>%
  mutate(.variable = case_when(
    is.na(Protein) ~ .variable,
    TRUE ~ paste(.variable, Protein)
  )
  ) %>% 
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < .5))) +
  geom_vline(xintercept = c(-.5, .5), linetype = "dashed") +
  scale_fill_manual(values = c("gray80", "skyblue")) +
  theme_bw() +
  stat_pointinterval(point_interval = "median_hdi", 
                     .width = c(0.89, 0.99)) 

# Pairwise contrasts ------------------------------------------------------
model_2 %>%
  recover_types(flow) %>%
  spread_draws(b[Protein]) %>%
  compare_levels(b, by = Protein) %>%
  ungroup() %>%
  mutate(Protein = reorder(Protein, b)) %>%
  ggplot(aes(y = Protein, x = b)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  stat_halfeye(.width = c(0.79, 0.89), fill = "skyblue") + 
  scale_fill_manual(values = "skyblue") +
  labs(y = "", x = "Difference in fitted value") +
  theme_bw() +
  theme()

# Visualize distribution of Sample intercepts -----------------------------
post_draws <- model_2 %>%
  recover_types(flow) %>% # converts index variables to original variables
  gather_draws(b[Protein], mu_sample, tau_sample)

median_mu_sample <- median(filter(post_draws, .variable == "mu_sample")$.value)
median_tau_sample <- median(filter(post_draws, .variable == "tau_sample")$.value)

normal_sim <- tibble(x = seq(0, 5000, 10),
                    y = dnorm(x, median_mu_sample, median_tau_sample)
)

ggplot(normal_sim, aes(x, 2e4*y)) + 
  geom_histogram(data = flow, aes(MFI), inherit.aes = FALSE) +
  geom_line(col = "blue") +
  theme_bw()

# Additional plots --------------------------------------------------------
plot_a <- flow %>%
  pivot_wider(names_from = Protein, values_from = MFI, id_cols = Sample) %>%
  pivot_longer(names_to = "WT",
               values_to = "OD-BLV",
               cols = c(`WT`, `MUTANT`)) %>%
  pivot_longer(names_to = "N121Q",
               values_to = "OD+BLV",
               cols = c(`WT+BLV`, `MUTANT+BLV`)) %>%
  filter((WT == "WT" & N121Q == "WT+BLV") |
           (WT == "MUTANT" & N121Q == "MUTANT+BLV")) %>%
  ggplot(aes(`OD-BLV`, `OD+BLV`, col = WT)) +
  geom_point() +
  labs(x = "OD, -biliverdin", y = "OD, +biliverdin") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  annotate(geom = "text", label = "N121Q", x = 5000, y = 5800, size = 5, col = "#377EB8") +
  annotate(geom = "text", label = "WT", x = 6500, y = 4000, size = 5, col = "#E41A1C") +
  scale_color_brewer(type = "qual", palette = "Set1", direction = -1) +
  theme_bw() +
  theme(legend.position = "none")

plot_b <- model_2 %>%
  recover_types(flow) %>%
  spread_draws(b[Protein]) %>%
  compare_levels(b, by = Protein) %>%
  ungroup() %>%
  filter(Protein == "WT+BLV - WT" | Protein == "MUTANT+BLV - MUTANT") %>%
  ggplot(aes(y = 1, x = b, fill = Protein)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  stat_halfeye(.width = c(0.85, 0.95), point_interval = "median_hdi") + 
  annotate(geom = "text", label = "WT", x = -950, y = 2.2, size = 3) +
  annotate(geom = "text", label = "N121Q", x = -10, y = 2.2, size = 3) +
  scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
  scale_y_discrete(labels = c("WT", "N121Q")) +
  labs(y = "", x = "\u0394OD (+/- biliverdin)") +
  theme_bw() +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.grid.major.y = element_blank())

plot_a + inset_element(plot_b, left = -0.02, bottom = 0.6, right = 0.5, top = 0.99)
ggsave("plots/Scatter plot with inlay cytometry.pdf", width = 8, height = 5, device = cairo_pdf)
