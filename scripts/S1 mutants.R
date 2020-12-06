# Load packages -----------------------------------------------------------
library(tidyverse)
library(rethinking)
library(tidybayes)
library(tidybayes.rethinking)
library(modelr)
library(patchwork)

# Read data ---------------------------------------------------------------
elisa <- read_csv("data/ELISA_data.csv")

# Summarise data ----------------------------------------------------------
elisa
summary(elisa)
glimpse(elisa)

# Convert Protein variable into factor ------------------------------------
elisa$Protein <- factor(
  elisa$Protein, 
  levels = c("S1", "AT-S1", "AT-S1+BLV", "MUTANT", "MUTANT+BLV")
  )

# Plot data ---------------------------------------------------------------
ggplot(elisa, aes(x = Protein, y = OD, group = Sample)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  theme_bw()

# Numerically encode categorical variables --------------------------------
model_data <- elisa %>%
  mutate(Protein = as.integer(Protein),
         Sample = as.integer(as.factor(Sample))) %>% 
  select(-SCO)

# Model data --------------------------------------------------------------
model_1 <- ulam(
  alist(
    # likelihood function
    OD ~ normal(mu, sigma),
    # linear model
    mu <- a[Sample] + b[Protein],
    # priors
    b[Protein] ~ normal(0, 1),
    sigma ~ log_normal(-0.5, 0.5),
    a[Sample] ~ normal(mu_sample, tau_sample),
    # hyperpriors
    mu_sample ~ normal(1, 5),
    tau_sample ~ exponential(1)
  ), 
  data = model_data, 
  chains = 4, 
  cores = 4, 
  iter = 2000,
  log_lik = TRUE
)

plot(model_1, depth = 1)
precis(extract.samples(model_1), depth = 2)
 
model_2 <- ulam(
  alist(
    # likelihood function
    OD ~ dgamma2(mu, scale),
    # linear model
    log(mu) <- a[Sample] + b[Protein],
    # priors
    b[Protein] ~ normal(0, 1),
    scale ~ exponential(5),
    a[Sample] ~ dgamma2(mu_sample, tau_sample),
    # hyperpriors
    mu_sample ~ normal(1, 0.2),
    tau_sample ~ exponential(5)
  ), 
  data = model_data, 
  chains = 4, 
  cores = 4, 
  iter = 10000,
  log_lik = TRUE
)

plot(model_2, depth = 2)
precis(extract.samples(model_2), depth = 2)

model_3 <- ulam(
  alist(
    # likelihood function
    OD ~ dgamma2(mu, scale),
    # linear model
    log(mu) <- b[Protein],
    # priors
    b[Protein] ~ normal(0, 1),
    scale ~ exponential(5)
  ), 
  data = model_data, 
  chains = 4, 
  cores = 4, 
  iter = 2000,
  log_lik = TRUE
)

plot(model_3, depth = 3)
precis(extract.samples(model_3), depth = 3)

# Comparing models --------------------------------------------------------
PSIS(model_1)
PSIS(model_2) # superior
PSIS(model_3)

# Visualize posterior ------------------------------------------------------
model_data %>%
  add_fitted_draws(model_2) %>%
  mutate(OD = .value,
         Protein = fct_recode(as.factor(Protein),
                              "S1" = "1", "AT-S1" = "2", "AT-S1+BLV" = "3", 
                              "MUTANT" = "4", "MUTANT+BLV" = "5")) %>%
  ggplot(aes(x = Protein, y = OD)) +
  stat_lineribbon(aes(y = .value), .width = seq(0.2, 0.9, 0.1)) +
  geom_point(stat = "summary", fun = "median", col = "red", aes(group = Protein)) +
  scale_fill_brewer(palette = "Blues") +
  theme_bw() 

# Parameter estimates -----------------------------------------------------
get_variables(model_2)

model_2 %>%
  recover_types(elisa) %>% # converts index variables to original variables
  gather_draws(b[Protein], mu_sample, tau_sample) %>%
  median_hdi(.width = 0.95) %>%
  mutate(.value = 10^(.value),
         .lower = 10^(.lower),
         .upper = 10^(.upper)
         )

# Plotting parameter uncertainty ------------------------------------------
model_2 %>%
  recover_types(elisa) %>% # converts index variables to original variables
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
  recover_types(elisa) %>%
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
  recover_types(elisa) %>% # converts index variables to original variables
  gather_draws(b[Protein], mu_sample, tau_sample)

median_mu_sample <- median(filter(post_draws, .variable == "mu_sample")$.value)
median_tau_sample <- median(filter(post_draws, .variable == "tau_sample")$.value)

gamma_sim <- tibble(x = seq(0, 3.5, 0.01),
                    y = dgamma2(x, median_mu_sample, 1.5)
)

ggplot(gamma_sim, aes(x, 200*y)) + 
  geom_histogram(data = elisa, aes(OD), inherit.aes = FALSE) +
  geom_line(col = "blue") +
  theme_bw()

# Additional plots --------------------------------------------------------
plot_a <- elisa %>%
  filter(Protein == "AT-S1" | Protein == "AT-S1+BLV" | 
         Protein == "MUTANT"|  Protein == "MUTANT+BLV") %>%
  pivot_wider(names_from = Protein, values_from = OD, id_cols = Sample) %>%
  pivot_longer(names_to = "WT",
               values_to = "OD-BLV",
               cols = c(`AT-S1`, `MUTANT`)) %>%
  pivot_longer(names_to = "N121Q",
               values_to = "OD+BLV",
               cols = c(`AT-S1+BLV`, `MUTANT+BLV`)) %>%
  filter((WT == "AT-S1" & N121Q == "AT-S1+BLV") |
           (WT == "MUTANT" & N121Q == "MUTANT+BLV")) %>%
  ggplot(aes(`OD-BLV`, `OD+BLV`, col = WT)) +
  geom_point() +
  labs(x = "OD, -biliverdin", y = "OD, +biliverdin") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  annotate(geom = "text", label = "WT", x = 2.8, y = 2.2, size = 5, col = "#E41A1C") +
  annotate(geom = "text", label = "N121Q", x = 2.2, y = 2.8, size = 5, col = "#377EB8") +
  scale_color_brewer(type = "qual", palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none")
  
plot_b <- model_2 %>%
  recover_types(elisa) %>%
  spread_draws(b[Protein]) %>%
  compare_levels(b, by = Protein) %>%
  ungroup() %>%
  filter(Protein == "AT-S1+BLV - AT-S1" | Protein == "MUTANT+BLV - MUTANT") %>%
  ggplot(aes(y = 1, x = b, fill = Protein)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  stat_halfeye(.width = c(0.85, 0.95), point_interval = "median_hdi") + 
  annotate(geom = "text", label = "WT", x = -0.114, y = 2.2, size = 3) +
  annotate(geom = "text", label = "N121Q", x = 0.021, y = 2.2, size = 3) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_y_discrete(labels = c("WT", "N121Q")) +
  labs(y = "", x = "\u0394OD (+/- biliverdin)") +
  theme_bw() +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.grid.major.y = element_blank())

plot_a + inset_element(plot_b, left = -0.02, bottom = 0.6, right = 0.5, top = 0.99)
ggsave("plots/Scatter plot with inlay.pdf", width = 8, height = 5, device = cairo_pdf)
