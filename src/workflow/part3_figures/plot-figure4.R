library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(here)
library(janitor)

rm(list = ls())
.wd <- getwd()
.datPF <- file.path(.wd, "out/osf/final")
.outPF <- file.path(.wd, "out/osf/final")

### select species to include in figure 4

# create function to data wrangle model results
read_model_results <- function(file_name, response){
  results <- read_csv(file.path(.datPF, file_name)) %>%
    mutate(response = rep(response, nrow(.))) %>%
    select(species, Estimate, LCL, HCL, sig_code, response)
}

# read in model results
area_ghm <- read_model_results("area_ghm_effects_2025-06-24.csv", "area_ghm")
area_sg <- read_model_results("area_sg_effects_2025-06-24.csv", "area_sg")
niche_ghm <- read_model_results("niche_ghm_effects_2025-06-24.csv", "niche_ghm")
niche_sg <- read_model_results("niche_sg_effects_2025-06-24.csv", "niche_sg")

### prediction results
species_list <- read_csv(file.path(.wd, "out/species_list.csv"))

## area results
pred_dat <- read_csv(file.path(.datPF, "area_change_prediction_2025-12-18.csv"))

# create data frame of prediction results
spl <- unique(pred_dat$species)

diff_out <- list()

i <- 1
for(i in 1:length(spl)){
  sp_dat <- pred_dat %>% 
    filter(species == spl[i])
  
  est_low <- sp_dat %>% 
    filter(ghm_case == "low" & sg_case == "low") %>% 
    pull(est_unscaled_exp)
  
  est_high <- sp_dat %>% 
    filter(ghm_case == "high" & sg_case == "high") %>% 
    pull(est_unscaled_exp)
  
  diff <- est_low-est_high
  
  tmp_out <- tibble(species = spl[i],
                    est_low = est_low,
                    est_high = est_high,
                    diff = diff,
                    model = sp_dat$model[1],
                    tot_sig = sp_dat$tot_sig[1])
  
  diff_out[[i]] <- tmp_out
}

area_diff_df <- do.call("rbind", diff_out) %>% 
  mutate(diff_km = -diff/1000000,
         prop = est_high/est_low,
         perc_num = -round((1-prop)*100, 0),
         percent_change = ((est_high-est_low)/est_low)*100) 

area_diff_non_sig <- area_diff_df %>%
  filter(tot_sig == "non-sig")

area_diff_df <- area_diff_df %>% 
  filter(tot_sig == "sig") %>% 
  left_join(., species_list, by = c("species" = "scientific_name")) %>%
  filter(species != "Numenius americanus")

## niche results
pred_dat <- read_csv(file.path(.datPF, "niche_change_prediction_2025-12-18.csv"))

spl <- unique(pred_dat$species)

diff_out <- list()
for(i in 1:length(spl)){
  sp_dat <- pred_dat %>% 
    filter(species == spl[i])
  
  est_low <- sp_dat %>% 
    filter(ghm_case == "low" & sg_case == "low") %>% 
    pull(est_unscaled_exp)
  
  est_high <- sp_dat %>% 
    filter(ghm_case == "high" & sg_case == "high") %>% 
    pull(est_unscaled_exp)
  
  diff <- est_low-est_high
  
  tmp_out <- tibble(species = spl[i],
                    est_low = est_low,
                    est_high = est_high,
                    diff = diff,
                    model = sp_dat$model[1],
                    tot_sig = sp_dat$tot_sig[1])
  
  diff_out[[i]] <- tmp_out
}

niche_diff_df <- do.call("rbind", diff_out) %>% 
  mutate(percent_change = ((est_high-est_low)/est_low)*100) 

niche_diff_non_sig <- niche_diff_df %>%
  filter(tot_sig == "non-sig")

niche_diff_df <- niche_diff_df %>% 
  filter(tot_sig == "sig") %>% 
  left_join(., species_list, by = c("species" = "scientific_name")) %>%
  filter(species != "Numenius americanus") # keep curlew out

# bind data by order of magnitude
area_diff <- area_diff_df %>%
  select(species, common_name, taxa, percent_change) %>%
  mutate(bin = cut(percent_change, 
                   breaks = c(-100000,-10000,-1000,-100,-10,-1,-0.1,0,0.1,1,10,100,1000,10000,1000000),
                   labels = c(-6.5,-5.5,-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5))) %>%
  filter(!is.na(bin)) %>%
  group_by(bin) %>%
  arrange(taxa) %>%
  mutate(y = seq(1:n()))  %>%
  ungroup() %>%
  mutate(x = as.numeric(as.character(bin))) 


niche_diff <- niche_diff_df %>%
  select(species, common_name, taxa, percent_change) %>%
  mutate(bin = cut(percent_change, 
                   breaks = c(-100000,-10000,-1000,-100,-10,-1,-0.1,0,0.1,1,10,100,1000,10000,1000000),
                   labels = c(-6.5,-5.5,-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5))) %>%
  filter(!is.na(bin)) %>%
  group_by(bin) %>%
  arrange(taxa) %>%
  mutate(y = seq(1:n()))  %>%
  ungroup() %>%
  mutate(x = as.numeric(as.character(bin))) 


max_val <- max(c(abs(min(area_diff$x)), max(area_diff$x),abs(min(niche_diff$x)), max(niche_diff$x)))

# plot results
p1 <- ggplot(data = area_diff) +
  geom_point(aes(x = x, y = y, color = taxa), size = 2) +
  scale_fill_manual(values = c("#1481BA","#cbd081")) +
  scale_color_manual(values = c("#1481BA","#cbd081")) +
  geom_vline(aes(xintercept = 0), linetype = "solid", size = 0.5, alpha = 0.8, color = "black") +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_line(colour = "#4a4e4d", linewidth =0.3, linetype='solid'),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text = element_text(size = 7, family = "Helvetica"),
    axis.title = element_text(size = 7, family = "Helvetica"),
    axis.ticks.x = element_line(color = "#4a4e4d")) +
  scale_y_continuous(breaks = seq(0,12, by = 2), expand = expansion(mult = c(0.05, 0.05))) +  
  scale_x_continuous(breaks = seq(-6,6, by = 1),
                     labels = c("-10K", "-1K", "-100", "-10", "-1", "-0.1", "0",
                                "0.1","1", "10", "100", "1K", "10K")) +
  coord_cartesian(xlim = c(-max_val,max_val)) +
  labs(x = 'Change in area size (%)',
       y = 'Species (n)')


p2 <- ggplot(data = niche_diff) +
  geom_point(aes(x = x, y = y, color = taxa), size = 2) +
  scale_fill_manual(values = c("#1481BA","#cbd081")) +
  scale_color_manual(values = c("#1481BA","#cbd081")) +
  geom_vline(aes(xintercept = 0), linetype = "solid", size = 0.5, alpha = 0.8, color = "black") +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_line(colour = "#4a4e4d", linewidth =0.3, linetype='solid'),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text = element_text(size = 7, family = "Helvetica"),
    axis.title = element_text(size = 7, family = "Helvetica"),
    axis.ticks.x = element_line(color = "#4a4e4d")) +
  scale_y_continuous(breaks = seq(0,10, by = 2), expand = expansion(mult = c(0.15, 0.15))) +  
  scale_x_continuous(breaks = seq(-6,6, by = 1),
                     labels = c("-10K", "-1K", "-100", "-10", "-1", "-0.1", "0",
                                "0.1","1", "10", "100", "1K", "10K")) +
  coord_cartesian(xlim = c(-max_val,max_val)) +
  labs(x = 'Change in niche size (%)',
       y = 'Species (n)')

# export figure 4
p <- p1/p2 +
  plot_layout(heights = c(1.85,1.15))

ggsave(file.path(.outPF, "fig4.pdf"), height = 70, width = 120, units = "mm")

### data summaries

# number spp with non-sig changes in A & N
ns_area <- paste0("Number species with NON-SIG area change: ", n_distinct(area_diff_non_sig$species), " species")
ns_niche <- paste0("Number species with NON-SIG niche change: ", n_distinct(niche_diff_non_sig$species), " species")

# for all sig spp, range of perc change in A & N
area_percent_change_range <- paste0("Area percent change RANGE ALL SIG SPECIES: ", paste(range(area_diff_df$percent_change), collapse = " – "))
niche_percent_change_range <- paste0("Niche percent change RANGE ALL SIG SPECIES: ", paste(range(niche_diff_df$percent_change), collapse = " – "))

### --- area sig mammals --- ###
area_diff_df_mammals <- area_diff_df %>%
  filter(taxa == "mammals")
  
# area sig mammals median & range
area_percent_change_median_mammals <- paste0("Area percent change MEDIAN MAMMALS: ", median(area_diff_df_mammals$percent_change))
area_diff_km_range_mammals <- paste0("Area diff (km) RANGE MAMMALS: ", paste(range(area_diff_df_mammals$diff_km), collapse = " – "))

### --- area sig birds --- ###
area_diff_df_birds <- area_diff_df %>%
  filter(taxa == "birds") 

# area sig birds median & range
area_percent_change_median_birds <- paste0("Area percent change MEDIAN BIRDS: ", median(area_diff_df_birds$percent_change))
area_diff_km_range_birds <- paste0("Area diff (km) RANGE BIRDS: ", paste(range(area_diff_df_birds$diff_km), collapse = " – "))

# niche sig median & range & mean
niche_percent_change_median_all <- paste0("Niche percent change MEDIAN ALL SIG: ", median(niche_diff_df$percent_change))
niche_percent_change_mean_all <- paste0("Niche percent change MEAN ALL SIG: ", mean(niche_diff_df$percent_change))
niche_percent_change_range_range_all <- paste0("Niche percent change RANGE ALL SIG: ", paste(range(niche_diff_df$percent_change), collapse = " – "))

### --- niche sig mammals --- ###
niche_diff_df_mammals <- niche_diff_df %>%
  filter(taxa == "mammals")

# niche sig mammals median & range & mean
niche_percent_change_median_mammals <- paste0("Niche percent change MEDIAN MAMMALS: ", median(niche_diff_df_mammals$percent_change))
niche_percent_change_mean_mammals <- paste0("Niche percent change MEAN MAMMALS: ", mean(niche_diff_df_mammals$percent_change))
niche_percent_change_range_mammals <- paste0("Niche percent change RANGE MAMMALS: ", paste(range(niche_diff_df_mammals$percent_change), collapse = " – "))

## --- niche sig birds --- ###
niche_diff_df_birds <- niche_diff_df %>%
  filter(taxa == "birds") 

# niche sig birds median & range & mean
niche_percent_change_median_birds <- paste0("Niche percent change MEDIAN BIRDS: ", median(niche_diff_df_birds$percent_change))
niche_percent_change_mean_birds <- paste0("Niche percent change MEAN BIRDS: ",mean(niche_diff_df_birds$percent_change))
niche_change_range_range_birds <- paste0("Niche percent change RANGE BIRDS: ", paste(range(niche_diff_df_birds$percent_change), collapse = " – "))

values <- c(ns_area, 
            ns_niche,
            area_percent_change_range,
            niche_percent_change_range,
            area_percent_change_median_mammals,
            niche_percent_change_median_mammals,
            niche_percent_change_range_mammals,
            niche_percent_change_mean_mammals,
            area_percent_change_median_birds,
            niche_percent_change_median_birds,
            niche_percent_change_mean_birds,
            niche_change_range_range_birds,
            area_diff_km_range_mammals,
            area_diff_km_range_birds,
            niche_percent_change_median_all,
            niche_percent_change_mean_all,
            niche_percent_change_range_range_all)

writeLines(values, file.path(.outPF, "figure4_summary.txt"), sep = "\n", useBytes = TRUE)
