#!/usr/bin/env Rscript 
# Plot model summary sheet (PDF) for area models. Select either the additive
# or interactive model based on significance. 


#---- Input Parameters ----#
if(interactive()) {
  
  .wd <- getwd()
  .datP <- file.path(.wd,'out/single_species_models_final')

} else {

  library(docopt)

  .wd <- getwd()
  source('src/funs/input_parse.r')
  .datP <- file.path(.wd,'out/single_species_models')
  .dbPF <- '/tmp/mosey_mod.db'
}


#---- Initialize Environment ----#
t0 <- Sys.time()

source('src/startup.r')

suppressWarnings(
  suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(ggthemes)
    library(patchwork)
    library(brms)
    library(grid)
    library(emmeans)
    library(parameters)
    library(gridExtra)
    library(sf)
    library(ggsflabel)
    library(janitor)
    library(glue)
    library(assertthat)
  }))

#Source all files in the auto load funs directory
list.files('src/funs/auto',full.names=TRUE) %>%
  walk(source)

palnew <- c("#F98177", "#8895BF")

#---- Load data ----#
message("***LOADING DATA***")

#-- Interaction Models --#
message('Loading interaction models...')
int_modlist <- list.files(path=file.path(.datP, "area_interactive"), full.names = F )
int_modlist_full <- list.files( path=file.path(.datP, "area_interactive"), full.names = T)
int_sp <- word(int_modlist, 1, sep = "_")

#-- Additive Models --#
message('Loading additive‚ models...')
add_modlist <- list.files(path=file.path(.datP, "area_additive"), 
                          full.names = F)
add_modlist_full <- list.files(path=file.path(.datP, "area_additive"), 
                               full.names = T)
add_sp <- word(add_modlist, 1, sep = "_")

# check that lists are same
int_sp == add_sp

# Init lists to store results
pred_out <- list()
sg_effects_out <- list()
ghm_effects_out <- list()

# define figure number sequence starting at 7
fig_nums <- 7:(7 + length(int_modlist_full) - 1)

# Loop over models
for(i in 1:length(int_modlist_full)){
  
  # Interactive Model
  if(int_modlist_full[i] != "NULL"){
    load(int_modlist_full[i]) # load model
    intmod <- out$model # store as object
    #sp <- out$species # extract sp
    
    fe <- parameters(intmod) #get fixed effects
    re <- posterior_summary(intmod, variable = c("sd_grp__Intercept", "sigma")) # get random effects
    
    if(fe$pd[fe$Parameter=="b_sg_norm:ghm_scale"] < 0.95 & fe$pd[fe$Parameter=="b_sg_norm:ghm_scale"] > 0.05){ # if the interacton effect is non-sig...
      
      #... load the additive model instead.
      if(add_modlist_full[i] != "NULL"){
        
        #-- ADDITIVE --#
        
        #- Model Basics -#        
        load(add_modlist_full[i]) # load model
        addmod <- out$model
        sp <- out$species

        # correct species names
        if (sp == "Accipiter gentilis"){
          sp <- "Astur atricapillus"
          }

        if (sp == "Spilogale putorius"){
          sp <- "Spilogale interrupta"
          }

        if (sp == "Cervus elaphus"){
          sp <- "Cervus canadensis"
          }

        out_add <- out
        fe_add <- fixef(out_add$model) #get fixed effects
        #fe_add <- parameters(out_add$model)
        mod <- "additive"
        variable <- "area size"

        caption <- glue("Figure S{fig_nums[i]}. Model diagnostics plots for {sp} {mod} {variable} models. A) Posterior predictive distribution (y rep)\n
                       compared to {variable} data distribution. B) Posterior predictive errors plot, showing error (y-y_rep) on x axis relative to \n
                       {variable} data distribution (y) on y axis. C: Trace plots showing  individual MCMC draws to evaluate convergence among\n
                       chains for each parameter in the model ('b_intercept' = global intercept, 'b_sg_norm' = effect of human mobility,\n
                       'b_ghm_scale' = effect of human modification, 'b_ndvi_scale' = effect of ndvi, 'b_tmax_scale' = effect of maximum daily \n
                       temperature, 'sd_grp_Intercept' = standard deviation of random intercept by individual-year, 'nu' = Student's t distribution\n
                       degrees of freedom parameter, 'ar[1]' = first-order autoregressive coefficient, 'b_sg_norm:ghm_scale' = interactive effect of\n
                       human modification and mobility [only included when interaction model selected]). D) Table reporting median and 95% credible\n
                       intervals for the posterior distributions, probability of direction ('pd'), Gelman-Rubin diagnostic ('Rhat'), and effective sample\n
                       size ('ESS'), for each parameter (names follow those described for panel C).")

        # convert into a ggplot object so it's compatible with
        caption_plot <- ggplot() +
                        geom_text(aes(-0.9, 0, label = caption), size = 3, hjust = 0) +
                        theme_void() + xlim(-1, 1)
        
        #- Extract parameter table as grob -#
        parameters_df <- parameters(addmod) %>% select(-c(CI, Component))
        param_table <- tableGrob(parameters_df)

        #- Posterior Predictive Plot -#
        (pp_dens <- pp_check(addmod)+ggtitle("Posterior Predictive Distribution"))
        
        #- Predictive Error Plot -#
        (pp_err <- pp_check(addmod, type='error_scatter_avg')+ggtitle("Posterior Predictive Errors"))  
      
        #- MCMC Trace Plot -#
        (trace <- mcmc_plot(addmod, type = "trace") +ggtitle("MCMC Traces"))    

        #---- Assemble Plots ---#
        
        # A = param table
        # B/C = ppdens & pperr
        # D = trace plots (need 2 rows so not vertically squished)
        # E = caption
        
        # Design layout
        design <- "AAAAAA
                   BBBCCC
                   DDDDDD
                   DDDDDD
                   EEEEEE
                   ######
                   ######"

        # Gather plots
        (model_out <- wrap_elements(wrap_elements(full = param_table) + ggtitle("Model Coefficient Table") +
                                    pp_dens+pp_err+
                                    trace+
                                    caption_plot+
                                    plot_layout(design = design)))

        # Write out plot
        ggsave(model_out, filename = glue(.wd, "/out/model_diagnostics_final/area/{sp}.pdf"), 
               width = 10, height = 20, device = cairo_pdf)
        
      } # if add model is not NULL
    } else {
      
      #-- INTERACTIVE --#
      #- Model Basics -#        
      load(int_modlist_full[i]) # load model
      intmod <- out$model
      sp <- out$species

      # correct species names
      if (sp == "Accipiter gentilis"){
        sp <- "Astur atricapillus"
        }

      if (sp == "Spilogale putorius"){
        sp <- "Spilogale interrupta"
        }

      if (sp == "Cervus elaphus"){
        sp <- "Cervus canadensis"
        }

      out_int <- out
      fe_add <- fixef(out_int$model) #get fixed effects
      mod <- "interactive"
      variable <- "area size"

      caption <- glue("Figure S{fig_nums[i]}. Model diagnostics plots for {sp} {mod} {variable} models. A) Posterior predictive distribution (y rep)\n
                       compared to {variable} data distribution. B) Posterior predictive errors plot, showing error (y-y_rep) on x axis relative to \n
                       {variable} data distribution (y) on y axis. C: Trace plots showing  individual MCMC draws to evaluate convergence among\n
                       chains for each parameter in the model ('b_intercept' = global intercept, 'b_sg_norm' = effect of human mobility,\n
                       'b_ghm_scale' = effect of human modification, 'b_ndvi_scale' = effect of ndvi, 'b_tmax_scale' = effect of maximum daily \n
                       temperature, 'sd_grp_Intercept' = standard deviation of random intercept by individual-year, 'nu' = Student's t distribution\n
                       degrees of freedom parameter, 'ar[1]' = first-order autoregressive coefficient, 'b_sg_norm:ghm_scale' = interactive effect of\n
                       human modification and mobility [only included when interaction model selected]). D) Table reporting median and 95% credible\n
                       intervals for the posterior distributions, probability of direction ('pd'), Gelman-Rubin diagnostic ('Rhat'), and effective sample\n
                       size ('ESS'), for each parameter (names follow those described for panel C).")

        # convert into a ggplot object so it's compatible
      caption_plot <- ggplot() +
                        geom_text(aes(-0.9, 0, label = caption), size = 3, hjust = 0) +
                        theme_void() + xlim(-1, 1)
      
      #- Extract parameter table as grob -#
      parameters_df <- parameters(addmod) %>% select(-c(CI, Component))
      param_table <- tableGrob(parameters_df)

      #- Posterior Predictive Plot -#
      (pp_dens <- pp_check(addmod)+ggtitle("Posterior Predictive Distribution"))
      
      #- Predictive Error Plot -#
      (pp_err <- pp_check(addmod, type='error_scatter_avg')+ggtitle("Posterior Predictive Errors"))  
    
      #- MCMC Trace Plot -#
      (trace <- mcmc_plot(addmod, type = "trace") +ggtitle("MCMC Traces"))
      
      #---- Assemble Plots ---#
      
      # Design layout
      design <- "AAAAAA
                 BBBCCC
                 DDDDDD
                 DDDDDD
                 EEEEEE
                 ######
                 ######"

      # Gather plots
      (model_out <- wrap_elements(wrap_elements(full = param_table) + ggtitle("Model Coefficient Table") +
                                  pp_dens+pp_err+
                                  trace+
                                  caption_plot+
                                  plot_layout(design = design)))

      # Write out plot
      ggsave(model_out, filename = glue(.wd, "/out/model_diagnostics_final/area/{sp}.pdf"), 
              width = 10, height = 20, device = cairo_pdf)
      
    } # else collect the interactions
  } else {#if int is NULL...

    message("int_modlist_full is NULL")

  } #else
  
}# i 

message("all done....")

