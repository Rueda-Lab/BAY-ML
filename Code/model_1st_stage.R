
library(dplyr)
library(brms)

main_path <- "/Users/solon/Cloud-Drive/Projects/ctDNA_original/ctDNA/"

##################################################################   
########################## Load data #############################
##################################################################   
load(paste0(main_path, "data_split/data_train_ichor.Rdata")) # train IchorCNA dataset

##################################################################   
################# Modelling ###################################### 
##################################################################
prior_custom <- c(set_prior("normal(0, 100)", class = "b"), # fixed effects prior 
                  set_prior("cauchy(0, 2)", class = "sd"), # random effects sd prior
                  set_prior("lkj(1)", class = "cor")) # prior for correlation coef for random effects

fit1_TF <- brm(ichorCNA_tr ~ time_ichor + ER.status + Her2.status + Treatment_new_final + Treatment_duration +
                               (1 + time_ichor | Patient.ID),
                           data = df_train_ichor, 
                           family = gaussian(), 
                           prior = prior_custom,
                           warmup = 2000,
                           iter = 10000,
                           chains = 2,
                           cores = getOption("mc.cores", 2), 
                           control = list(adapt_delta = 0.99, max_treedepth = 12))

fit1_noTF <- brm(ichorCNA_tr ~ time_ichor + ER.status + Her2.status + Treatment_duration +
                      (1 + time_ichor | Patient.ID),
                  data = df_train_ichor, 
                  family = gaussian(), 
                  prior = prior_custom,
                  warmup = 2000,
                  iter = 10000,
                  chains = 2,
                  cores = getOption("mc.cores", 2), 
                  control = list(adapt_delta = 0.99, max_treedepth = 12))


##################################################################   
################# Post-processing ################################
##################################################################

################################ diagnostics, summaries, plots
# summarise output
summary(fit1_TF)

# export to word
export_summs(fit1_TF, error_format = "[{conf.low}, {conf.high}]", to.file = "Word", file.name = "1st_stage.docx")

# plot posterior distributions and chains
#plot(fit1_TF, pars = "^sd")
#plot(fit1_TF, pars = "^sigma")

# plot posterior intervals
#mcmc_plot(fit1_TF)


