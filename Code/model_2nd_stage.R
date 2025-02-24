
library(brms)
library(xtable)

main_path <- "/Users/solon/Cloud-Drive/Projects/ctDNA_original/ctDNA/"

##################################################################   
########################## Load files ############################
################################################################## 
# data
load(paste0(main_path, "data_split/data_train_CT.Rdata")) # train CT dataset

# rename 1st stage model
model_1st_stage <- fit1_TF

##################################################################   
########################## Processing ############################
################################################################## 
# filter Patient.IDs to correspond to stage 1 model
data_train_CT_filtered <- data_train_CT %>% 
    filter(Patient.ID %in% unique(model_1st_stage$data$Patient.ID))

# extract random effects
data_extract_model_1st_stage <- as.data.frame(model_1st_stage$fit@sim$samples[[1]])

# get posterior means for random effects
r_intercept <- data_extract_model_1st_stage %>% select(r_Patient.ID.DT040.Intercept. : r_Patient.ID.DT363.Intercept.) %>% summarise(estim_inter = colMeans(.))
r_slope <- data_extract_model_1st_stage %>% select(r_Patient.ID.DT040.time_ichor.: r_Patient.ID.DT363.time_ichor.) %>% summarise(estim_slope = colMeans(.))

combine_model_1st_stage <- cbind(Patient.ID = unique(model_1st_stage$data$Patient.ID), r_intercept, r_slope)

# combine datasets
data_train_CT_final <- merge(data_train_CT_filtered, combine_model_1st_stage, by = "Patient.ID")

##################################################################   
################# Modelling ###################################### 
##################################################################
prior_custom <- c(set_prior("normal(0, 100)", class = "b"), # fixed effects prior 
                  set_prior("cauchy(0, 2)", class = "sd"))

fit2_CT_TF <- brm(Progression ~ time + ER.status + Her2.status + Treatment_new_final + estim_inter + estim_slope + (1 | Patient.ID),
                         data = data_train_CT_final, 
                         family = bernoulli(link = "logit"), 
                         prior = prior_custom,
                         warmup = 5000,
                         iter = 100000,
                         chains = 2,
                         cores = getOption("mc.cores", 2),
                         control = list(adapt_delta = 0.95, max_treedepth = 12))

fit2_CT_no_TF <- brm(Progression ~ time + ER.status + Her2.status + Treatment_new_final + (1 | Patient.ID),
                        data = data_train_CT_final, 
                        family = bernoulli(link = "logit"), 
                        prior = prior_custom,
                        warmup = 5000,
                        iter = 100000,
                        chains = 2,
                        cores = getOption("mc.cores", 2),
                        control = list(adapt_delta = 0.95, max_treedepth = 12))

############################
# reporting results
#str(summary(fit2_CT_TF))
xtable(summary(fit2_CT_TF)$fixed[,c(1, 3, 4)])
xtable(summary(fit2_CT_TF)$spec_pars[,c(1, 3, 4)])

##################################################################   
################# Post-processing ################################
##################################################################
summary(fit2_CT_TF)

# export to word 
export_summs(fit2_CT_TF, error_format = "[{conf.low}, {conf.high}]", to.file = "Word", file.name = "2nd stage.docx")

