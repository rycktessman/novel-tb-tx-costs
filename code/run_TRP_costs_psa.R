library(dplyr)
library(tidyr)

load("data/params_psa.Rda")
load("data/cost_params_psa.Rda")
load("data/regimen_chars_psa.Rda")
load("data/cost_chars_psa.Rda")
source("code/TRP22_functions.R")

#specify country and RS or RR TB
country <- "philippines" #india, southafrica, philippines
tb_type <- "DS"

country <- Sys.getenv('country')
tb_type <- Sys.getenv('tb_type')

print(country)
print(tb_type)

##################################################
# 1. Model Parameters                            #
##################################################
#pull regimen prioritization characteristics
fullpars <- regimens[[tb_type]]

#pull regimen cost characteristics
fullpars <- c(fullpars, cost_chars[[country]][[tb_type]])

#pull non-regimen-varying treatment model parameters - propogate over all regimens
fullpars <- c(fullpars, params)
n.reg <- length(fullpars$names)
n.sims <- nrow(regimens[["DS"]][["efficacy"]])

#keep non-regimen-varying cost parameters separate (same names)
costs <- cost_params[[country]]

#some costs and params vary by DS vs. DR TB
if(tb_type=="DS") {
  costs$dst <- costs$dst_ds
  costs$support <- costs$support_ds
  costs$hosp <- costs$hosp_ds
  costs$oop <- costs$oop_ds
  costs$indirect <- costs$indirect_ds
  costs$indirect_hosp <- costs$indirect_ds_hosp
  costs$indirect_outpatient <- costs$indirect_ds_outpatient
  costs$p_dst <- costs$p_dst_ds
  costs$p_hosp <- costs$p_hosp_ds
  costs$p_liver <- costs$p_liver_ds
  costs$p_pancreas <- costs$p_pancreas_ds
  costs$p_anemia <- costs$p_anemia_ds
  costs$p_neutropenia <- costs$p_neutropenia_ds
  costs$p_qtcf <- costs$p_qtcf_ds
  costs$p_renal <- costs$p_renal_ds
  costs$p_vision <- costs$p_vision_ds
  costs$p_arthralgia <- costs$p_arthralgia_ds
  costs$p_neuro_short <- costs$p_neuro_short_ds
  costs$p_neuro_long <- costs$p_neuro_long_ds
  costs$drugs1 <- matrix(c(costs$drugs_hrze, rep(0, length(costs$drugs_hrze))), ncol=2) #SOC only - line 1 = HRZE
  costs$drugs2 <- matrix(c(costs$drugs_lrze, rep(0, length(costs$drugs_lrze))), ncol=2) #SOC only - line 2 = LevoRZE
  costs$p_resist <- costs$inh_resist
} else if (tb_type=="DR") {
  costs$dst <- costs$dst_dr
  costs$support <- costs$support_dr
  costs$hosp <- costs$hosp_dr
  costs$oop <- costs$oop_dr
  costs$indirect <- costs$indirect_dr
  costs$indirect_hosp <- costs$indirect_dr_hosp
  costs$indirect_outpatient <- costs$indirect_dr_outpatient
  costs$p_dst <- costs$p_dst_dr
  costs$p_hosp <- costs$p_hosp_dr
  costs$p_liver <- costs$p_liver_dr
  costs$p_pancreas <- costs$p_pancreas_dr
  costs$p_anemia <- costs$p_anemia_dr
  costs$p_neutropenia <- costs$p_neutropenia_dr
  costs$p_qtcf <- costs$p_qtcf_dr
  costs$p_renal <- costs$p_renal_dr
  costs$p_vision <- costs$p_vision_dr
  costs$p_arthralgia <- costs$p_arthralgia_dr
  costs$p_neuro_short <- costs$p_neuro_short_dr
  costs$p_neuro_long <- costs$p_neuro_long_dr
  costs$drugs1 <- matrix(c(costs$drugs_bpalm, rep(0, length(costs$drugs_bpalm))), ncol=2) #SOC only - line 1 = BPaLM
  costs$drugs2 <- matrix(c(costs$drugs_bpal, rep(0, length(costs$drugs_bpal))), ncol=2) #SOC only - line 2 = BPaL
  costs$p_resist <- costs$moxi_resist
}

#load serial interval distribution
serial_interval <- read.csv("data/serial_interval.csv") %>% select(year, p_tb_inc) %>%
  rename("prob"="p_tb_inc")


###################################################
# 2. Run Models                                   #
################################################### 

dur_soc_adj <- unique(regimens[[tb_type]][["duration"]][[1]]/6) #to convert from monthly to full course - essentially how many wks in a month
out_all <- list() #loop over novel regimens
for(i in 1:n.reg) {
  name <- fullpars$names[[i]]
  print(name)
  out_all_scen <- list()
  for(j in 1:n.sims) {
    #cat('\r', paste(round(j/n.sims * 100), "% done", sep = " ")) # display the progress of the simulation
    out <- cost_model_psa(fullpars, costs, tb_type, regimen_num=i, index=j) #main model
    out <- calc_secondary_psa(out, fullpars, costs, serial_interval, index=j) #retreatments, non-cures, secondary cases
    
    #additional drug costs that would be incurred - based on % ppl taking them over time (not relevant for SOC)
    #first account for monthly dispensing schedule
    on_tx_drugs <- rep(0, length(out$on_tx_time))
    on_tx_drugs[(seq(1, fullpars$duration[j, i], by=4))] <- out$on_tx_time[(seq(1, fullpars$duration[j, i], by=4))]
    if(fullpars$duration[j,i]%%4!=0) {
      on_tx_drugs[[max(seq(1, fullpars$duration[j,i], by=4))]] <- 0 #if ever > 4 weeks per month (last dose)
    }
    treat_drugs <- on_tx_drugs*(1+costs$wastage[j])
    retreat_drugs <- on_tx_drugs*out$outcomes[["non_cures_alive"]]*(1+costs$wastage[j])
    secondary_drugs <- on_tx_drugs*out$outcomes[["secondary_cases"]]*costs$cdr[j]*(1+costs$wastage[j])
    
    out[["drugs_time"]] <- list("treat_drugs"=treat_drugs,
                                "retreat_drugs"=retreat_drugs,
                                "secondary_drugs"=secondary_drugs)
    
    out_all_scen[[j]] <- out
  }
  out_all[[name]] <- out_all_scen
}

#combine into list of dataframes for each strategy
out <- out_all
out_all <- list()
for(i in names(out)) {
  print(i)
  out_all[[i]] <- list()
  for(j in names(out[[1]][[1]])) {
    print(j) 
    out_tmp <- sapply(1:n.sims, function(x) unlist(out[[i]][[x]][[j]]), simplify=T, USE.NAMES=T)
    if(length(dim(out_tmp))==0) {
      out_tmp <- as.vector(out_tmp)
    } else {
      out_tmp <- as.data.frame(t(out_tmp))
    }
    if(j=="drugs_time") {
      out_tmp <- pivot_longer(out_tmp %>% mutate(id=1:n.sims), 
                              cols=-id, 
                              names_to=c(".value", "time"),
                              names_pattern = "(\\D+)(\\d+)")
      out_tmp <- out_tmp %>% group_by(id) %>% select(-time) %>%
        summarise_all(sum) %>% select(-id)
    }
    out_all[[i]][[j]] <- out_tmp
  }
}


###################################################
# 3. Calculate price thresholds                   #
################################################### 

#cost categories used for each calculation
cost_cats1 <- c("outpatient", "inpatient", "labs", "aes", "support", "drugs", "oop", "indirect")
cost_cats1_hs <- c("outpatient", "inpatient", "labs", "aes", "support", "drugs")

#calculate costs
for(i in names(out_all)) {
  print(i)
  #1. Short-term affordability: per-person cost-neutral price threshold
  out_all[[i]][["cost_patient"]] <- rowSums(out_all[[i]][["costs"]][,cost_cats1])
  out_all[[i]][["cost_patient_hs"]] <- rowSums(out_all[[i]][["costs"]][,cost_cats1_hs])

  #2. Medium-term affordability: 5 year cost-neutral price threshold
  out_all[[i]][["cost_5yr"]] <- rowSums(out_all[[i]][["costs"]])
  out_all[[i]][["cost_5yr_hs"]] <- rowSums(out_all[[i]][["costs_hs"]])
}

#calculate thresholds - skip SOC
for(i in 2:n.reg) {
  #1. Short-term affordability: per-person cost-neutral price threshold
  out_all[[i]][["cost_neutral_pp"]] <- ((out_all[["soc"]][["cost_patient"]] - out_all[[i]][["cost_patient"]])/
    out_all[[i]][["drugs_time"]][["treat_drugs"]])*fullpars$duration[[i]]/dur_soc_adj #convert from monthly to full course
  out_all[[i]][["cost_neutral_pp_hs"]] <- ((out_all[["soc"]][["cost_patient_hs"]] - out_all[[i]][["cost_patient_hs"]])/
    out_all[[i]][["drugs_time"]][["treat_drugs"]])*fullpars$duration[[i]]/dur_soc_adj
  
  #2. Medium-term affordability: 5 year cost-neutral price threshold
  out_all[[i]][["cost_neutral_5yr"]] <- ((out_all[["soc"]][["cost_5yr"]] - out_all[[i]][["cost_5yr"]])/
    rowSums(out_all[[i]][["drugs_time"]]))*fullpars$duration[[i]]/dur_soc_adj
  out_all[[i]][["cost_neutral_5yr_hs"]] <- ((out_all[["soc"]][["cost_5yr_hs"]] - out_all[[i]][["cost_5yr_hs"]])/
    rowSums(out_all[[i]][["drugs_time"]]))*fullpars$duration[[i]]/dur_soc_adj
  
  #3. Medium-term cost-effectiveness: lifetime cost-effective price threshold
  out_all[[i]][["inc_costs"]] <- rowSums(out_all[[i]][["disc_costs"]]) - rowSums(out_all[["soc"]][["disc_costs"]])
  out_all[[i]][["inc_dalys"]] <- rowSums(out_all[["soc"]][["disc_dalys"]]) - rowSums(out_all[[i]][["disc_dalys"]])
  out_all[[i]][["cost_cea"]] <- ((costs$wtp*out_all[[i]][["inc_dalys"]] - out_all[[i]][["inc_costs"]])/
    rowSums(out_all[[i]][["drugs_time"]]))*fullpars$duration[[i]]/dur_soc_adj
  
  out_all[[i]][["inc_costs_hs"]] <- rowSums(out_all[[i]][["disc_costs_hs"]]) - rowSums(out_all[["soc"]][["disc_costs_hs"]])
  out_all[[i]][["cost_cea_hs"]] <- ((costs$wtp*out_all[[i]][["inc_dalys"]] - out_all[[i]][["inc_costs_hs"]])/
                                  rowSums(out_all[[i]][["drugs_time"]]))*fullpars$duration[[i]]/dur_soc_adj
}

#for SOC, "threshold" (to show in sensitivity analysis) would just be current price - just need to calc mix of SOC given resistance
out_all[["soc"]][["cost_neutral_pp"]] <- sum(costs$drugs1)*(1-costs$p_dst) +
  sum(costs$drugs1)*costs$p_dst*(1-fullpars[["resistance"]][[1]]*costs$p_resist) + 
  sum(costs$drugs2)*costs$p_dst*fullpars[["resistance"]][[1]]*costs$p_resist
out_all[["soc"]][["cost_neutral_pp_hs"]] <- out_all[["soc"]][["cost_neutral_pp"]]
out_all[["soc"]][["cost_neutral_5yr"]] <- out_all[["soc"]][["cost_neutral_pp"]]
out_all[["soc"]][["cost_neutral_5yr_hs"]] <- out_all[["soc"]][["cost_neutral_pp"]]
out_all[["soc"]][["cost_cea"]] <- out_all[["soc"]][["cost_neutral_pp"]]
out_all[["soc"]][["cost_cea_hs"]] <- out_all[["soc"]][["cost_neutral_pp"]]


#combine cost threshold estimates into dataframe
thresholds <- data.frame("name"=unlist(lapply(1:n.reg, function(x) rep(names(out_all)[[x]], n.sims))),
                         "id"=rep(1:n.sims, n.reg),
                         "short_neutral"=unlist(lapply(1:n.reg, function(x) out_all[[x]][["cost_neutral_pp"]])),
                         "medium_neutral"=unlist(lapply(1:n.reg, function(x) out_all[[x]][["cost_neutral_5yr"]])),
                         "cea"=unlist(lapply(1:n.reg, function(x) out_all[[x]][["cost_cea"]])))
thresholds_hs <- data.frame("name"=unlist(lapply(1:n.reg, function(x) rep(names(out_all)[[x]], n.sims))),
                            "id"=rep(1:n.sims, n.reg),
                            "short_neutral"=unlist(lapply(1:n.reg, function(x) out_all[[x]][["cost_neutral_pp_hs"]])),
                            "medium_neutral"=unlist(lapply(1:n.reg, function(x) out_all[[x]][["cost_neutral_5yr_hs"]])),
                            "cea"=unlist(lapply(1:n.reg, function(x) out_all[[x]][["cost_cea_hs"]])))
#save to file
write.csv(thresholds, paste0("output_psa/cost_thresholds_", tb_type, "_", country, ".csv"), row.names=F)
write.csv(thresholds_hs, paste0("output_psa/cost_thresholds_hs_", tb_type, "_", country, ".csv"), row.names=F)
save(out_all, file=paste0("output_psa/out_costs_all_", tb_type, "_", country, ".Rda"))

