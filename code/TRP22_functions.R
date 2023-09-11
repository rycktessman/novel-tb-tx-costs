#version of treatment_model for PSA
treatment_model_psa <- function(fullpars, tb_type="DS", regimen_num=1, index=1) {
  pars <- fullpars[[tb_type]] #filter by DS vs. DR-TB parameters
  with(pars, {
    
    # model starts from when tb is diagnosed 
    # regimen is chosen according to "eligibility" criterion
    startmat <- c(1-eligibility[index, regimen_num], eligibility[index, regimen_num])*(regimen_num!=1) +
      c(1, 0)*(regimen_num==1) # probability of starting SOC (1) or novel (2) regimen
    names(startmat) <- c("SOC", "novel")
    
    #initial loss (0 doses of tx started, doesn't vary by regimen)
    on_tx <- startmat*(1-initialloss[index])
    
    #multiply by regimen "efficacy" (= % cured)
    cured <- on_tx*unlist(efficacy[index, c(1, regimen_num)])
    
    #model LTFU (doesn't vary by regimen) over regimen "duration"
    durs <- unlist(duration[index, c(1, regimen_num)])
    max_dur <- max(durs)
    cured <- t(outer(cured, ((1-ltfu[index])^(0:max_dur))))
    #zero out rows that exceed the duration of the regimen (+1 for week of initiation)
    cured <- cured*cbind(c(rep(1, duration[index, 1]+1), 
                           rep(0, max_dur-duration[index, 1])),
                         c(rep(1, duration[index, regimen_num]+1), 
                           rep(0, max_dur-duration[index, regimen_num])))
    #subtract so that cured represents % that only took regimen for that long
    cured <- cured-rbind(cured[2:(max_dur+1),], c(0,0))
    #multiply by relative efficacy to calc % cured by time of LTFU
    #may need to switch to hazard ratio approach instead
    cured <- model_ltfu(cured, durs, max_dur, as.vector(ltfu_rel_eff[index,]), ltfu_durs)
    cured <- colSums(cured) #no need to track dropout over time anymore
    
    #model "adherence" (varies by regimen)
    #adherence-efficacy relationship determined by "forgiveness" (varies by regimen)
    #first calculate % by # doses taken per week (while still on tx)
    adherence_soc <- as.numeric(adherence[index,,1])
    adherence_novel <- as.numeric(adherence[index,,regimen_num])
    cured <- cbind(adherence_soc, adherence_novel) %*% diag(cured)
    #determine which adherence categories are poor vs. adequate - based on forgiveness
    adequate_adhere_soc <- which(adhere_grps==1-forgive[index, 1])
    adequate_adhere_novel <- which(adhere_grps==1-forgive[index, regimen_num])
    #next multiply by relative efficacy
    cured <- cured*cbind(c(0, rep(rel_eff_forgive[index], adequate_adhere_soc-1), 
                           rep(1, length(adhere_grps)-adequate_adhere_soc)),
                         c(0, rep(rel_eff_forgive[index], adequate_adhere_novel-1),
                           rep(1, length(adhere_grps)-adequate_adhere_novel)))
    #sum over adherence categories and regimens
    cured <- colSums(cured)
    cured <- sum(cured)
    
    return(cured)
  })
}

#model impact of LTFU over regimen duration - used within treatment_model function
model_ltfu <- function(cured, durs, max_dur, ltfu_rel_eff, ltfu_durs) { #version of relapsefracs in original code
    # use relapse %s at 2, 4, and 6 months (defined in pars) for a 6 month regimen, and 100% relapse at 0, 
    # interpolate linearly 0-2, 2-4, and 4-6 to get relapse % after any fraction of treatment course.
    
    # first calculate % of duration completed corresponding to each timestep
    fracs <- cbind(c((0:(durs[1]-1))+1/2, durs[1], 
                     rep(durs[1], max_dur-durs[1])),
                   c((0:(durs[2]-1))+1/2, durs[2], 
                     rep(durs[2], max_dur-durs[2])))%*%diag(1/durs)
    
    # next calculate % efficacy achieved (and resulting % cured) based on % duration completed
    if(is.matrix(ltfu_rel_eff)) {
      loops <- nrow(ltfu_rel_eff)-1
    } else {
      loops <- length(ltfu_rel_eff)-1
    }
    for(i in 1:loops) {
      cured[fracs >= ltfu_durs[i] & fracs < ltfu_durs[i+1]] <-
        ((cured*ltfu_rel_eff[i])*(ltfu_durs[i+1] - fracs))[fracs >= ltfu_durs[i] & fracs < ltfu_durs[i+1]]/(1/(length(ltfu_durs)-1)) +
        ((cured*ltfu_rel_eff[i+1])*(fracs-ltfu_durs[i]))[fracs >= ltfu_durs[i] & fracs < ltfu_durs[i+1]]/(1/(length(ltfu_durs)-1)) 
    }
    return(cured)
}

#version of cost_model for PSA
cost_model_psa <- function(fullpars, costs, tb_type="DS", regimen_num=2, index=1) { #here, 1=SOC, 2=novel
  pars <- fullpars
  with(pars, {
    
    # model starts from when tb is diagnosed 
    # regimen is chosen according to "eligibility" criterion
    startmat <- c(1-eligibility[index, regimen_num], eligibility[index, regimen_num])*(regimen_num!=1) +
      c(1, 0)*(regimen_num==1) # probability of starting SOC (1) or novel (2) regimen
    names(startmat) <- c("SOC", "novel")
    
    #calculate TB dalys by multiplying over regimen duration - incurred regardless of tx continuation/cure
    dalys_tb <- sum(startmat*duration[index, c(1, regimen_num)]*costs$dalys_tb[index])/52
    dalys_posttb <- costs$dalys_posttb[index] #post TB dalys are incurred no matter what
    
    #initial loss (0 doses of tx started, doesn't vary by regimen)
    on_tx <- startmat*(1-initialloss[index])

    #model LTFU (doesn't vary by regimen) over regimen "duration"
    durs <- unlist(duration[index, c(1, regimen_num)])
    max_dur <- max(durs)
    on_tx <- t(outer(on_tx, ((1-ltfu[index])^(0:max_dur))))
    #zero out rows that exceed the duration of the regimen (+1 for week of initiation)
    on_tx <- on_tx*cbind(c(rep(1, duration[index, 1]+1), 
                           rep(0, max_dur-duration[index, 1])),
                         c(rep(1, duration[index, regimen_num]+1), 
                           rep(0, max_dur-duration[index, regimen_num])))
    
    #apply costs here - % of ppl still on tx over time 
    outpatient_costs <- on_tx*outpatient[,c(1, regimen_num)]*costs$outpatient[index]
    inpatient_costs <- on_tx*hosp[,c(1, regimen_num)]*costs$hosp[index]*costs$p_hosp[index] 
    lab_costs <- on_tx*(dst[,c(1, regimen_num)]*costs$p_dst[index]*costs$dst[index] + 
                          smear[,c(1, regimen_num)]*costs$smear[index] +
                          culture[,c(1, regimen_num)]*costs$culture[index] + 
                          xpert[,c(1, regimen_num)]*costs$xpert[index] +
                          cxr[,c(1, regimen_num)]*costs$cxr[index] + 
                          liver_test[,c(1, regimen_num)]*costs$liver_test[index] +
                          full_blood[,c(1, regimen_num)]*costs$full_blood[index] + 
                          ecg[,c(1, regimen_num)]*costs$ecg[index] +
                          neuro_screen[,c(1, regimen_num)]*costs$neuro_screen[index])
    ae_costs <- on_tx*(ae_inc[,c(1, regimen_num)]*costs$p_liver[index]*costs$liver_ae[index] + 
                         ae_inc[,c(1, regimen_num)]*costs$p_pancreas[index]*costs$pancreatitis[index] +
                         ae_inc[,c(1, regimen_num)]*costs$p_anemia[index]*costs$anemia[index] + 
                         ae_inc[,c(1, regimen_num)]*costs$p_neutropenia[index]*costs$neutropenia[index] +
                         ae_inc[,c(1, regimen_num)]*costs$p_qtcf[index]*costs$qtcf[index] + 
                         ae_inc[,c(1, regimen_num)]*costs$p_renal[index]*costs$renal_ae[index] +
                         ae_inc[,c(1, regimen_num)]*costs$p_vision[index]*costs$vision_ae[index] + 
                         ae_inc[,c(1, regimen_num)]*costs$p_arthralgia[index]*costs$arthralgia[index] +
                         ae_inc[,c(1, regimen_num)]*costs$p_neuro_short[index]*costs$neuro_ae[index])
    support_costs <- on_tx*support[,c(1, regimen_num)]*costs$support[index]
    drug_costs <- (on_tx*costs$drugs1*(1-costs$p_resist[index]*resistance[c(1, regimen_num)]) + #non-resistant get first line
      on_tx*costs$drugs1*costs$p_resist[index]*(resistance[c(1, regimen_num)]*(1-costs$p_dst[index]*dst[1,c(1, regimen_num)])) + #undetected resistant get first line
      on_tx*costs$drugs2*costs$p_resist[index]*(resistance[c(1, regimen_num)]*costs$p_dst[index]*dst[1,c(1, regimen_num)]))* #detected resistant get second line
      (1+costs$wastage[index]) #add wastage to all drug components
    oop_costs <- on_tx*costs$oop[index,] 
    indirect_costs <- on_tx*costs$indirect[index,] + 
      on_tx*hosp[,c(1, regimen_num)]*costs$indirect_hosp[index]*costs$p_hosp[index] + 
      on_tx*outpatient[,c(1, regimen_num)]*costs$indirect_outpatient[index]
    
    outpatient_costs <- sum(outpatient_costs)
    inpatient_costs <- sum(inpatient_costs)
    lab_costs <- sum(lab_costs)
    ae_costs <- sum(ae_costs)
    support_costs <- sum(support_costs)
    drug_costs <- sum(drug_costs)
    oop_costs <- sum(oop_costs)
    indirect_costs <- sum(indirect_costs)
    costs_out <- list("outpatient"=outpatient_costs,
                      "inpatient"=inpatient_costs,
                      "labs"=lab_costs,
                      "aes"=ae_costs,
                      "support"=support_costs,
                      "drugs"=drug_costs,
                      "oop"=oop_costs,
                      "indirect"=indirect_costs
    )
    
    on_tx_novel <- on_tx[,2] #save this for calculating the drug price thresholds - weekly price is multiplied by this
    
    #calculate AE dalys here - all AEs last for 1 month except short-term neuropathy (3 months)
    dalys_ae <- on_tx*(ae_inc[,c(1, regimen_num)]*costs$p_liver[index]*costs$dalys_liver[index]/12 + 
                         ae_inc[,c(1, regimen_num)]*costs$p_pancreas[index]*costs$dalys_pancreatitis[index]/12 +
                         ae_inc[,c(1, regimen_num)]*costs$p_anemia[index]*costs$dalys_anemia[index]/12 + 
                         ae_inc[,c(1, regimen_num)]*costs$p_neutropenia[index]*costs$dalys_neutropenia[index]/12 +
                         ae_inc[,c(1, regimen_num)]*costs$p_qtcf[index]*costs$dalys_qtcf[index]/12 + 
                         ae_inc[,c(1, regimen_num)]*costs$p_renal[index]*costs$dalys_renal[index]/12 +
                         ae_inc[,c(1, regimen_num)]*costs$p_vision[index]*costs$dalys_vision[index]/12 + 
                         ae_inc[,c(1, regimen_num)]*costs$p_arthralgia[index]*costs$dalys_arthralgia[index]/12 +
                         ae_inc[,c(1, regimen_num)]*costs$p_neuro_short[index]*costs$dalys_neuro_short[index]/4)
    #assume only those on tx for at least a month incur long-term disability from it
    dalys_ae <- sum(dalys_ae) + sum(on_tx[5,]*neuro_long[c(1, regimen_num)]*
                                      costs$p_neuro_long[index]*costs$disc_dalys_neuro_long[index]) #only long-term AE; already discounted (short-term AEs don't get discounted)
    dalys_out <- list("tb"=dalys_tb,
                      "posttb"=dalys_posttb,
                      "aes"=dalys_ae)

   
    
    #subtract so that on_tx represents % that only took regimen for that long
    on_tx <- on_tx-rbind(on_tx[2:(max_dur+1),], c(0,0))
    
    #multiply by relative efficacy to calc % cured by time of LTFU
    #may need to switch to hazard ratio approach instead
    cured <- model_ltfu(on_tx, durs, max_dur, as.vector(ltfu_rel_eff[index,]), ltfu_durs)
    #multiply by regimen "efficacy" (= % cured)
    cured <- t(t(cured)*unlist(efficacy[index, c(1, regimen_num)]))
    cured <- colSums(cured) #no need to track dropout over time anymore
    
    #model "adherence" (varies by regimen)
    #adherence-efficacy relationship determined by "forgiveness" (varies by regimen)
    #first calculate % by # doses taken per week (while still on tx)
    adherence_soc <- as.numeric(adherence[index,,1])
    adherence_novel <- as.numeric(adherence[index,,regimen_num])
    cured <- cbind(adherence_soc, adherence_novel) %*% diag(cured)
    #determine which adherence categories are poor vs. adequate - based on forgiveness
    adequate_adhere_soc <- which(adhere_grps==1-forgive[index, 1])
    adequate_adhere_novel <- which(adhere_grps==1-forgive[index, regimen_num])
    #next multiply by relative efficacy
    cured <- cured*cbind(c(0, rep(rel_eff_forgive[index], adequate_adhere_soc-1), 
                           rep(1, length(adhere_grps)-adequate_adhere_soc)),
                         c(0, rep(rel_eff_forgive[index], adequate_adhere_novel-1),
                           rep(1, length(adhere_grps)-adequate_adhere_novel)))
    #sum over adherence categories and regimens
    cured <- colSums(cured)
    cured <- sum(cured)
    
    out <- list("cured"=cured,
                "costs"=costs_out,
                "dalys"=dalys_out,
                "on_tx_time"=on_tx_novel)
    
    return(out)
  })
  
}

#version of calc_secondary for PSA
calc_secondary_psa <- function(out, fullpars, costs, serial_interval, index=1) {
  
  out$outcomes <- list() #to save outcomes other than costs and dalys
  out$costs_hs <- out$costs[c("outpatient", "inpatient", "labs", "aes", "support", "drugs")] #health-system only costs
  
  #secondary cases and retreatments
  non_cures <- 1 - out$cured
  non_cures_die <- non_cures*costs$cfr_tx[index] + #% of tx failure resulting in [avoidable] mortality
    non_cures*(1-costs$cfr_tx[index])*costs$cfr[index] #CFR applied to those that don't die during tx
  non_cures_alive <- non_cures - non_cures_die
  secondary_cases <- non_cures_alive*costs$R_eff[index]
  out$outcomes[["non_cures_alive"]] <- non_cures_alive
  out$outcomes[["secondary_cases"]] <- secondary_cases
  
  #secondary and retreatment costs - what gets saved is in 1st 5 years only
  treat_costs <- sum(unlist(out$costs))
  retreat_costs <- non_cures_alive*treat_costs
  secondary_costs <- secondary_cases*(1-costs$cfr[index])*treat_costs #everyone who doesn't die gets treated eventually
  secondary_costs_retreat <- secondary_cases*(1-costs$cfr[index])*retreat_costs
  out$costs[["retreat_costs"]] <- retreat_costs
  out$costs[["secondary_costs"]] <- secondary_costs*sum(serial_interval %>% filter(year<5) %>% pull(prob)) + #for 5 yr horizon, only secondary cases in 1st 5 years
    secondary_costs_retreat*sum(serial_interval %>% mutate(year=year+costs$delay_retreat[index]) %>%
                                  filter(year<5) %>% pull(prob))
  
  treat_costs_hs <- sum(unlist(out$costs_hs))
  retreat_costs_hs <- non_cures_alive*treat_costs_hs
  secondary_costs_hs <- secondary_cases*(1-costs$cfr[index])*treat_costs_hs #everyone who doesn't die gets treated eventually
  secondary_costs_retreat_hs <- secondary_cases*(1-costs$cfr[index])*retreat_costs_hs
  out$costs_hs[["retreat_costs"]] <- retreat_costs_hs
  out$costs_hs[["secondary_costs"]] <- secondary_costs_hs*sum(serial_interval %>% filter(year<5) %>% pull(prob)) + #for 5 yr horizon, only secondary cases in 1st 5 years
    secondary_costs_retreat_hs*sum(serial_interval %>% mutate(year=year+costs$delay_retreat[index]) %>%
                                     filter(year<5) %>% pull(prob))
  
  #discounted costs - could be higher because over longer horizon
  out$disc_costs <- out$costs #treat_costs don't get discounted
  out$disc_costs[["retreat_costs"]] <- retreat_costs*(1/(1+costs$disc_rate[index]))^(costs$delay_retreat[index])
  #for secondary cases, apply serial interval distribution
  lag_secondary_costs <- secondary_costs*(1/(1+costs$disc_rate[index]))^(serial_interval$year + 1/(-log(1-costs$cdr[index]))) #add treatment delay - 1 over case detection rate
  lag_secondary_retreat_costs <- secondary_costs_retreat*(1/(1+costs$disc_rate[index]))^
    (serial_interval$year + costs$delay_retreat[index] + 1/(-log(1-costs$cdr[index]))) #further add the retx delay
  out$disc_costs[["secondary_costs"]] <- sum(lag_secondary_costs*serial_interval$prob) + 
    sum(lag_secondary_retreat_costs*serial_interval$prob)
  
  out$disc_costs_hs <- out$costs_hs #treat_costs don't get discounted
  out$disc_costs_hs[["retreat_costs"]] <- retreat_costs_hs*(1/(1+costs$disc_rate[index]))^costs$delay_retreat[index]
  #for secondary cases, apply serial interval distribution
  lag_secondary_costs_hs <- secondary_costs_hs*(1/(1+costs$disc_rate[index]))^(serial_interval$year + 1/(-log(1-costs$cdr[index])))
  lag_secondary_retreat_costs_hs <- secondary_costs_retreat_hs*(1/(1+costs$disc_rate[index]))^
    (serial_interval$year + costs$delay_retreat[index] + 1/(-log(1-costs$cdr[index]))) #further add the retx delay
  out$disc_costs_hs[["secondary_costs"]] <- sum(lag_secondary_costs_hs*serial_interval$prob) + 
    sum(lag_secondary_retreat_costs_hs*serial_interval$prob)
  
  #deaths and disability
  treat_dalys <- sum(unlist(out$dalys)) #apply for secondary cases - include postTB dalys
  treat_dalys_nopost <- sum(unlist(out[["dalys"]][c("tb", "aes")])) #don't 2x count postTB dalys for retx's
  non_cures_dalys <- non_cures_die*costs$disc_life_exp[index] + 
    non_cures_alive*treat_dalys_nopost + #those that treated incur tx DALYs, minus postTB to avoid 2x counting
    non_cures*(1-costs$cfr_tx[index])*costs$dalys_tb[index]*costs$time_symptom[index]/12 #both later deaths and tx initiators incur DALYs
  out$dalys[["non_cures"]] <- non_cures_dalys
  secondary_die <- secondary_cases*costs$cfr[index]
  secondary_dalys <- secondary_die*costs$disc_life_exp[index] + 
    secondary_cases*(1-costs$cfr[index])*treat_dalys +
    secondary_cases*(1-costs$cfr[index])*non_cures_dalys +
    secondary_cases*costs$dalys_tb[index]*costs$time_symptom[index]/12 #both deaths and tx initiators incur DALYs
  out$dalys[["secondary"]] <- secondary_dalys
  out$outcomes[["non_cures_die"]] <- non_cures_die
  out$outcomes[["secondary_die"]] <- secondary_die
  
  #version with dalys stratified by reason
  out$dalys2 <- list()
  out$dalys2[["tb"]] <- out$dalys[["tb"]] + non_cures_alive*out$dalys[["tb"]] +
    secondary_cases*(1-costs$cfr[index])*out$dalys[["tb"]] +
    non_cures*(1-costs$cfr_tx[index])*costs$dalys_tb[index]*costs$time_symptom[index]/12 +
    secondary_cases*costs$dalys_tb[index]*costs$time_symptom[index]/12
  out$dalys2[["aes"]] <- out$dalys[["aes"]] + 
    non_cures_alive*out$dalys[["aes"]] +
    secondary_cases*(1-costs$cfr[index])*out$dalys[["aes"]]
  out$dalys2[["posttb"]] <- out$dalys[["posttb"]] + 
    secondary_cases*(1-costs$cfr[index])*out$dalys[["posttb"]]
  out$dalys2[["deaths"]] <- non_cures_die*costs$disc_life_exp[index] + secondary_die*costs$disc_life_exp[index]
  out$dalys2[["secondary_retreats"]] <- secondary_cases*(1-costs$cfr[index])*non_cures_dalys #hard to disaggregate here so just make this a new category
  
  #discounted dalys
  out$disc_dalys <- out$dalys #treat_dalys don't get discounted
  #non-cures & secondary cases
  out$disc_dalys[["non_cures"]] <- non_cures_die*costs$disc_life_exp[index]*(1/(1+costs$disc_rate[index]))^(costs$delay_retreat[index]) + 
    non_cures_alive*treat_dalys_nopost*(1/(1+costs$disc_rate[index]))^(costs$delay_retreat[index]) +
    (non_cures*(1-costs$cfr_tx[index])*costs$dalys_tb[index]*costs$time_symptom[index]/12)*
    (1/(1+costs$disc_rate[index]))^(costs$delay_retreat[index] - costs$time_symptom[index]/12)
  #for secondary cases, apply serial interval distribution
  lag_secondary_dalys <- secondary_die*costs$disc_life_exp[index]*(1/(1+costs$disc_rate[index]))^(serial_interval$year + -log(1-costs$cdr[index])) + 
    secondary_cases*(1-costs$cfr[index])*treat_dalys*(1/(1+costs$disc_rate[index]))^(serial_interval$year + -log(1-costs$cdr[index])) +
    (secondary_cases*costs$dalys_tb[index]*costs$time_symptom[index]/12)*(1/(1+costs$disc_rate[index]))^
    (serial_interval$year + -log(1-costs$cdr[index]) - costs$time_symptom[index]/12)
  lag_secondary_dalys_retreat <- secondary_cases*(1-costs$cfr[index])*non_cures_dalys*(1/(1+costs$disc_rate[index]))^
    (serial_interval$year + costs$delay_retreat[index] + 1/(-log(1-costs$cdr[index]))) #further add the retx delay
  out$disc_dalys[["secondary"]] <- sum(lag_secondary_dalys*serial_interval$prob) +
    sum(lag_secondary_dalys_retreat*serial_interval$prob)
  
  return(out)
}
