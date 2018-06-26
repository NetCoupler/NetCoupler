
# Summary statistics and multiple testing adjustment for net.coupler.out with survival object#

mult.stat.surv <- function(sum_netout, adjust_method, rule_1 = 12, rule_1_cut = 0.1, rule_2 = 5, rule_2_cut = 0.05, rule_3 = 15, rule_4 = 14, ass_rule1=16, round_number) {

  # sum_netout: output of getExp.coef.out after adding original metabolite names
  # adjust_method: select adjustment method(s), e.g. "fdr" or c("fdr","bonferroni"), see p.adjust.methods for available options
  # rule= "set DE=1 or =2 if ..."
  # rule_1: determines first argument of rule (e.g. adjusted p-value of full model < rule_1_cut)
  # rule_2: determines second argument of rule (e.g. upperP < rule_2_cut)
  # rule_3: determines first variable of third argument of rule (e.g. highEst)
  # rule_4: determines second variable of third argument of rule (e.g. lowEst)
  # possible arguments for rule_...: adjusted p-values of full model=12=rule_1
  #                                 largest p-value across all adjacency models=upperP=5=rule_2
  #                                 highest estimated beta across all adjacency models=highEst=15=rule_3
  #                                 lowest estimated beta across all adjacency models=lowEst=14=rule_4
  # default rule: if(adjusted p< 0.1 & upperP<0.05 & (sign(highEst)==sign(lowEst)) & lowEst>0){DE=1}
  #              if(adjusted p< 0.1 & upperP<0.05 & (sign(highEst)==sign(lowEst)) & lowEst<0){DE=2}else{DE=0}
  # ass_rule1: determines first variable of second argument of association rule (e.g. bestGuess=16)
  # default association rule: if(adjusted p<0.1 & bestGuess>0){Assoc=1}
  #                          if(adjusted p<0.1 & bestGuess<0){Assoc=2}else{Assoc=0}
  # round_number: actual round number in ambiguous metabolites loop

  # return: for each metabolite: maximum number of models, mean hazard ratio, maximum HR, minimum HR, maximum p-value, minimum p-value, as well as number of covariates, HR and p-value for model with full adjustment:
  sum_netout_sum <- dplyr::left_join(sum_netout %>% dplyr::group_by(Outcome) %>% dplyr::summarise(avgHR = mean(HR), minHR = min(HR), maxHR = max(HR), upperP = max(P), lowerP = min(P), Nbmds = max(Nbmds)),
    sum_netout %>% dplyr::group_by(Outcome) %>% dplyr::filter(nchar(Covariables) == max(nchar(Covariables))) %>% dplyr::select(NCov = NCov, Outcome = Outcome, Metabolite = Metabolite, ind_HR = HR, ind_P = P),
    by = "Outcome"
  )

  # multiple testing correction:
  p_adjust_M <- p.adjust.methods[p.adjust.methods %in% adjust_method] # select multiple testing correction method(s)
  p_adj <- sapply(p_adjust_M, function(meth) {
    p.adjust(sum_netout_sum$ind_P, meth)
  }) # calculate adjusted p-values for ind_P

  sum_netout_sum_adj <- bind_cols(sum_netout_sum, data.frame(p_adj))
  colnames(sum_netout_sum_adj)[12] <- adjust_method

  # add beta-coefficients:
  sum_netout_sum_adj <- sum_netout_sum_adj %>% dplyr::mutate(avgEst = log(avgHR), lowEst = log(minHR), highEst = log(maxHR), bestGuess = log(ind_HR))

  # is there an effect of specified metabolite on time-to-event?: determine effect-indicator (DE!=0: direct effect of metabolite on time-to-event, DE=0: ambiguous if adjusted ind_p<0.1):
  sum_netout_sum_adj <- mutate(sum_netout_sum_adj,
    DE = derivedFactor(
      "1" = (sum_netout_sum_adj[, rule_1] < rule_1_cut & sum_netout_sum_adj[, rule_2] < rule_2_cut & (sign(sum_netout_sum_adj[, rule_3]) == sign(sum_netout_sum_adj[, rule_4])) & sum_netout_sum_adj[, rule_4] > 0),
      "2" = (sum_netout_sum_adj[, rule_1] < rule_1_cut & sum_netout_sum_adj[, rule_2] < rule_2_cut & (sign(sum_netout_sum_adj[, rule_3]) == sign(sum_netout_sum_adj[, rule_4])) & sum_netout_sum_adj[, rule_4] < 0),
      .method = "first", .default = 0
    )
  )

  # is there any effect of specified metabolite on time-to-event?: determine association-indicator (Assoc!=0: there is an effect (no differenciation between direct or indirect), Assoc=0: there is no effect)
  sum_netout_sum_adj <- mutate(sum_netout_sum_adj,
    Assoc = derivedFactor(
      "1" = (sum_netout_sum_adj[, rule_1] < rule_1_cut & sum_netout_sum_adj[, ass_rule1] > 0),
      "2" = (sum_netout_sum_adj[, rule_1] < rule_1_cut & sum_netout_sum_adj[, ass_rule1] < 0),
      .method = "first", .default = 0
    )
  )

  # collect ambiguous effects:
  amb <- sum_netout_sum_adj %>% filter(Assoc != 0 & DE == 0)

  # collect direct effects:
  direct <- sum_netout_sum_adj %>% filter(Assoc != 0 & DE != 0)

  # collect summary statistics including round information
  sum_netout_sum_adj_FV <- sum_netout_sum_adj %>% dplyr::mutate(round = as.numeric(round_number), Assoc_FV = as.character(Assoc), DE_FV = as.character(DE))

  return(list(sum_netout_sum_adj = sum_netout_sum_adj, amb = amb, direct = direct, sum_netout_sum_adj_FV = sum_netout_sum_adj_FV))
}
