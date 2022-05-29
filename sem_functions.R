##############################################################################
# Script containing all functions for running Structural Equations Modelling
# Will be called by either main R script or markdown.
##############################################################################


exploratory_plots <- function(df_subset){
    cov_mat <- cov(df_subset, use = "pairwise.complete.obs")
    corr_mat <- cov2cor(cov_mat)
    p <- corrplot::corrplot(corr_mat,
                       is.corr = T,       # whether is a correlation matrix
                       method = "circle",     # magnitude of covariances as circles
                       type = "upper",        # remove the bottom of the covariance matrix
    )
    p
}

print_fit <- function(fit, caption, label){
    fit_print <- standardizedsolution(fit) %>% mutate_if(is.numeric, round, digits=3)
    if (label=="three-fact-mediation"){
      indices <- c(1:3, 4:8)
    }else{
      indices <- c(1:3, 4:7)
    }
    print(kbl(fit_print[,indices], booktabs = T, escape=T, caption=caption, label=label, linesep = "") %>%
              kable_styling(latex_options = "HOLD_position"))
}

print_mi <- function(fit, caption, label){
    # Get modification indices
    mi<-inspect(fit, "mi")
    mi.sorted <- mi[order(-mi$mi), ] %>% mutate_if(is.numeric, round, digits=3)
    print(kbl(mi.sorted[1:10,1:4], booktabs = T, escape=T, caption=caption, label=label, linesep = "") %>%
              kable_styling(latex_options = "HOLD_position"))
}

fit_trust_model <- function(df_subset){
    trust_model <- '
        # Trust in government
        gov_trust =~ trstprl + trstprt + trstplt  + trstep + trstun
    '

    fit_trust <- cfa(trust_model,
                          data = df_subset)
    #print(fit_summary$FIT)
    #print(fit_summary$PE)
    #fit_print <- parameterEstimates(fit_trust)
    return(fit_trust)
}

# Remove trstplt because it has many mi indices
fit_trust_model_adj <- function(df_subset){
    trust_model_adj <- '
        # Trust in government
        gov_trust =~ trstprl + trstprt + trstep + trstun
        trstep ~~ trstun
    '
    fit_trust_adj <- cfa(trust_model_adj,
                          data = df_subset)
    #fit_summary <- summary(fit_trust_adj, fit.measures=T, standardized=T)
    #print(fit_summary$FIT)
    #print(fit_summary$PE)
    return(fit_trust_adj)
}

compare_one_factor <- function(fit_trust, fit_trust_model_adj){
    m1 <- fitMeasures(fit_trust, c("logl","AIC", "BIC", "chisq", "df", "pvalue", "cfi", "tli","rmsea"), output = "matrix")
    m2 <- fitMeasures(fit_trust_adj, c("logl","AIC", "BIC", "chisq", "df", "pvalue", "cfi", "tli","rmsea"), output = "matrix")
    df_compare <- data.frame("One-factor"=round(m1[,1],3), "One-factor-covariances"=round(m2[,1],3))
    names(df_compare) <- c("One-factor CFA", "One-factor with covariances")
    caption <- paste("Comparison of One-factor CFA models")
    label <- paste("one-fact-compare")
    print(kbl(df_compare, booktabs = T, escape=T, caption=caption, label=label, linesep = "") %>%
              kable_styling(latex_options = "HOLD_position"))
}

fit_politics_model <- function(df_subset){
    # not using polintr
    politics_model <- '
        # Specify latent factors
        # Trust in government
        gov_trust =~ trstprl + trstprt + trstplt  + trstep + trstun
        # Trust in people
        people_trust =~ ppltrst + pplfair + pplhlp
        # political ability
        political_able =~ psppipla + actrolga + psppsgva + cptppola
    '

    fit_politics <- cfa(politics_model,
                          data = df_subset)
    return(fit_politics)
}


fit_politics_collapsed <- function(df_subset){
    # In the first column (labeled Std.lv), only the latent variables are standardized.
    # In the second column (labeled Std.all), both latent and observed variables are standardized.
    # The latter is often called the ‘completely standardized solution’.

    # Not sure if we really need three factors here, so compare nested models
    model_collapsed <- '
        # Specify latent factors
        # Trust
        trust =~ trstprl + trstprt + trstplt + trstep + trstun + ppltrst + pplfair + pplhlp
        # political ability
        political_able =~ psppipla + actrolga + psppsgva + cptppola
        # Add residual covariance
    '
    fit_collapsed <- cfa(model_collapsed,
                          data = df_subset)
    return(fit_collapsed)
}

compare_poli <- function(fit_politics, fit_collapsed){
    anova_out <- anova(fit_politics, fit_collapsed)
    caption <- paste("Comparison of 2- and 3-factor CFA models for trust and political ability")
    label <- paste("collapse-compare")
    print(kbl(anova_out, booktabs = T, escape=T, caption=caption, label=label, linesep = "") %>%
              kable_styling(latex_options = "HOLD_position"))
    # This indicates that the models are indeed different

    m1 <- fitMeasures(fit_politics, c("logl","AIC", "BIC", "chisq", "df", "pvalue", "cfi", "tli","rmsea"), output = "matrix")
    m2 <- fitMeasures(fit_collapsed, c("logl","AIC", "BIC", "chisq", "df", "pvalue", "cfi", "tli","rmsea"), output = "matrix")
    df_compare <- data.frame("Three-factor-CFA"=round(m1[,1],3), "Two-factor-CFA"=round(m2[,1],3))
    names(df_compare) <- c("Three-factor CFA", "Two-factor CFA")
    caption <- paste("Comparison of 2- and 3-factor CFA models")
    label <- paste("three-fact-compare")
    print(kbl(df_compare, booktabs = T, escape=T, caption=caption, label=label, linesep = "") %>%
              kable_styling(latex_options = "HOLD_position"))
}


fit_improved_politics_model <- function(df_subset){
    # Remove trstplt and add in residual covariances for a few variables.
    # Note that adding in psppipla ~~ psppsgva gives a negative variance estimate
    # and when I account for the covariance actrolga ~~ cptppola the high mi for
    # psppipla ~~ psppsgva goes away anyway.
    politics_model_adj_1 <- '
        # Specify latent factors
        # Trust in government
        gov_trust =~ trstprl + trstprt + trstep + trstun
        # Trust in people
        people_trust =~ ppltrst + pplfair + pplhlp
        # political ability
        political_able =~ psppipla + actrolga + psppsgva + cptppola
        # Add residual covariance
        actrolga ~~ cptppola
        trstep ~~ trstun
    '

    fit_politics_adj_1 <- cfa(politics_model_adj_1,
                          data = df_subset)
    # Well this model decreases the chi sequare
    # And also improves other  metrics.
    return(fit_politics_adj_1)
}

fit_improved_politics_model_2 <- function(df_subset){
    politics_model_adj_2 <- '
        # Trust in government
        gov_trust =~ trstprl + trstprt + trstep + trstun
        # Trust in people
        people_trust =~ ppltrst + pplfair + pplhlp + trstep
        # political ability
        political_able =~ psppipla + actrolga + psppsgva + cptppola
        # Add residual covariance
        actrolga ~~ cptppola
        trstep ~~ trstun
        trstprt ~~ trstep
    '

    fit_politics_adj_2 <- cfa(politics_model_adj_2,
                          data = df_subset)
    return(fit_politics_adj_2)
}

compare_poli_adj <- function(fit_politics, fit_politics_adj_1, fit_politics_adj_2){
    m0=fitMeasures(fit_politics, c("logl","AIC", "BIC", "chisq", "df", "pvalue", "cfi", "tli","rmsea", "srmr"), output = "matrix")
    m1=fitMeasures(fit_politics_adj_1, c("logl","AIC", "BIC", "chisq", "df", "pvalue", "cfi", "tli","rmsea", "srmr"), output = "matrix")
    m2=fitMeasures(fit_politics_adj_2, c("logl","AIC", "BIC", "chisq", "df", "pvalue", "cfi", "tli","rmsea", "srmr"), output = "matrix")
    df_compare <- data.frame("fit1"=round(m0[,1],3), "fit2"=round(m1[,1],3), "fit3"=round(m2[,1],3))
    names(df_compare) <- c("Original fit", "First modified fit", "Second modified fit")
    # Well that doesn't make much difference.
    caption <- paste("Comparison of modified 3-factor CFA models")
    label <- paste("mod-three-fact-compare")
    print(kbl(df_compare, booktabs = T, escape=T, caption=caption, label=label, linesep = "") %>%
              kable_styling(latex_options = "HOLD_position"))
}

graph_model <- function(fit){
    lay <- get_layout(
    "ppltrst", "pplfair", "pplhlp", "", "", "",
    "", "people_trust", "trstprl", "trstprt", "trstep", "trstun",
    "", "polintr", "", "gov_trust", "", "",
    "", "", "political_able", "", "", "",
    "", "psppipla",  "actrolga", "psppsgva", "cptppola", "",
    rows = 5)
  p <- graph_sem(model = fit,
          layout = lay,
          angle = 170
          #label = "est_std"
  )
  p
}

fit_mediation_model <- function(df_subset){
    # So then add in mediation model.
    politics_mediation_model <- '
        # Trust in government
        gov_trust =~ trstprl + trstprt + trstep + trstun
        # Trust in people
        people_trust =~ ppltrst + pplfair + pplhlp
        # political ability
        political_able =~ psppipla + actrolga + psppsgva + cptppola
        # Add residual covariance
        actrolga ~~ cptppola
        trstep ~~ trstun

        ## Direct effect(s) ##
        political_able ~ c1*polintr + c2*people_trust

        ## Mediator ##
        ## Path A
        gov_trust ~ a1*people_trust + a2*polintr
        ## Path B
        political_able ~ b1*gov_trust
        ## Indirect effect (a*b) ##
        ab_trust := a1*b1
        ab_interest := a2*b1
        ## Total effect ##
        total_interest := c1 + ab_interest
        total_trust := c2 + ab_trust
    '
    fit_mediation <- cfa(politics_mediation_model,
                          data = df_subset)
    # Well this is sort of an interesting idea, see drawing in notebook.
    # It says that there is an influential effect of people_trust on gov_trust
    # and that there is a significant effect of political interest on political_able
    # There is also a mediation effect from polintr -> gov_trust -> political_able
    return(fit_mediation)
}

graph_mediation <- function(fit){
  lay <- get_layout(
    "ppltrst", "pplfair", "pplhlp", "", "", "",
    "people_trust", "", "trstprl", "trstprt", "trstep", "trstun",
    "", "polintr", "", "gov_trust", "", "",
    "", "", "political_able", "", "", "",
    "", "psppipla",  "actrolga", "psppsgva", "cptppola", "",
    rows = 5)
  p <- graph_sem(model = fit,
          layout = lay,
          angle = 170
          #label = "est_std"
  )
  p
}

print_med_fit_measures <- function(fit_med){
    m1 <- fitMeasures(fit_med, c("logl","AIC", "BIC", "chisq", "df", "pvalue", "cfi", "tli","rmsea"), output = "matrix")
    df_compare <- data.frame("med"=round(m1[,1],3))
    names(df_compare) <- c("Mediation model")
    caption <- paste("Fit measures for mediation model")
    label <- paste("med-fit")
    print(kbl(df_compare, booktabs = T, escape=T, caption=caption, label=label, linesep = "") %>%
              kable_styling(latex_options = "HOLD_position"))
}


# Not including this in analysis
fit_multigroup_model <- function(df_subset){

    politics_model_adj_2 <- '
        # Specify latent factors
        # Trust in government
        gov_trust =~ trstprl + trstprt + trstep + trstun
        # Trust in people
        people_trust =~ ppltrst + pplfair + pplhlp
        # political ability
        political_able =~ psppipla + psppsgva + cptppola
        # Add residual covariance
        actrolga ~~ cptppola
        trstep ~~ trstun
    '
    # fit multigroup model.
    df_subset$gndr <- factor(df_subset$gndr,
                          levels = c("1", "2"),         # levels
                          labels = c("Male", "Female")) # labels

    fit_multi <- cfa(politics_model_adj_2,
                          data = df_subset , group="gndr")
    return(fit_multi)
}
