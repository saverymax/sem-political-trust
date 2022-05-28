# Data from
# https://www.europeansocialsurvey.org/downloadwizard/

library("dplyr")
#library("psych")
#library('stringr')
library("lavaan")
library("corrplot")
library("tidySEM")
library("kableExtra")

#data_file <- "ESS1-9e01_1/ESS1-9e01_1_all_dates.csv"
##user_dir <- "d:/asus_documents/ku_leuven/courses/structural_equations/project/"
#df <- read.csv(paste(user_dir, data_file, sep=""))
#df <- read.csv(data_file, sep="")
#df <- read.csv(data_file, sep=",")
#View(df)

# Load data into global env
source("sem_utils.R")
source("sem_functions.R")
df_subset <- data_processing(df)
exploratory_plots(df_subset)
fit_trust <- fit_trust_model(df_subset)
caption <- "One-factor CFA for trust in government"
label <- "one-fact-trust"
print_fit(fit_trust, caption, label)
caption <- "Modification indices for one-factor model. Top 10 are shown."
label <- "one-fact-mi"
print_mi(fit_trust, caption, label)

# model adjusted by mi
fit_trust_adj <- fit_trust_model_adj(df_subset)
caption <- paste("One-factor CFA for trust in government, allowing for residual covariances")
label <- paste("one-fact-trust-adj")
print_fit(fit_trust_adj, caption, label)
compare_one_factor(fit_trust, fit_trust_adj)

# Add in poltical variables and trust in people
fit_politics <- fit_politics_model(df_subset)
caption <- "Three-factor model for trust in government, trust in people, and political ability"
label <- "three-fact-poli"
print_fit(fit_politics, caption, label)

# 3 -> 2 factor
fit_coll <- fit_politics_collapsed(df_subset)
caption <- "Collapsed model, where trust becomes 1 factor"
label <- "two-fact-poli"
print_fit(fit_coll, caption, label)

# Then compare them
compare_poli(fit_politics, fit_coll)
caption <- "Modification indices for three-factor model. Top 10 are shown."
label <- "three-fact-mi"
print_mi(fit_politics, caption, label)

# Change fit by mi
fit_adj_1 <- fit_improved_politics_model(df_subset)
caption <- paste("Modified three-factor model for trust in government, trust in people, and political ability. Residual covariances are added and trstplt is removed.")
label <- paste("three-fact-mod-1")
print_fit(fit_adj_1, caption, label)
caption <- paste("After the first set of model changes, we consider again the modification
                     indices for the three-factor model. Top 10 are shown.")
label <- paste("three-fact-mi-2")
print_mi(fit_adj_1, caption, label)

# Change again
fit_adj_2 <- fit_improved_politics_model_2(df_subset)
caption <- paste("Modified three-factor model for trust in government, trust in people, and political ability. Variable actrolga is removed, as well as its associated covariance.")
label <- paste("three-fact-mod-2")
print_fit(fit_adj_2, caption, label)
compare_poli_adj(fit_adj_1, fit_adj_2)

graph_model(fit_adj_2)

# Then run mediation
fit_med <- fit_mediation_model(df_subset)
caption <- "Three-factor medidation model, considering structural relationship
    between trust in people and trust in government"
label <- "three-fact-mediation"
print_fit(fit_med, caption, label)
