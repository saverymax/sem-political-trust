# Set up data processing variables

data_processing <- function(df){
    df_round_select <- df[df$essround==9,]
    # Works with round 9
    select_vars <- c("gndr", "trstprl", "trstlgl", "trstplc", "trstplt", "trstprt", "trstep", "trstun",
                     "psppipla", "actrolga", "psppsgva", "cptppola",
                     "psppsgva", "actrolga", "ppltrst", "pplfair", "pplhlp", "polintr")


    # filter out 77. 88, 99 for refusal, don't know, or no answer.
    # Weird notation in dplyr across
    df %>% filter(essround==9) %>% select(select_vars) %>% filter(across(.cols=everything(), .fns=~ . <= 10)) -> df_subset
    print(df_subset[is.na(df_subset)])
    return(df_subset)
}

data_file <- "ESS1-9e01_1/ESS1-9e01_1_all_dates.csv"
#user_dir <- "d:/asus_documents/ku_leuven/courses/structural_equations/project/"
#df <- read.csv(paste(user_dir, data_file, sep=""))
df <- read.csv(data_file, sep=",")
