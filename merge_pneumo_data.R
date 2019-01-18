#SV Scarpino
#Nov. 2018
#Data downloaded from: https://data.cdc.gov/browse?q=pneumococcal%20&sortBy=relevance
#See citation for pos. association with RSV https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001776#s3

#set working dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))#sets working directory to source file location

#libraries (not included in limits_acc_functions.R)
library(lubridate)
library(MMWRweek)
data(state)

#########
#Globals#
#########
write_new <- FALSE #set to TRUE to save a new csv
time_stamp <- as.numeric(Sys.time())

###########
#acc funcs#
###########
split_years <- function(x){
  x.split <- unlist(strsplit(x = x, split = "_"))
  x.split.1 <- x.split[length(x.split)]
  x.split.2 <- unlist(strsplit(x = x.split.1, split = ".csv"))[1]
  return(x.split.2)
}

############
#Merge Data#
############
#1. List files and get years
files <- list.files("Raw_Data/")

years <- unlist(lapply(files, split_years))

#2. Make data matrix

#get year weeks
years_weeks <- c()
weeks_weeks <- c()
raw_data <- list()
colnames_use <- c()
for(i in files){
  raw.data.i <- read.csv(paste0("Raw_Data/", i))
  use.year.i <- which(colnames(raw.data.i) == "MMWR.Year")
  use.week.i <- which(colnames(raw.data.i) == "MMWR.Week")
  use.reporting.i <- which(colnames(raw.data.i) == "Reporting.Area")
  if(length(use.year.i) != 1 | length(use.week.i) != 1 | length(use.reporting.i) != 1){
    stop()
  }
  
  years.i <- unique(raw.data.i[,use.year.i])
  if(length(years.i) != 1){
    stop()
  }
  weeks.i <- unique(raw.data.i[,use.week.i])
  
  years_weeks <- c(years_weeks, rep(years.i, length(weeks.i)))
  weeks_weeks <- c(weeks_weeks, weeks.i)
  
  colnames_use <- c(colnames_use, colnames(raw.data.i)[4])
  
  raw_data[[i]] <- raw.data.i 
}

merge_data <- matrix(NA, ncol = 53, nrow = length(years_weeks))
colnames(merge_data) <- c("MMWR_Year", "MMWR_Week", "Date", toupper(state.name))
merge_data <- data.frame(merge_data)

merge_data$MMWR_Year <- years_weeks
merge_data$MMWR_Week <- weeks_weeks
merge_data$Date <- MMWRweek2Date(MMWRyear = years_weeks, MMWRweek = weeks_weeks)
merge_data$Date <- strptime(merge_data$Date, format = "%Y-%m-%d")
merge_data <- merge_data[order(merge_data$Date, decreasing = FALSE), ]

#3. Merge data
for(i in 1:length(files)){
  data.i <- raw_data[[i]]
  for(s in 4:ncol(merge_data)){
    use.s <- which(data.i$Reporting.Area == colnames(merge_data)[s])
    use.y <- which(merge_data$MMWR_Year == years[i])
    mt.week <- match(merge_data$MMWR_Week[use.y], data.i$MMWR.Week[use.s])
    merge_data[use.y,s] <- data.i[use.s,colnames_use[i]][mt.week]
  }
}

#################
#State pop sizes#
#################
#data from https://www.census.gov/data/tables/2017/demo/popest/state-total.html

state_pops <- read.csv("State_Pops/nst-est2017-alldata.csv", stringsAsFactors = FALSE)

pop_cols <- grep(pattern = "POPESTIMATE", x = colnames(state_pops), ignore.case = TRUE)
years_pops <- unlist(lapply(colnames(state_pops)[pop_cols], function(x) unlist(strsplit(x, split = "POPESTIMATE"))[2]))
years_pops_use <- which(years_pops %in% years)
pop_cols_use <- pop_cols[years_pops_use]

pop_data <- matrix(NA, ncol = ncol(merge_data), nrow = nrow(merge_data))
colnames(pop_data) <- colnames(merge_data)
pop_data <- data.frame(pop_data)
pop_data$MMWR_Year <- merge_data$MMWR_Year
pop_data$MMWR_Week <- merge_data$MMWR_Week
pop_data$Date <- merge_data$Date

tab_years <- table(pop_data$MMWR_Year)

state_pop_names_for_mt <- toupper(state_pops$NAME)
state_pop_names_for_mt <- gsub(pattern = "[ ]", replacement = ".", x = state_pop_names_for_mt)
for(i in 4:ncol(pop_data)){
  use.i <- which(state_pop_names_for_mt == colnames(pop_data)[i])
  if(length(use.i) != 1){
    stop()
  }
  data.i <- as.numeric(state_pops[use.i,pop_cols_use])
  x <- 1:4
  mod.i <- lm(data.i~x)
  pred.18.i <- predict(mod.i, newdata = data.frame(x = 5))
  data.i <- c(data.i, pred.18.i)
  pop_data[,i] <- rep(data.i, tab_years)
}

national_cases <- rowSums(merge_data[,-c(1:3)], na.rm = TRUE)
natioanl_pops <- rowSums(pop_data[,-c(1:3)], na.rm = TRUE)
national_rate <- national_cases/natioanl_pops
national_rate_out <- data.frame(merge_data[,1:3], national_rate)
colnames(national_rate_out) <- c(colnames(merge_data)[1:3], "prop_invasive_pneumococcal_disease")

###########
#Save Data#
###########
if(write_new == TRUE){
  file_name_merge <- paste0(time_stamp, "_", min(years), "_", max(years), "_NNDSS_Table_II_invasive_pneumococcal_disease.csv")
  write.csv(merge_data, file = file_name_merge, row.names = FALSE, quote = FALSE)
  
  file_name_pop <- paste0(time_stamp, "_", min(years), "_", max(years), "_US_census_state_pop_data.csv")
  write.csv(pop_data, file = file_name_pop, row.names = FALSE, quote = FALSE)
  
  file_name_rate <- paste0(time_stamp, "_", min(years), "_", max(years), "_prop_NNDSS_Table_II_invasive_pneumococcal_disease.csv")
  write.csv(national_rate_out, file = file_name_rate, row.names = FALSE, quote = FALSE)
}