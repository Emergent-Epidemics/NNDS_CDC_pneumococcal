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

###########
#Save Data#
###########
if(write_new == TRUE){
  file_name <- paste0(time_stamp, "_", min(years), "_", max(years), "_NNDSS_Table_II_invasive_pneumococcal_disease.csv")
  write.csv(merge_data, file = file_name, row.names = FALSE, quote = FALSE)
}