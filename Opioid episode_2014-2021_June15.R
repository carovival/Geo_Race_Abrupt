#install.packages("arsenal")
library(readxl) 
library(data.table)
library(lubridate) 
library(tidyverse)
library(survival)
library(AdhereR)
library(ggpubr)
library(parallel)
library(arsenal)

memory.limit(size=1000000000)



################################################################################

# Read in PPP data
df_ppp = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA_analysis\\[P2] Upstream model\\MME_calculation\\Data\\Data_2014_2021\\RX_Patient_MME_2014to2021.csv")
# 337,273,622 obs

listpatient = unique(df_ppp$PATIENT_ID)
# 47,911,161 patients



# Clean up the RX file

# Check distribution of DAYS_SUPPLY_CNT
summary(df_ppp$DAYS_SUPPLY_CNT)

df_ppp[!is.na(DAYS_SUPPLY_CNT)] %>% dplyr::summarize(
  quants = quantile(DAYS_SUPPLY_CNT, probs = c(0.5, 0.75, 0.90, 0.95, 0.99, 0.995,0.995)))

# quants
# 1     15
# 2     30
# 3     30
# 4     30
# 5     30
# 6     60
# 7     60

nrow(df_ppp[DAYS_SUPPLY_CNT==999,]) #104
nrow(df_ppp[DAYS_SUPPLY_CNT==0,]) #4745

# Remove the records with 999, NA (14), and 0 days
ppp_rev = df_ppp[!is.na(DAYS_SUPPLY_CNT) & DAYS_SUPPLY_CNT!=999 & DAYS_SUPPLY_CNT!=0,]
# 337,268,759 observations

# check duplicated records
ppp_rev2 = ppp_rev[!duplicated(ppp_rev)]
# 336,892,984 obs after removing duplicates
nrow(ppp_rev2)/nrow(ppp_rev)*100
# 0.12% duplicated removed

#summary(ppp_rev2$DAYS_SUPPLY_CNT)
fwrite(ppp_rev2, "C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\ppp_rev2.csv")



################################################################################
################################################################################

## DEFINE TREATMENT EPISODE 
ppp_rev2 = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\ppp_rev2.csv",
                 select=c("PATIENT_ID","SVC_DT","DAYS_SUPPLY_CNT"))
# 336,892,984 obs

listpatient = unique(ppp_rev2$PATIENT_ID)
# 47,910,093 patients


# Let's try every 1 million patients and save their LRx data
n=1000000
nr=length(listpatient)
dat=split(listpatient, rep(1:ceiling(nr/n), each=n, length.out=nr))

for (i in 1:(ceiling(nr/n)))
{
  LRx = ppp_rev2[PATIENT_ID %in% data.frame(dat[i])[,1],]
  fwrite(LRx, 
         file=paste("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split LRx\\LRx",i,".csv",sep=""))
  
}


# Read in data and double check
# filenum_rx = c(seq(1,48,1))
# datanames_rx = paste("LRx", filenum_rx, sep='')
# 
# for(i in 1:length(datanames_rx))
# {
#   assign (datanames_rx[i], 
#     fread(file = paste("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split LRx\\", datanames_rx[i],".csv", sep="")
#   ))
# } 


# Functions to run opioid episodes
episode_function = function(dataset,idx){
  
  # Convert service to date
  dataset$SVC_DT_date = format(ymd(dataset$SVC_DT), "%m/%d/%Y")
  
  opioid_episode=compute.treatment.episodes(
    dataset,
    ID.colname="PATIENT_ID",
    event.date.colname="SVC_DT_date",
    event.duration.colname="DAYS_SUPPLY_CNT",
    carryover.within.obs.window = FALSE, 
    carry.only.for.same.medication = FALSE,
    medication.change.means.new.treatment.episode = FALSE,
    maximum.permissible.gap = 30, 
    maximum.permissible.gap.unit = "days", 
    followup.window.start = 0, 
    followup.window.start.unit = "days",
    # The max duration would be 3650 days
    followup.window.duration = 365 * 10,  
    followup.window.duration.unit = "days",
    date.format = "%m/%d/%Y")
  
  # Save the episode results
  fwrite(opioid_episode, 
         file=paste("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split episode\\Opioid_episode_",idx,"_2022Jun15.csv",sep=""))
  
}


# Run the episode function to build episode for each 1M patients
# Run in batches
for(i in 1:10)
{
  assign ("dataset", 
          fread(file = paste("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split LRx\\LRx", i,".csv", sep="")))
  
  episode_function(dataset,i)
} 

# check=fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split episode\\Opioid_episode_1_2022Jun15.csv")
# check[1:10,]
# 
# check2 = dataset[1:100,]
################################################################################

# Read in all the episode results

filenum = c(seq(1, 48, 1))
filenames_episode = paste("Opioid_episode_", filenum, "_2022Jun15.csv", sep="")
datanames = paste("episode", filenum, sep="")

for(i in 1:48)
{
  assign (datanames[i], 
          fread(file = paste("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split episode\\", 
                             filenames_episode[i], sep="") ))
} 

rm(datanames, filenames_episode, filenum,i)

# Combine all the episode data
df_list = mget(ls(pattern = "episode*"))

episode = dplyr::bind_rows(df_list)
# 119,494,139 episodes in total

setorder(episode,PATIENT_ID, episode.ID)

length(unique(episode$PATIENT_ID))
# This should be 47,910,093 patients, checked!

# summary(episode$episode.duration)
# episode[!is.na(episode.duration)] %>% dplyr::summarize(
#   quants = quantile(episode.duration, probs = c(0.5, 0.75, 0.90, 0.95, 1)))

# Compute the real ended dates
episode$episode.end.date = episode$episode.end-1

# Save the episodes data
fwrite(episode, "C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split episode\\Opioid_episode_all_2022Jun15.csv")



################################################################################
################################################################################

# Chop the period to every 30 days

episode=fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split episode\\Opioid_episode_all_2022Jun15.csv")
# 119,494,139 episodes in total


# Let's try to split every 1 million records of episode
n=1000000
nr=nrow(episode)
dat.episode0=split(episode, rep(1:ceiling(nr/n), each=n, length.out=nr))

for (i in 1:(ceiling(nr/n)))
{
  dat.episode = data.frame(dat.episode0[i])
  colnames(dat.episode)=c("PATIENT_ID","episode.ID","episode.start","end.episode.gap.days","episode.duration","episode.end","episode.end.date")
  fwrite(dat.episode, 
         file=paste("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split time segments\\Split_episode_",i,".csv",sep=""))
  
}



# This is the function that can split the episode to every 30 days
split_function = function(dat){
  
  dat$num_int = ceiling(dat$episode.duration / 30)
  
  f_rep = function(rowindex) {rep(rowindex, each = dat[rowindex, ]$num_int)}
  g_rep = function(rowindex) {1:dat[rowindex, ]$num_int}
  
  rowindex_rep = unlist(Map(f_rep, 1:nrow(dat)))
  
  dat2 = dat[rowindex_rep, ]
  dat2$seq = unlist(Map(g_rep, 1:nrow(dat)))
  
  dat2$split_end_date = dat2$episode.end.date-(dat2$seq-1)*30
  dat2$split_start_date = dat2$split_end_date-29
  dat2[which(dat2$seq == dat2$num_int), ]$split_start_date = dat2[which(dat2$seq == dat2$num_int), ]$episode.start
  dat2$split_duration = as.integer(dat2$split_end_date - dat2$split_start_date)+1
  
  return(dat2)
}


ppp_rev2 = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Data 2014_2021\\ppp_rev2.csv")
# 336,892,984 obs

# Patient info
df_patient = fread("\\\\rfawin.partners.org\\mgh-itafisma\\IQVIA\\Patient.gz")


# Run in batches
for(i in 61:120)
{
  assign ("dat", 
          fread(file = paste("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split time segments\\Split_episode_",i,".csv", sep="")))
  
  split_episode=split_function(dat)
  
  # The list of patients
  listpatient = unique(split_episode$PATIENT_ID) 
  
  # Get only relevant LRx records
  RX_epi = ppp_rev2[PATIENT_ID %in% listpatient]
  
  RX_epi$supply.start = as.IDate(ymd(RX_epi$SVC_DT), "%Y-%m-%d")
  RX_epi$supply.end = RX_epi$supply.start + RX_epi$DAYS_SUPPLY_CNT-1
  
  # Merge with split data 
  RX_episode = merge (split_episode, RX_epi, by="PATIENT_ID",allow.cartesian=TRUE)
  RX_episode_red = RX_episode[!(supply.end<split_start_date | supply.start>split_end_date),]
  
  # There are time-segments that are totally removed because of the 30 days gap, fill them back up
  full = split_episode[,c("PATIENT_ID", "episode.ID", "seq")]
  full = full[!duplicated(full)]
  
  remained = RX_episode_red[,c("PATIENT_ID", "episode.ID", "seq")]
  remained = remained[!duplicated(remained)]
  
  removed = anti_join(full, remained)
  
  tofill = merge(split_episode, removed)
  
  # Add back the rows
  RX_episode_fill = rbind(RX_episode_red, tofill, fill=TRUE)
  
  # Since several rows are filled up and patient info will be missing, therefore re-do that part
  RX_episode_fill = RX_episode_fill[,-c("PAT_BRTH_YR_GRP","PAT_BRTH_YR_NBR","PAT_GENDER_CD","PAT_ETHNICITY")]
  RX_epiall = merge(RX_episode_fill,df_patient,by="PATIENT_ID", all.x=TRUE)
  
  
  # Count how many days contribution
  RX_epiall$contributed.time=0
  RX_epiall[supply.start<split_start_date & supply.end<=split_end_date,]$contributed.time =
    RX_epiall[supply.start<split_start_date & supply.end<=split_end_date,]$supply.end - 
    RX_epiall[supply.start<split_start_date & supply.end<=split_end_date,]$split_start_date +1
  
  RX_epiall[supply.start>=split_start_date & supply.end<=split_end_date,]$contributed.time = 
    RX_epiall[supply.start>=split_start_date & supply.end<=split_end_date,]$DAYS_SUPPLY_CNT
  
  RX_epiall[supply.end>split_end_date & supply.start>=split_start_date,]$contributed.time = 
    RX_epiall[supply.end>split_end_date & supply.start>=split_start_date,]$split_end_date - 
    RX_epiall[supply.end>split_end_date & supply.start>=split_start_date,]$supply.start+1
  
  RX_epiall[supply.start<split_start_date & supply.end>split_end_date,]$contributed.time =
    RX_epiall[supply.start<split_start_date & supply.end>split_end_date,]$split_duration 
  
  # Calculating MME and DMME
  setorder (RX_epiall, PATIENT_ID, episode.ID, seq, supply.start)
  
  # Calculate MME
  # This is the data file that contains all PATIENT, episode, time segments, and relevant opioid prescriptions
  RX_epiall$contributed.MME = (RX_epiall$DSPNSD_QTY/RX_epiall$DAYS_SUPPLY_CNT)* RX_epiall$contributed.time * RX_epiall$Strength_Per_Unit * RX_epiall$MME_Conversion_Factor
  
  # Now aggregate by PATIENT_ID, episode ID, and 30-day time segment
  # Aggregate the MME for each time segment
  MME_time = aggregate(x=RX_epiall$contributed.time,
                       by=list(RX_epiall$PATIENT_ID,RX_epiall$episode.ID,RX_epiall$seq),
                       FUN=function(x) c(contributed.sum.time=sum(x)))
  colnames(MME_time)=c("PATIENT_ID", "episode.ID", "seq", "contributed.sum.time")
  
  MME_MME = aggregate(x=RX_epiall$contributed.MME,
                      by=list(RX_epiall$PATIENT_ID,RX_epiall$episode.ID,RX_epiall$seq),
                      FUN=function(x) c(contributed.sum.MME=sum(x)))
  colnames(MME_MME)=c("PATIENT_ID", "episode.ID", "seq", "contributed.sum.MME")
  
  summary_time_mme = merge(split_episode, MME_time, by=c("PATIENT_ID", "episode.ID", "seq"))
  summary_time_mme = merge(summary_time_mme, MME_MME, by=c("PATIENT_ID", "episode.ID", "seq"))
  
  # Calculate DMME
  summary_time_mme$DMME.time.segment = summary_time_mme$contributed.sum.MME/summary_time_mme$split_duration
  
  # Merge patient info
  summary_time_mme = merge(summary_time_mme,df_patient,by="PATIENT_ID", all.x=TRUE)
  
  fwrite(summary_time_mme, 
         file=paste("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split time segments\\Time_segments_",i,".csv", sep=""))
  
} 




# Double check
# check = time_seg %>% group_by(PATIENT_ID, episode.ID) %>% filter(row_number()==1)
# sum(check$num_int)
# nrow(time_seg)

# summary(RX_epiall$contributed.time)
# Should have 0 records here:
# RX_epiall[contributed.time>DAYS_SUPPLY_CNT,]
# summary(RX_epiall$contributed.MME)
# summary_time_mme[!is.na(DMME.time.segment)] %>% dplyr::summarize(
#  quants = quantile(DMME.time.segment, probs = c(0.5, 0.75, 0.90, 0.95, 1)))



# Read in all the results

filenum = c(seq(1, 120, 1))
filenames_episode = paste("Time_segments_", filenum, ".csv", sep="")
datanames = paste("Segment_", filenum, sep="")

for(i in 61:120)
{
  assign (datanames[i], 
          fread(file = paste("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split time segments\\", 
                             filenames_episode[i], sep="") ))
} 

rm(datanames, filenames_episode, filenum,i)

# Combine all the episode data
df_list = mget(ls(pattern = "Segment*"))

Segmentall = dplyr::bind_rows(df_list)
# 140,580,414 obs in part1
# 139,423,733 obs in part2
# 280,004,147 obs in total

setorder(Segmentall,PATIENT_ID, episode.ID,seq)

length(unique(Segmentall$PATIENT_ID))
# 24,052,450 in part1
# 23,857,644 in part2
# This should be 47,910,093 patients, checked!


# Save the full time segment data
# fwrite(Segmentall, "C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split time segments\\Time_segments_all_part2.csv")


part1 = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split time segments\\Time_segments_all_part1.csv")
part2 = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split time segments\\Time_segments_all_part2.csv")

# all = rbind(part1, part2)
# 
# setorder(all,PATIENT_ID, episode.ID,seq)
# 
# length(unique(all$PATIENT_ID))
# # 47,910,093 unique patients
# 
# fwrite(all, "C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Split time segments\\Time_segments_all.csv")
# 
