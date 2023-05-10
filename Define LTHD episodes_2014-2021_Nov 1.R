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


# Read in all the episode data
episode=fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA_analysis\\[P2] Upstream model\\MME_calculation\\Data\\Split episode\\Opioid_episode_all_2022Jun15.csv")
# 119,494,139 episodes in total

ppp_rev2 = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA_analysis\\[P2] Upstream model\\MME_calculation\\Data\\Data_2014_2021\\ppp_rev2.csv")
# 336,892,984 obs

# Patient info
df_patient = fread("\\\\rfawin.partners.org\\mgh-itafisma\\IQVIA\\Patient.gz",
                   select=c("PATIENT_ID","PAT_BRTH_YR_GRP","PAT_BRTH_YR_NBR","PAT_GENDER_CD","PAT_ETHNICITY"))

# Restrict on at least 90 days
episode2 = episode[episode.duration>=90,]
# 9,813,672 episodes had more than 90 days


epi.pat = merge(episode2,df_patient,by="PATIENT_ID", all.x=TRUE)

table(epi.pat$PAT_BRTH_YR_GRP,useNA="always")
table(epi.pat$PAT_GENDER_CD,useNA="always")
table(epi.pat$PAT_ETHNICITY,useNA="always")

# Restrict the sample to complete age, sex, ethnicity information
epi.pat2 = epi.pat[PAT_BRTH_YR_GRP %in% c("0-11","12-20","21-34","35-54","55+"),]
# 9,808,737 episodes after removing missing in age

epi.pat3 = epi.pat2[PAT_GENDER_CD %in% c("F","M"),]
# 9,784,658 episodes after removing missing in sex

epi.pat4 = epi.pat3[PAT_ETHNICITY %in% c("ASIAN - OTHER","BLACK","HISPANIC","WHITE"),]
# 4,426,023 episodes after removing missing for race

table(epi.pat4$PAT_BRTH_YR_GRP,useNA="always")
table(epi.pat4$PAT_GENDER_CD,useNA="always")
table(epi.pat4$PAT_ETHNICITY,useNA="always")

# The list of patients
listpatient = unique(epi.pat4$PATIENT_ID) 
# 2,286,751 patients

# Get only relevant LRx records
RX_epi = ppp_rev2[PATIENT_ID %in% listpatient,]
# 107,450,104 RXs

RX_epi2 = RX_epi[,c("PATIENT_ID","PRODUCT_ID","SVC_DT","DAYS_SUPPLY_CNT","DSPNSD_QTY","MKTED_PROD_NM",
                    "STRNT_DESC","Strength_Per_Unit","MME_Conversion_Factor")]

RX_epi2$supply.start = as.IDate(ymd(RX_epi2$SVC_DT), "%Y-%m-%d")
RX_epi2$supply.end = RX_epi2$supply.start + RX_epi2$DAYS_SUPPLY_CNT-1

rm(episode,episode2,epi.pat,epi.pat2,epi.pat3,df_patient,ppp_rev2)
gc()



# # Let's try with a small sample
# testsample = epi.pat4[c(1:50),c("PATIENT_ID","episode.ID","episode.start","episode.duration","episode.end","episode.end.date")]
# 
# testID = unique(testsample$PATIENT_ID)
# 
# # Merge with the RX data
# testRX = RX_epi2[PATIENT_ID %in% testID,]

RX_episode = merge(epi.pat4,RX_epi2, by="PATIENT_ID", allow.cartesian=TRUE)
RX_episode_red = RX_episode[!(supply.end<episode.start | supply.start>episode.end.date),]


# Calculating DMME
setorder(RX_episode_red, PATIENT_ID, episode.ID, supply.start)

RX_episode_red$DMME = (RX_episode_red$DSPNSD_QTY/RX_episode_red$DAYS_SUPPLY_CNT)* RX_episode_red$Strength_Per_Unit * RX_episode_red$MME_Conversion_Factor

# Double check the number of episodes
check = RX_episode_red[,c("PATIENT_ID","episode.ID")]
nrow(check[!duplicated(check)])
#4,426,023 episodes matched with the sample restriction results

#####################################

# Take out all segments with DMME>=90
high90 = RX_episode_red[DMME>=90,c("PATIENT_ID","episode.ID","episode.start","episode.duration","episode.end","episode.end.date",
                                       "DAYS_SUPPLY_CNT","supply.start","supply.end","DMME")]

# Double check the number of episodes had at least some days with dose >=90 DMME
check = high90[,c("PATIENT_ID","episode.ID")]
nrow(check[!duplicated(check)])
# 848,769 episodes

# Get the last supply end date
high90[, last_supply_end_date := shift(supply.end, n=1, type = 'lag', fill = NA), by = list(PATIENT_ID,episode.ID)] 
high90[,gap := supply.start-last_supply_end_date]

high90$flag=1
high90[gap>1 | is.na(gap),]$flag=0

high90$continue = cumsum(high90$flag==0)

# Get the duration for continued high dose episodes 
setorder(high90, PATIENT_ID, episode.ID, supply.start,continue)

# Should have this many continued high dose segments
max(high90$continue)
# 5,058,763

# Earliest and last 
duration = high90[,c("PATIENT_ID","episode.ID","episode.start","episode.end","supply.start","supply.end","continue")]
first_high = duration %>% group_by(continue) %>% arrange(PATIENT_ID,episode.ID,continue) %>% filter(row_number()==1)
first_high = first_high[,c("PATIENT_ID","episode.ID","episode.start","episode.end","continue","supply.start")]

last_high = duration %>% group_by(continue) %>% arrange(PATIENT_ID,episode.ID,continue) %>% filter(row_number()==n())
last_high = last_high[,c("continue", "supply.end")]

duration_com = merge(first_high,last_high,by="continue")
duration_com$duration = difftime(duration_com$supply.end, duration_com$supply.start, units="days")+1
duration_com$day=as.numeric(duration_com$duration)

# MAX length of continued days with >=90 DMME for each episode
setorder(duration_com, PATIENT_ID, episode.ID, duration)

# Keep the last record as the longest days
longest_perepisode = duration_com %>% group_by(PATIENT_ID,episode.ID) %>% arrange(PATIENT_ID, episode.ID) %>% filter(row_number()==n())
# 848,769 episodes
longest_perepisode = longest_perepisode[,c("PATIENT_ID","episode.ID","episode.start","episode.end","day")]


# Calculate the total number of days with high dose
tot = aggregate(x=duration_com$day,
                          by=list(duration_com$PATIENT_ID,duration_com$episode.ID),
                          FUN=function(x) c(Total.HDdays=sum(x)))
colnames(tot)=c("PATIENT_ID","episode.ID","Total.HDdays")
setorder(tot, PATIENT_ID,episode.ID)
# 848,769 episodes

longest.sum.perepisode = merge(longest_perepisode,tot,by=c("PATIENT_ID","episode.ID"))
setorder(longest.sum.perepisode, PATIENT_ID,episode.ID)
colnames(longest.sum.perepisode)=c("PATIENT_ID","episode.ID","episode.start","episode.end2","Consecutive.HDdays","Total.HDdays")

summary(longest_perepisode$day)
quantile(longest_perepisode$day, c(.25, .50, .75,.90,.95)) 
hist(longest_perepisode$day,xlim=c(0,460))

nrow(longest.sum.perepisode[longest.sum.perepisode$day>=81,]) # 278,490
nrow(longest.sum.perepisode[longest.sum.perepisode$day>=90,]) # 245,789
nrow(longest.sum.perepisode[longest.sum.perepisode$day>=14,]) # 551,458

#fwrite(longest.sum.perepisode,"C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA_analysis\\[P2] Upstream model\\Racial_Disparities_Geo\\Data\\longest and total high dose per episode.csv")

# # Plot 
par(mar=c(7, 7, 7, 7) + 0.1)
hist(longest_perepisode[longest_perepisode$day<=400,]$day, breaks = seq(from=0, to=400, by=1), xlim=c(0, 410),
      xlab="", ylab="", main=NA, font=2,las=1)
par(new=TRUE)
plot(ecdf(longest_perepisode[longest_perepisode$day<=400,]$day),col="red", xlim=c(0, 410),do.points=FALSE,
      xlab="Max consecutive days per episode with DMME>=90", ylab="", main=NA,xaxt='n', yaxt='n',lwd=2)
axis(4, ylim=c(0,1), col="red",col.axis="red",las=1, font=2)


par(mar=c(7, 7, 7, 7) + 0.1)
hist(longest_perepisode[longest_perepisode$day<=100,]$day, breaks = seq(from=0, to=100, by=1), xlim=c(0, 110),
     xlab="", ylab="", main=NA, font=2,las=1)


##########################################################
# Keep the records from the original time segments datasets

part1 = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA_analysis\\[P2] Upstream model\\MME_calculation\\Data\\Split time segments\\Time_segments_all_part1.csv")
part2 = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA_analysis\\[P2] Upstream model\\MME_calculation\\Data\\Split time segments\\Time_segments_all_part2.csv")

part1.red = part1[episode.duration>=90,]
part2.red = part2[episode.duration>=90,]

seg = rbind(part1.red,part2.red)
setorder(seg,PATIENT_ID, episode.ID,seq)


select14 = longest_perepisode[longest_perepisode$day>=14,c("PATIENT_ID","episode.ID")]
# 551,458
select90 = longest_perepisode[longest_perepisode$day>=90,c("PATIENT_ID","episode.ID")]
# 245,789

# Restrict the segment data to episodes with 90+ days with 90+ DMME
L90H90days = merge(seg,select90,by=c("PATIENT_ID","episode.ID"))

### Double check
check = L90H90days[,c("PATIENT_ID","episode.ID")]
nrow(check[!duplicated(check)])
# 245,789
table(L90H90days$PAT_BRTH_YR_GRP)
table(L90H90days$PAT_GENDER_CD)
table(L90H90days$PAT_ETHNICITY)


# Restrict the segment data to episodes with 14+ days with 90+ DMME
L90H14days = merge(seg,select14,by=c("PATIENT_ID","episode.ID"))

### Double check
check = L90H14days[,c("PATIENT_ID","episode.ID")]
nrow(check[!duplicated(check)])
# 551,458
table(L90H14days$PAT_BRTH_YR_GRP)
table(L90H14days$PAT_GENDER_CD)
table(L90H14days$PAT_ETHNICITY)


fwrite(L90H90days,"C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA_analysis\\[P2] Upstream model\\Racial_Disparities_Geo\\Data\\L90H90days.csv")
fwrite(L90H14days,"C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA_analysis\\[P2] Upstream model\\Racial_Disparities_Geo\\Data\\L90H14days.csv")


