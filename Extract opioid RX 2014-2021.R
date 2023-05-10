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



### STEP1: FIND OPIOID PRODUCT ID

# Lookup lists
lookup_drug_class = data.table(read_excel("\\\\rfawin.partners.org\\mgh-itafisma\\IQVIA\\lookup_lists.xlsx", sheet = "opioids_pain"))

# Filter down to opioid prescription only from the lookup lists
list_opioid_mkted_prod_name = lookup_drug_class$MKTED_PROD_NM

# Product
d_product = fread("\\\\rfawin.partners.org\\mgh-itafisma\\IQVIA\\Product.gz")

# Get the product IDs from the product data where the MKTED_PROD_NM matches
opioid_productID_NM = d_product[MKTED_PROD_NM %in% list_opioid_mkted_prod_name,c("PRODUCT_ID", "MKTED_PROD_NM","USC_CD", "USC_DESC", "STRNT_DESC")]

list_opioid_PRODUCT_ID = opioid_productID_NM$PRODUCT_ID
# 6197 unique product_id



### STEP2: FILTER RX DATA WITH PRODUCT ID

# Selectively read in Rx data

filenum_rx = c(seq(2014,2021,1))
datanames_rx = paste("FACT_RX_", filenum_rx, sep='')

for(i in 1:length(datanames_rx))
{
  assign (
    datanames_rx[i], 
    fread(file = paste("\\\\rfawin.partners.org\\mgh-itafisma\\IQVIA\\Rx_data\\", datanames_rx[i], sep=""),
          select = c("SVC_DT", "PATIENT_ID","PRODUCT_ID","DAYS_SUPPLY_CNT","DSPNSD_QTY"), quote="")
    [PRODUCT_ID %in% list_opioid_PRODUCT_ID]
  )
} 


FACT_RX = rbind (FACT_RX_2014, FACT_RX_2015, FACT_RX_2016, FACT_RX_2017,
                 FACT_RX_2018, FACT_RX_2019, FACT_RX_2020, FACT_RX_2021)

RX_product = merge(FACT_RX,opioid_productID_NM, by="PRODUCT_ID",all.x=TRUE )

table(RX_product$USC_DESC)


fwrite(RX_product, "C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\RX_product_opioids_2014_2021.csv")
# 351,517,666 LRx


### STEP3: MERGE RX, PRODUCT, PATIENT, MME data, and SAMPLE RESTRICTIONS

df_opioids = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\RX_product_opioids_2014_2021.csv")
# 351,517,666 LRx
list_patient = unique(df_opioids$PATIENT_ID)
# 48,311,539 unique patients

# Product info
df_product = fread("\\\\rfawin.partners.org\\mgh-itafisma\\IQVIA\\Product.gz", select=c("PRODUCT_ID", "NDC_CD", "DOSAGE_FORM_NM"))
df_prescription_product = merge(df_opioids,df_product, by="PRODUCT_ID", all.x=TRUE)

# Patient info
df_patient = fread("\\\\rfawin.partners.org\\mgh-itafisma\\IQVIA\\Patient.gz")[PATIENT_ID %in% list_patient]
df_rx_patient = merge(df_prescription_product,df_patient,by="PATIENT_ID", all.x=TRUE)
# 351,517,666 records

# fwrite(df_rx_patient, "C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\df_rx_patient.csv")
# 351,517,666 LRx

# ls()
# rm("df_opioids" , "df_patient", "df_prescription_product", "df_product", "list_patient")


# Sample restriction 1: Keep opioid flag 1
df_rx_patient_1 = df_rx_patient[OPIOID_FLAG==1,]
# 351,403,952 records

# Sample restriction 2: Remove injectables and crude/bulk
df_rx_patient_2 = df_rx_patient_1[!(USC_DESC %in% c("MORPHINE/OPIUM, INJECTABLE","SYNTH NARCOTIC, INJECTABLE","CRUDE/BULK MEDICINAL")),]
# 351,049,715 records

# Sample restriction 3: Based on USC codes
table(df_rx_patient_2$USC_DESC)
table(df_rx_patient_2$USC_CD)
# 2214      2222        2232        9150       34210     34220     34240     34290      34380 
# 70861942  24641337    241770343   11621631   1454     1658033    13450     481522       3 

rm("df_rx_patient" , "df_rx_patient_1")

df_rx_patient_3 = df_rx_patient_2[USC_CD %in% c(2214,2222,2232),]
# 337,273,622 obs

# fwrite(df_rx_patient_3, "C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\df_rx_patient_3.csv")


# Check NDC code;No missing value
check = df_rx_patient_3[,c("PATIENT_ID","NDC_CD")]
sum(is.na(check))

# ls()
# rm("df_rx_patient" , "df_rx_patient_1", "df_rx_patient_2", "check")

### STEP4: READ IN MME CONVERSION FILE

convert = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\MME conversion table.csv")

# Merge the Rx file with the converstion table
df_rx_mme = merge(df_rx_patient_3, convert, by.x = "NDC_CD", by.y = "NDC_Numeric", all.x=TRUE)
# 337,273,622 obs

# Remove extra variables after some checking
df_rx_mme2 = df_rx_mme[,-c("MOUD_FLAG","NALOXONE_FLAG","OPIOID_FLAG",
                           "NDC","PRODNME","GENNME","Master_Form","Class","Drug","LongShortActing","DEAClassCode")]


# The records could not be linked
missconvert = df_rx_mme2[is.na(MME_Conversion_Factor),]
# 1,592,351 (0.47%) records cannot be linked

missprod = missconvert[,c("MKTED_PROD_NM", "STRNT_DESC","Strength_Per_Unit", "UOM", "MME_Conversion_Factor")]
missprod = missprod[!duplicated(missprod)]
# 87 combinations of drug name and strength 

# Save the file for manual check and fix
# fwrite(missprod, "C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Missproduct_toFill_2014_2021.csv")


# The records that can be linked
nomissconvert = df_rx_mme2[!is.na(MME_Conversion_Factor),]
# 335,681,271 obs

# Use this as reference for filling up the missing value
nomissprod = nomissconvert[,c("MKTED_PROD_NM", "STRNT_DESC", "Strength_Per_Unit", "UOM", "MME_Conversion_Factor")]
nomissprod = nomissprod[!duplicated(nomissprod)]
# 384 cominations of drug name and strength


# Use the manual fixed data file to fill up the missing value
tofill = fread("C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\Missproduct_Conversion_Factor_2014_2021.csv")

missconvert2 = missconvert[,-c("Strength_Per_Unit","UOM","MME_Conversion_Factor")]

missconvert_fill = merge (missconvert2, tofill, by=c("MKTED_PROD_NM","STRNT_DESC"), all.x=TRUE)

# Comine the nomiss one and the filled up missed rows
RX_convert_full = rbind(nomissconvert, missconvert_fill)
setorder(RX_convert_full, PATIENT_ID, SVC_DT)
# 337,273,622 obs, match the initial dataset

# Check the conversion factor
summary(RX_convert_full$MME_Conversion_Factor)
# No missing value

# Save the full dataset
fwrite(RX_convert_full, "C:\\Users\\hd820\\OneDrive - Mass General Brigham\\Data Repository\\IQVIA analysis\\[P2] Upstream model\\MME calculation\\Data\\RX_Patient_MME_2014to2021.csv")
