########################
########################
### DATA 2023 MARIDI ###
########################
########################
################
## LOAD PACKAGES
################
library(readxl)
library(dplyr)
library(ggplot2)
###############
## READ IN DATA
###############
DATA <- read_excel("Seroprevalence_dataset.xlsx")
View(DATA)

###################
## DATA EXPLORATION
###################

class(DATA) # Checking upload
list(DATA) #First look at the data
colnames(DATA) #oversight variables
summary(DATA) #A summary of the data set, each variable is explained in the accompanying codebook.
#Median age: 6 (IQR: 3-8)

###Preparing data
#Making factorial factors:
names <- c('village' ,'sex', 'level_of_education', 'level_of_primary', 'ivermectin', 'dermititis', 
           'epilepsy', 'ivermectin_questionnaire', 'ov16rdt') #Selection of variables to be turned into factors
DATA[,names] <- lapply(DATA[,names] , factor) #Changing them all together to factors using lapply

str(DATA) #Checking changes

###Exploring variables
##Dermatitis
table(DATA$dermititis) #Showing proportions
# 0   1 
#159  93 

#Percentage dermatitis
93/(159+93) #36.9%

##Epilepsy
table(DATA$epilepsy)
#  0   1 
#249   4 

#Percentage epilepsy
4/(249+4) #1.6%

##Sex
table(DATA$sex)
# 1   2   (Female/Male)
#137 116 

#Percentage males
(116)/(116+137) #45.8%

##Age
table(DATA$age)
table(DATA$age,DATA$ov16rdt)
median(DATA$age)
summary(DATA$age)


##################
## OV16 PREVALENCE
##################
#####Overall prevalence OV16
OV16_ALL <- as.data.frame(table(DATA$ov16rdt)) #Summarizing the OV16 variable in a data frame
prevalence_all <- (prod(OV16_ALL[1:1, "Freq"])/(prod(OV16_ALL[1:1, "Freq"])+ prod(OV16_ALL[2:2, "Freq"])))*100 #Calculating the prevalence using the data frame (OV16 positive/All tested cases times 100)
print(prevalence_all) #Printing the calculated prevalence percentage
#ALL: 30.64516%
prop.test(prod(OV16_ALL[1:1, "Freq"]),(prod(OV16_ALL[1:1, "Freq"])+ prod(OV16_ALL[2:2, "Freq"]))) #Performing a proportion test based on chi-square and fisher tests to determine the 95% confidence interval
#CI: 25.1-36.9%

#####Prevalence per village OV16
villages <- c(levels(DATA$village)) #creating a separate variable based on the levels of the 'village' variable to use in a looped based prevalence calculation
for (i in villages ) #For loop based on the previously created separate variable
{tmp <- subset.data.frame(DATA, village == i) #Summarizing the initial data frame based on a single village in a separate temporary data frame
tmp2 <- as.data.frame(table(tmp$ov16rdt)) #Summarizing the OV16 variable separated based on the temporary village data frame
print(i) #printing the village used for the calculation
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100) #Printing the calculated prevalence percentage for that particular village
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"]))) #Performing a proportion test on the OV16 positive present in the selected village and the total number of persons tested within that village
print(proptest) #Printing the proportion test
}

#"Gabat":     4% (0.7-14.9%)
#"KazanaI":   44% (30.3-58.7%)
#"KazanaII":  50% (36.1-63.9%)
#"Matara":    34.6% (22.3-49.2%)
#"Tarawa":    22% (12.0-36.3%)

#####Prevalence per age OV16
for (i in 3:9 ) #creating a separate variable based on the ages testes to use in a looped based prevalence calculation
{tmp <- subset.data.frame(DATA, age == i)  #For loop based on the previously created separate variable
tmp2 <- as.data.frame(table(tmp$ov16rdt))
print(i)
print(prod(tmp2[1:1, "Freq"])) #Printing the amount of OV16 positive cases within the age
print(length(tmp$patient_id[!is.na(tmp$ov16rdt)]))  #Printing the amount of cases tested within the age
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100)
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
}

# 3: 8.823529% <- 3/34 (2.3-24.8%)
# 4: 28.57143% <- 14/49 (16.7-42.7%)
# 5: 36.36364% <- 12/33 (21.0-54.9%)
# 6: 50%       <- 18/36 (34.5-65.5%)
# 7: 23.33333% <- 7/30 (10.6-42.7%)
# 8: 30%       <- 9/30 (17.3-52.5%)
# 9: 35.29412% <- 12/34 (20.3-53.5%)

#####Prevalence per age group OV16
##Age group 3-6
#OVERALL
tmp <- subset.data.frame(DATA, age<=6) #Creating a separate data frame based on age group below or equal to 6
table(tmp$ov16rdt)
tmp2 <- as.data.frame(table(tmp$ov16rdt))
print(length(tmp$patient_id))
print(length(tmp$patient_id[!is.na(tmp$ov16rdt)]))
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100)
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
#ALL: 30.7% <- 47/153 (23.7-38.8%)

#PER VILLAGE
villages <- c(levels(DATA$village))
for (i in villages )
{tmp <- subset.data.frame(DATA, village == i)
tmp1 <- subset.data.frame(tmp, age<=6) #Creating a separate data frame based on age group below or equal to 6 of a separate village
tmp2 <- as.data.frame(table(tmp1$ov16rdt))
print(i)
print(prod(tmp2[1:1, "Freq"]))
print(length(tmp1$patient_id[!is.na(tmp1$ov16rdt)]))
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100)
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
}

#"Gabat":     3.448276%   -> 1/29 (0.2-19.6%)
#"KazanaI":   46.66667%   -> 14/30 (28.8-65.4%)
#"KazanaII":  43.75%      -> 14/32 (26.0-60.6%)
#"Matara":    35.48387%   -> 11/31 (19.8-54.6%)
#"Tarawa":    23.33333%   -> 7/30 (10.6-42.7%)

##Age group 7-9
tmp <- subset.data.frame(DATA, age >6) #Creating a separate data frame based on age group above 6
table(tmp$ov16rdt)
tmp2 <- as.data.frame(table(tmp$ov16rdt))
print(length(tmp$patient_id))
print(length(tmp$patient_id[!is.na(tmp$ov16rdt)]))
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100)
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
#ALL: 29.78723% <- 22/95 (21.7-40.9%)


#PER VILLAGE
villages <- c(levels(DATA$village))
for (i in villages )
{tmp <- subset.data.frame(DATA, village == i)
tmp1 <- subset.data.frame(tmp, age>6) #Creating a separate data frame based on age group above 6 of a separate village
tmp2 <- as.data.frame(table(tmp1$ov16rdt))
table(tmp1$ov16rdt)
print(i)
print(prod(tmp2[1:1, "Freq"]))
print(length(tmp1$patient_id[!is.na(tmp1$ov16rdt)]))
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100)
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
}

#"Gabat":     4.761905%   -> 1/21 (0.2-25.9%)
#"KazanaI":   40%         -> 8/20 (20.0-63.6%)
#"KazanaII":  66.66667%   -> 9/13 (38.9-89.6)
#"Matara":    33.33333%   -> 7/21 (15.5-57.0%)
#"Tarawa":    20%         -> 4/20 (0.7-44.3%)

##Statistics: Is there a difference in seroprevalence between the age groups?
prop.test(x = c(47,22), n = c(153,95)) #(OV16+ age below/equal to 6,OV16+ age above 6)(tested age below/equal to 6, tested age above 6 )
#P-value = 0.2518


#####Prevalence per sex OV16
##FEMALE
#OVERALL
tmp <- subset.data.frame(DATA, sex == 1)
table(tmp$ov16rdt)
tmp2 <- as.data.frame(table(tmp$ov16rdt))
print(prod(tmp2[1:1, "Freq"]))
print(length(tmp$patient_id[!is.na(tmp$ov16rdt)]))
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100)
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
#ALL: 31.85185% <- 43/135 (24.3-40.5%)

##MALE
#OVERALL
tmp <- subset.data.frame(DATA, sex == 2)
tmp2 <- as.data.frame(table(tmp$ov16rdt))
table(tmp$ov16rdt)
print(prod(tmp2[1:1, "Freq"]))
print(length(tmp$patient_id[!is.na(tmp$ov16rdt)]))
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100)
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
#ALL: 29.20354% <- 33/113 (21.2-38.6%)

##Statistics:Is there a difference in the overall seroprevalence between sexes?
prop.test(x = c(43,33), n = c(135,113))
#P-value = 0.7548

#PER VILLAGE
village <- c(levels(DATA$village))
sex <- c(levels(DATA$sex))
for (i in village) #For loop based on two variables, here village (i) and sex (j)
{for (j in sex)
{tmp <- subset.data.frame(DATA, village == i)
tmp1 <- subset.data.frame(tmp, sex == j)
tmp2 <- as.data.frame(table(tmp1$ov16rdt))
print(i)
print (j)
print(prod(tmp2[1:1, "Freq"]))
print((prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"])))
print(((prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"])))*100))
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
}}

##Female;
#"Gabat":     4.545455%   -> 1/22 (0.2-24.9%)
#"KazanaI":   41.17647%   -> 14/34 (25.2-59.2%)
#"KazanaII":  46.15385%   -> 12/26 (27.1-66.3%)
#"Matara":    32.35294%   -> 11/34 (18.0-50.6%)
#"Tarawa":    26.31579%   -> 5/19 (10.1-51.4%)

##Male:
#"Gabat":     3.571429%   -> 1/28 (0.2-20.2%)
#"KazanaI":   50%         -> 8/16 (28.0-72.0%)
#"KazanaII":  55%         -> 11/20 (32.0-76.2%)
#"Matara":    38.88889%   -> 7/18 (18.3-64.0%)
#"Tarawa":    19.35484%   -> 6/31 (8.1-38.1%)

#####Prevalence dermatitis - OV16
table(DATA$dermititis[!is.na(DATA$ov16rdt)],DATA$ov16rdt[!is.na(DATA$ov16rdt)] )
#    1   2
#0  35 120
#1  40  52

##prevalence:
#Dermatitis:
40/(40+52) #43.5%
#No dermatitis:
35/(35+120) #22.6%

##Statistics:Is there a difference in seroprevalence between having dermatitis or not?
prop.test(x = c(35,40), n = c(120+35,40+52)) #p-value = 0.000933


#####Prevalence epilepsy
table(DATA$epilepsy[!is.na(DATA$ov16rdt)],DATA$ov16rdt[!is.na(DATA$ov16rdt)] )
#    1   2
#0  76 168
#1   0   4

##Prevalence:
#Without epilepsy
76/(76+168) #31.1%

##Statistics:Is there a difference in seroprevalence between having epilepsy or not?
prop.test(x = c(0,76), n = c(4,76+168)) #--> STATISTICALLY WRONG


#####Prevalence Ivermectin
table(DATA$ivermectin[!is.na(DATA$ov16rdt)],DATA$ov16rdt[!is.na(DATA$ov16rdt)] )
#    1   2
#0  48 124
#1  28  48

##Prevalence:
#not taken ivermectin
48/(76+168) #19.7%
#taken ivermectin
28/(28+46) #37.8%

##Statistics:Is there a difference in seroprevalence between having taken ivermectin or not?
prop.test(x = c(48,28), n = c(48+124,28+48)) #p-value: 0.209

#####################
## DATA VISUALISATION
#####################

##Data upload
peragepervillage<- read_excel("Seroprevalence_dataset.xlsx", sheet = "Prevalence_Age_Village") #Table including all previously calculated prevalences
View(peragepervillage)

peragepervillage$age <- as.factor(peragepervillage$age) #Changing the age to a factorial variable to prepare for easy ggplot creation

p <- ggplot(data=peragepervillage, aes(x=village, y=prevalence_2023, fill=age)) + #ggplot based in the previously uploaded data, villages on X axis, prevalence on Y axis and split by age factors
  geom_bar(stat="identity", color="black", position=position_dodge())+ #Choosing for a bar-plot with black identity and making sure no overlap is present
  theme_minimal()+ #Choosing the theme of the barplot
  labs(x="Village",y="Prevalence",  #Changing the titals of the axes
       title="OV16 Seroprevalence per village and per age") #Giving the graph an appropriate title

p <- p + scale_fill_grey() + # Change the colors manually to grey scale
  theme_bw() + #Choosing a black and white theme for the general plot
  scale_x_discrete(labels = c('Hai Gabat','Kazana 1','Kazana 2', 'Hai Matara', 'Hai Tarawa')) #Changing names to those appropriate for the paper

print(p) #Preview

##Saving plot as pdf within your directory
pdf("Seroprevalence_village_age.pdf") 
print(p)
dev.off()

##Saving plot as EPS within your directory
setEPS()
postscript("Seroprevalence_village_age.eps")
print(p)
dev.off()

###########################
###########################
## DATA 2019 COMPARISON ###
###########################
###########################
###############
## READ IN DATA
###############
MERGED_DATA <- read_excel("Seroprevalence_dataset.xlsx",  sheet = "2019_2023")
View(MERGED_DATA)

##Creating an extra variable to divide the villages in high and low transmission zones
MERGED_DATA$zone <- ifelse(MERGED_DATA$village == "MATARA", "High", ifelse(MERGED_DATA$village == "KAZANA2","High", ifelse(MERGED_DATA$village == "KAZANA1","High","Low")))
table(MERGED_DATA$zone)

##Changing variables to factorial factors
names <- c('zone','survey','epilepsy', 'village')
MERGED_DATA[,names] <- lapply(MERGED_DATA[,names] , factor)

##################
## OV16 PREVALENCE
##################
#####Prevalence overall between surveys OV16
table(MERGED_DATA$ov16rdt_ORIGINAL,MERGED_DATA$survey)

#   2019 2023
#1   35   76
#2  109  172

##Statistics: is there a difference in OV16 seroprevalence between the surveys?
prop.test(x = c(35,76), n = c(35+109,76+172)) #p-value = 0.2199

#####Prevalence per transmission zone OV16
table(MERGED_DATA$zone,MERGED_DATA$survey,MERGED_DATA$ov16rdt_ORIGINAL)

zone <- c(levels(MERGED_DATA$zone))
survey <- c(levels(MERGED_DATA$survey))
for (i in zone) #For loop based on two variables, here zone (i) and survey (j)
{for (j in survey)
{tmp <- subset.data.frame(MERGED_DATA, zone == i)
tmp1 <- subset.data.frame(tmp, survey == j)
tmp2 <- as.data.frame(table(tmp1$ov16rdt_ORIGINAL))
print(i)
print (j)
print(prod(tmp2[1:1, "Freq"]))
print((prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"])))
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100)
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
}}

#High/2019: 32.96703%   -> 30/91 (23.7-43.7%)
#High/2023: 42.56757%   -> 63/148 (34.6-51.0%)
#Low/2019: 9.433962%   -> 5/53 (3.5-21.4%)
#Low/2023: 13%  -> 13/100 (7.4-21.6%)

##Statistics:Is there a difference in seroprevalence between HTZs?
prop.test(x = c(30,63), n = c(30+61,63+85)) 
#P-value: 0.180

##Statistics:Is there a difference in seroprevalence between LTZs?
prop.test(x = c(5,13), n = c(5+48,13+87))
#P-value: 0.698


#####Prevalence per village OV16
table(MERGED_DATA$ov16rdt_ORIGINAL, MERGED_DATA$village, MERGED_DATA$survey)
levels(MERGED_DATA$village)
village <- c("KAZANA1","KAZANA2","MATARA","TARAWA") #NO Gabath because of too low sample size
survey <- c(levels(MERGED_DATA$survey))
for (i in village)
{for (j in survey)
{tmp <- subset.data.frame(MERGED_DATA, village == i)
tmp1 <- subset.data.frame(tmp, survey == j)
tmp2 <- as.data.frame(table(tmp1$ov16rdt_ORIGINAL))
print(i)
print (j)
print(prod(tmp2[1:1, "Freq"]))
print((prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"])))
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100)
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
}}

#2019:
# KAZANA1   -> 14/40 (21.1-51.7%)
# KAZANA2   -> 12/24 (31.4-68.6%)
# MATARA    -> 4/27 (4.9-34.6%)
# TARAWA    -> 5/29 (6.5-36.5%)


#2023:
# KAZANA1   -> 22/50 (30.3-58.7%)
# KAZANA2   -> 23/46 (36.1-63.9%)
# MATARA    -> 18/52 (22.3-49.2%)
# TARAWA    -> 11/50 (12.0-36.3%)

##Statistics:Is there a difference in seroprevalence between individual villages?
#Gabat: Imposible to test
#Kazana1:
prop.test(x=c(14,22),n=c(40,50)) #p-value = 0.516
# KAZANA2:
prop.test(x=c(12,23),n=c(24,46)) #p-value = 1
# MATARA:
prop.test(x=c(4,18),n=c(27,52)) #p-value = 0.1101
# TARAWA:   
prop.test(x=c(5,11),n=c(29,50)) #p-value = 0.8283

#####Prevalence epilepsy
table(MERGED_DATA$epilepsy,MERGED_DATA$survey,MERGED_DATA$ov16rdt_ORIGINAL)

#OV16 positive
#   2019 2023
#0   33   76
#1    2    0

#OV16 negative
#   2019 2023
#0  102  168
#1    7    4

#####Prevalence per age OV16
table(MERGED_DATA$age,MERGED_DATA$survey,MERGED_DATA$ov16rdt_ORIGINAL)

survey <- c(levels(MERGED_DATA$survey))
for (i in 3:9)
{for (j in survey)
{tmp <- subset.data.frame(MERGED_DATA, age == i)
tmp1 <- subset.data.frame(tmp, survey == j)
tmp2 <- as.data.frame(table(tmp1$ov16rdt_ORIGINAL))
print(i)
print (j)
print(prod(tmp2[1:1, "Freq"]))
print((prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"])))
print(prod(tmp2[1:1, "Freq"])/(prod(tmp2[1:1, "Freq"])+ prod(tmp2[2:2, "Freq"]))*100)
proptest <- prop.test(prod(tmp2[1:1, "Freq"]),(prod(tmp2[1:1, "Freq"])+prod(tmp2[2:2, "Freq"])))
print(proptest)
}}

#2019:
# 3   -> 3/24 (3.3-33.5%)
# 4   -> 5/31 (6.1-34.5%)
# 5   -> 5/17 (11.4-56.0%)
# 6   -> 6/24 (10.6-47.1%)
# 7   -> 11/30 (20.5-56.1%)
# 8   -> 2/10 (3.5-55.8%)
# 9   -> 3/7 (11.8-79.8%)

#2023:
# 3   -> 3/34 (2.3-24.8%)
# 4   -> 14/50 (16.7-42.7%)
# 5   -> 12/33 (21.0-54.9%)
# 6   -> 18/36 (34.5-65.5%)
# 7   -> 7/30 (10.6-42.7%)
# 8   -> 10/31 (17.3-51.5%)
# 9   -> 12/34 (20.3-53.5%)

##Statistics: Does the seroprevalence differ between ages and surveys?
prop.test(x = c(3,14,12,18,7,10,12), n = c(3+31,14+36,12+21,18+18,7+23,10+21,12+22))
#P-value: 0.015

ignore_fit <- lm(ov16rdt_ORIGINAL~ age, data = MERGED_DATA[MERGED_DATA$survey == "2023",])
summary(ignore_fit)
#P-value: 0.068


#####################
## DATA VISUALISATION
#####################
prevalence_2023_2019 <- read_excel("Seroprevalence_dataset.xlsx",  sheet = "Prevalence_Age_Survey", col_types = c("numeric", "numeric", "numeric"))
View(prevalence_2023_2019)

prevalence_2023_2019$AGE <- as.factor(prevalence_2023_2019$AGE)
prevalence_2023_2019$YEAR <- as.factor(prevalence_2023_2019$YEAR)
prevalence_2023_2019$PREVALENCE <- prevalence_2023_2019$PREVALENCE*100

p_2 <- ggplot(data=prevalence_2023_2019, aes(x=AGE, y=PREVALENCE, fill=YEAR)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+
  labs(x="Age",y="Prevalence", 
       title="OV16 Seroprevalence per age per survey")+
  geom_text(aes(label = round(PREVALENCE, digits = 2)),  colour = "black", size = 3,
            vjust = 1.5, position = position_dodge(.9))

p_2 <- p_2 + scale_fill_grey(start=0.8,end =0.5) +
  theme_bw()

print(p_2)

#Printing as PDF
pdf("Seroprevalence_survey_age.pdf")
print(p_2)
dev.off()

#Printing as EPS
setEPS()
postscript("Seroprevalence_survey_age.eps")
print(p_2)
dev.off()