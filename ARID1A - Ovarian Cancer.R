# Library activation
library("cgdsr")
library("ggfortify")
library("survival")
library("eeptools")
library("data.table")
library("rms")

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)

# Get list of cancer studies at server
getCancerStudies(mycgds)

# The Cancer Genome Atlas (TCGA)</a> 
# Serous Ovarian Cancer project. 489 cases.<br> 
# <a href="https://tcga-data.nci.nih.gov/docs/publications/ov_2011/">Raw data via the TCGA Data Portal</a>.

mycancerstudy = getCancerStudies(mycgds)[110, 1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[ 1, 1]

#Get Clinical Data for caselist
myclinicalprofile = getClinicalData(mycgds, mycaselist)

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[ 9, 1]

###########################################

#List of genes
a <- getProfileData(mycgds,c('ARID1A'),mygeneticprofile,mycaselist)
#Removing the NAs
a <- na.omit(a)


### Merge with clinical data

armerge <- merge(a, myclinicalprofile, by="row.names")
armerge <- na.omit(armerge)


#High ARID1A expression Group
arsort <- armerge[order(-armerge$ARID1A), ]
arhigh <- arsort[1:35, ]

#Low ARID1A expression Group
arsort <- armerge[order(armerge$ARID1A), ]
arlow <- arsort[1:35, ]


############################################
## Add survival and ExprStatus objects

arhigh$ExprStatus <- with(arhigh, 'High')
arlow$ExprStatus <- with(arlow, 'Low')
merge.ar <- rbind(arhigh, arlow)
merge.ar$SurvObj <- with(merge.ar, Surv(OS_MONTHS, OS_STATUS == 'DECEASED'))


############################################
## Calculation of the p-values

diffar = survdiff(Surv(time=OS_MONTHS, OS_STATUS == 'DECEASED') ~ ExprStatus, data=merge.ar)
pval.ar <- pchisq(diffar$chisq, length(diffar$n)-1, lower.tail = FALSE)
pval.ar


############################################
## Kaplan-Meier survival plot

ar <- survfit(Surv(time=OS_MONTHS, OS_STATUS == 'DECEASED') ~ ExprStatus, data=merge.ar)
autoplot(ar, CI=TRUE, pval=TRUE, plotTable=TRUE, xlab = "Months", ylab = "Survival", legendLabs = NULL, divideTime = 1, returnTable = TRUE)

