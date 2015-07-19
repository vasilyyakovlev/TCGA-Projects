# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)

# Get list of cancer studies at server
getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study: 
## "http://cancergenome.nih.gov/">The Cancer Genome Atlas (TCGA)</a> 
## Breast Invasive Carcinoma project. 825 cases.<br><i>Nature 2012.</i> 
## <a href="https://tcga-data.nci.nih.gov/docs/publications/brca_2012/">
## Raw data via the TCGA Data Portal</a>.
mycancerstudy = getCancerStudies(mycgds)[14, 1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1, 1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[ 9, 1]

# Get data slices for a specified list of genes, genetic profile and case list
y <- getProfileData(mycgds,c('BRCA1', 'NOS1', 'NOS2', 'NOS3'),mygeneticprofile,mycaselist)
# Removing the NA data
y2 <- na.omit(y)

# Low NOSs expression group
ysort <- y2[ order(-y2[,4]), ]
ysort2 <- ysort[1:45, ]
ysort3 <- ysort2[ order(-ysort2[,3]), ]
ysort4 <- ysort3[1:22, ]
ysort5 <- ysort4[ order(-ysort4[,2]), ]
ylow <- ysort5[1:8, ]
# High NOSs Expression group
ysort6 <- y2[ order(y2[,4]), ]
ysort7 <- ysort6[1:60, ]
ysort8 <- ysort7[ order(ysort7[,3]), ]
ysort9 <- ysort8[1:40, ]
ysort10 <- ysort9[ order(ysort9[,2]), ]
ysort11 <- ysort10[1:10, ]
yhigh <- ysort5[1:8, ]
ylow <- ysort10[1:10, ]
z <- getProfileData(mycgds,c('BRCA1', 'NOX3', 'NOX4'),mygeneticprofile,mycaselist)
z2 <- na.omit(z)
#High NOXs expression Group
zsort <- z2[ order(-z2[,3]), ]
zsort2 <- zsort[1:40, ]
zsort3 <- zsort2[ order(-zsort2[,2]), ]
zhigh <- zsort3[1:19, ]
#Low NOXs expression Group
zsort4 <- z2[ order(z2[,3]), ]
zsort5 <- zsort4[1:50, ]
zsort6 <- zsort5[ order(zsort5[,2]), ]
zlow <- zsort6[1:10, ]

#Generation of Plots
boxplot(ylow, main="Low NOSs Expression Group", ylim=c(-3, 3), ylab="Relative mRNA Expression(normali", col=(c("blue","red", "red", "red")))
boxplot(yhigh, main="High NOSs Expression Group", ylim=c(-3, 3), ylab="Relative mRNA Expression", col=(c("blue","red", "red", "red")))
boxplot(zlow, main="Low NOXs Expression Group", ylim=c(-3, 3), ylab="Relative mRNA Expression", col=(c("blue","gold", "gold")))
boxplot(zhigh, main="High NOXs Expression Group", ylim=c(-3, 3), ylab="Relative mRNA Expression", col=(c("blue","gold", "gold")))

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)

