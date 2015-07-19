# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)

# Get list of cancer studies at server
getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[22, 1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1, 1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[10 , 1]
mygeneticprofile
# Get data slices for a specified list of genes, genetic profile and case list
y <- getProfileData(mycgds,c('BRCA1', 'TP53INP1', 'IFNGR1', 'SOCS1'),mygeneticprofile,mycaselist)
y2 <- na.omit(y)
ysort <- y2[ order(y2[,1]), ]
ylow <- ysort[1:12, ]
yhigh <- ysort[185:195, ]
boxplot(ylow)
boxplot(yhigh)
# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)

