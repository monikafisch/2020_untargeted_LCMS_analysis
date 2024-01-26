#######################################################
########### Untargeted Metabolomic Analysis ########### 
#######################################################
### 1. Filter/Wrangle MS-DIAL output
### 2. Explore data with Venn Diagrams
#######################################################
#######################################################
#
# Peak Height data exported from MS-DIAL as a .txt
# .txt file edited in Excel: 
### Remove first three rows of junk
### Change the name of column1 to "AlignmentID" 
### Most Importantly: change all peak heigh values to scientific notation with 6 decimal places
### save-as a .csv file and load here:
library(data.table)
msdat <- fread(file.choose())

#check the class of each column:
msdat[ ,lapply(.SD, class)] #all peak height data columns should be numeric!

## Basic anatomy of the MS-DIAL output:
## rows = features or ions
## columns ~1-29 = info about each feature
## columns ~30-n = peak height values for each of your samples

# create a reference table that lists column names and column numbers
# because we will be refering to columns by there number
msdat.colnames <- data.table(colnames(msdat)) 

##################################
######## FILTER RAW DATA #########
##################################
# keep only the features that have a Peak Height value >100000 in any sample (designate the sample columns that you want!)
# In this example Peak Height values are in columns #32-47
msdat.sub1 <- msdat[apply(msdat[, c(32:47)], MARGIN = 1, function(x) any(x > 100000)), ]
# Note: MS-DIAL also filters your data based on a threshold, i.e. 100000, so this step may be redundant

# keep only the features that have a Peak Height value <100000 in all Blank samples
# In this example, Blanks are in columns #35-37
msdat.sub2 <- msdat.sub1[apply(msdat.sub1[, c(35:37)], MARGIN = 1, function(x) all(x < 100000)), ]

# remove features that dont have an MS2
# Note: depending on the version of MS-DIAL, this column name may vary!
msdat.sub3 <- msdat.sub2[!(msdat.sub2$`MS/MS assigned`=="FALSE"), ]

# Export:
write.csv(msdat.sub3, "your-file-path-here")


# look up info about specific features!
msdat.sub3[AlignmentID == "4257", 'Spectrum reference file name']
# "Where the AlignmentID is 4257, tell me what is on the same row in the 'Spectrum reference file name' column"

################################################
##### Remove features in Extract Controls ######
## and build presence/absence tables for Venn ##
################################################

# Create new datatables for each sample
# "grab columns that contain "A_"...", etc.
msdat_A <- msdat.sub3[ ,grepl("A_", names(msdat.sub3)), with=FALSE]
msdat_B <- msdat.sub3[ ,grepl("B_", names(msdat.sub3)), with=FALSE]

# Add AlignmentID column to each datatable
msdat_A[ ,AlignmentID:=msdat.sub3$AlignmentID]
msdat_B[ ,AlignmentID:=msdat.sub3$AlignmentID]

# Remove all features that are present in EtAc Extract Controls 
msdat_A.sub1 <- msdat_A[apply(msdat_A[ ,grepl("ExtractControl", names(msdat_A)), with=FALSE], MARGIN = 1, function(x) all(x < 100000)), ] 
msdat_B.sub1 <- msdat_B[apply(msdat_B[ ,grepl("ExtractControl", names(msdat_B)), with=FALSE], MARGIN = 1, function(x) all(x < 100000)), ] 

# Use melt() to re-arrange the table to make it easier to build presence/absence table
msdat_A.sub1.melt <- melt(msdat_A.sub1, id.vars="AlignmentID")
msdat_B.sub1.melt <- melt(msdat_B.sub1, id.vars="AlignmentID")

# Create sample column with unique sample names for each set of there reps
msdat_A.sub1.melt[grepl("SoilSample", variable), sample:="SampleA"]
msdat_A.sub1.melt[grepl("ExtractControl", variable), sample:="ExtractControlA"]
msdat_B.sub1.melt[grepl("SoilSample", variable), sample:="SampleB"]
msdat_B.sub1.melt[grepl("ExtractControl", variable), sample:="ExtractControlB"]

# Check for NAs, a.k.a. samples that didn't get a name in the "sample" column
# "Empty data.table" means there are no NA's! Good Job!
# If there are NA's present, go back to the previous step and give those samples a name!
msdat_A.sub1.melt[is.na(sample), ]
msdat_B.sub1.melt[is.na(sample), ]

# Count the number of reps for each sample that have a value >100000:
# Create a new column with these values called "count_samples_gr_100k"
msdat_A.sub1.melt[ ,count_samples_gr_100k:=sum(value>100000), by=list(AlignmentID, sample)]
msdat_B.sub1.melt[ ,count_samples_gr_100k:=sum(value>100000), by=list(AlignmentID, sample)]

# Assign presence/absence: if at least 2 of reps of a sample are >100000 =1, otherwise =0
# Create a new column with these values called "count_samples_gr_100k_pa"
msdat_A.sub1.melt[ ,count_samples_gr_100k_pa:=ifelse(count_samples_gr_100k>=2, 1, 0), by=list(sample, AlignmentID)]
msdat_B.sub1.melt[ ,count_samples_gr_100k_pa:=ifelse(count_samples_gr_100k>=2, 1, 0), by=list(sample, AlignmentID)]

# Use dcast to reverse the effect of melt()
# Creates a table of presence/absence information for each AlignementID (rows) by samples (columns)
VennDat_A <- dcast.data.table(unique(msdat_A.sub1.melt[ ,list(AlignmentID, sample, count_samples_gr_100k_pa)]), 
                                  AlignmentID~sample, value.var = "count_samples_gr_100k_pa")
VennDat_B <- dcast.data.table(unique(msdat_B.sub1.melt[ ,list(AlignmentID, sample, count_samples_gr_100k_pa)]), 
                              AlignmentID~sample, value.var = "count_samples_gr_100k_pa")
        
# Merge the above tables together, by matching the AlignmentIDs in each table
# Order here matters! Always start with the largest table, then merge smaller tables onto it.
merge1 <- merge(msdat.sub3, VennDat_A, by="AlignmentID", all.x=TRUE)
MSDIAL_plus_VennDat_Example <- merge(merge1, VennDat_B, by="AlignmentID", all.x=TRUE)

colnames(MSDIAL_plus_VennDat_Example)

# Now you have a filtered table with all the info about every feature and
# weather or not it was determined to be "present" or "absent" in each sample!
# Export: 
write.csv(MSDIAL_plus_VennDat_Example, "your-file-path-here")


######################################
#### BUILD VENN DIAGRAMS MANUALLY ####
######################################
#
# NOTE: This is particularly handy when you have MANY samples (i.e. greater than ~3-5)
# With many samples, a proper Venn Diagram can be overwhelming to digest, 
# but with the method below, you can extract exactly what you want, for example:
# the features present in Samples #1-5, but non present in Samples #7, 10, and 15
#
# If you really want to make a complicated 4+ sample Venn Diagram,
# Use this manual method to generate the numbers and draw the Venn Diagram by hand
#
# Mannually generate the numbers for each section of a Venn Diagram
# Change 1's to 0's to get all possible combinations
# For exmple, this will give you the number of features present ("1")
# in both SampleA AND SampleB:
nrow(MSDIAL_plus_VennDat_Example[SampleA == 1 & 
                                 SampleB == 1, ])

# Export a table with only the features of interest 
# (i.e. features present in both samples)
temp <- MSDIAL_plus_VennDat_Example[SampleA == 1 & 
                                    SampleB == 1]
write.csv(temp, "your-file-path-here")


#########################################################
#### USING THE EULERR PACKAGE TO DRAW A VENN DIAGRAM ####
#########################################################
# NOT recommended if you want to compare more than 3 samples

library(eulerr)

# Create a new datatable with the samples you want to compare in a Venn Diagram
DataForVenn <- data.table(cbind(MSDIAL_plus_VennDat_Example$SampleA, 
                                MSDIAL_plus_VennDat_Example$SampleB))

# Remove NAs, if there are any
DataForVenn <- na.omit(DataForVenn)

# Use the eulerr package to draw the Venn Diagram!
TheVenn <- euler(DataForVenn)
plot(TheVenn, quantities = TRUE, labels = c("Lovely Sample A", "Gross Sample B"))

## Circles are scaled to the total number of features in each sample
## "quantities = TRUE" adds the values in each circle 
## (these value can be double checked using the manual method above!)
# Export as a .pdf to manipulate colors, fonts, etc. in Adobe Illustrator! 

# There are a bunch of options in the eulerr package..
# Read all about it here!
?euler

# Here's one example of some more options...
plot(TheVenn, quantities = TRUE, 
     lty = "solid", 
     fill = c("red", "blue"), 
     cex = 2,
     fontface = 2, 
     labels = c("Lovely Sample A", "Gross Sample B"))
