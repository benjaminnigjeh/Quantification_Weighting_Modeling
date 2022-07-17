# input the CSV file generated from Skyline containig Replicate Name,	Protein,	
# Peptide,	Transition,	Retention Time,	Area,	and Background columns

rm(list=ls())
setwd(dir = "D:/Quan10273")
Input_DF <- read.csv(file.choose(), header = T)

# make sure that the area background and retention time values are numerical

Input_DF$Area <- as.character(Input_DF$Area)
Input_DF$Area <- as.numeric(Input_DF$Area)

Input_DF$Retention.Time <-as.character(Input_DF$Retention.Time)
Input_DF$Retention.Time <- as.numeric(Input_DF$Retention.Time)

Input_DF$Background <- as.character(Input_DF$Background)
Input_DF$Background <- as.numeric(Input_DF$Background)


# remove the unspecified targets from the data
Subset_DF <- subset(Input_DF, Replicate.Name !="")


# replace the missing data from the specified targets with the value of "1"
Subset_DF$Retention.Time[is.na(Subset_DF$Retention.Time)] <- 1
Subset_DF$Area[is.na(Subset_DF$Area)] <- 1
Subset_DF$Background[is.na(Subset_DF$Background)] <- 1

#acquire replicate names 
Subset_DF$Replicate.Name <- as.character(Subset_DF$Replicate.Name)

#provide the AUC information for each replicate 
Unstack_DF <- unstack(Subset_DF, form = Area ~ Replicate.Name)

#provide the retention time for each replicate
Unstack_Retention_Time_DF <- unstack(Subset_DF, form = Retention.Time ~ Replicate.Name)


#calculate the number of runs in the sample
n <- length(unique(Subset_DF$Replicate.Name))

# assign levels for each run and it should be numeric
Subset_DF$Replicate.Name <- c(1:n)
Subset_DF$Replicate.Name <- as.numeric(Subset_DF$Replicate.Name)

# generate an ID dataframe to extract unique peptide protein and transitions
IDs_DF <- Subset_DF[Subset_DF$Replicate.Name == "1",]

# add the names of transitions proteins and peptides to the dataframe
Unstack_DF$Transition <- IDs_DF$Transition
Unstack_DF$Protein <- IDs_DF$Protein
Unstack_DF$Peptide <- IDs_DF$Peptide

# add the names of transitions proteins and peptides to the retention time dataframe
Unstack_Retention_Time_DF$Transition <- IDs_DF$Transition
Unstack_Retention_Time_DF$Protein <- IDs_DF$Protein
Unstack_Retention_Time_DF$Peptide <- IDs_DF$Peptide

# estimate the variation level for each transtion based on the control runs
Unstack_DF$mean <- apply(Unstack_DF[, c(1,2,n-1,n)], 1, mean)
Unstack_DF$STDev <- apply(Unstack_DF[, c( 1, 2, n-1, n)], 1, sd) 
Unstack_DF$variation <- ((Unstack_DF$STDev)/(Unstack_DF$mean))*100

# estimate the retention variation level for each transtion based on the control runs
Unstack_Retention_Time_DF$mean <- apply(Unstack_Retention_Time_DF[, c(1,2,n-1,n)], 1, mean)
Unstack_Retention_Time_DF$STDev <- apply(Unstack_Retention_Time_DF[, c( 1, 2, n-1, n)], 1, sd) 
Unstack_Retention_Time_DF$variation <- ((Unstack_Retention_Time_DF$STDev)/(Unstack_Retention_Time_DF$mean))*100


#First phase of data analysis: how the quality of our data is

# Plot the retention time variation versus retention time variation
pdf("Retention_Time_Variation.pdf")
plot(Unstack_Retention_Time_DF$STDev ~ Unstack_Retention_Time_DF$mean, xlab = "Retention time / minutes", ylab = "StDev / minutes", ylim=c(0, 1))
dev.off()

# Plot the signal variation versus retention time variation
pdf("AUC_Variation.pdf")
Unstack_DF$Retention.Time <- Unstack_Retention_Time_DF$mean
plot(Unstack_DF$variation ~ Unstack_DF$Retention.Time, xlab = "retntion time / minutes", ylab = "variations", ylim=c(0,50))
dev.off()

# plot TICs for control runs
pdf("CTL_TIC_Variation.pdf")
n <- as.numeric(n)
TICs <- c("CTL1" = sum(Unstack_DF[, 1]),"CTL2" =  sum(Unstack_DF[, 2]), "CTL3" =  sum(Unstack_DF[, n-1]), "CTL4" =  sum(Unstack_DF[, n]))
barplot(TICs, ylim = c(0, 1e10))
dev.off()


# Plot the linear regression of control runs without any normalization


pdf("CTL_Linear_Regression.pdf")
par(mfrow=c(3,1))
plot(Unstack_DF[, 1] ~ Unstack_DF[, 2], xlim = c(0, 8e+08), ylim = c(0, 8e+08))
abline(lm(Unstack_DF[, 1] ~ Unstack_DF[, 2]), col = "blue")
plot(Unstack_DF[, 1] ~ Unstack_DF[, n-1], xlim = c(0, 8e+08), ylim = c(0, 8e+08))
abline(lm(Unstack_DF[, 1] ~ Unstack_DF[, n-1]), col = "blue")
plot(Unstack_DF[, 1] ~ Unstack_DF[, n], xlim = c(0, 8e+08), ylim = c(0, 8e+08))
abline(lm(Unstack_DF[, 1] ~ Unstack_DF[, n]), col = "blue")
dev.off()



# begin with plotting the jpeg file for the vriation distribution
pdf("Model_Distribution.pdf")


# plot the distribution results from the raw data 
test1 <- density(Unstack_DF$variation)
plot(test1, col = "black", main = "", xlab = "Variation percentage", xlim=c(-10, 100), ylim=c(0.00, 0.06))

# add a vertical line at 20% variation as reference
abline(v = 20, col= "blue")



# generate a new copy for the TIC normalization model

Unstack_DF_TIC <- Unstack_DF


# TIC normalize the entire raw files
i <- 1
while (i <= n) {
  Unstack_DF_TIC <- transform(Unstack_DF_TIC, temp = Unstack_DF_TIC[, i]/sum(Unstack_DF_TIC[, i]))
  Unstack_DF_TIC[, i] <- Unstack_DF_TIC$temp
  i <- i + 1
}


# estimation of the variations after TIC normalization
Unstack_DF_TIC$mean <- apply(Unstack_DF_TIC[, c(1,2,n-1,n)], 1, mean)
Unstack_DF_TIC$STDev <- apply(Unstack_DF_TIC[, c( 1, 2, n-1, n)], 1, sd) 
Unstack_DF_TIC$variation <- ((Unstack_DF_TIC$STDev)/(Unstack_DF_TIC$mean))*100

#plot the result from TIC normalization
test2 <- density(Unstack_DF_TIC$variation)
lines(test2, col = "red")




# for background reduction we need to begin from the subsetted raw data
Subset_BackGround_Reduced_DF <- Subset_DF

#substract the background from the area
Subset_BackGround_Reduced_DF$Area <- abs(Subset_BackGround_Reduced_DF$Area - Subset_BackGround_Reduced_DF$Background)

#replace the missing values with 1
Subset_BackGround_Reduced_DF$Retention.Time[is.na(Subset_BackGround_Reduced_DF$Retention.Time)] <- 1
Subset_BackGround_Reduced_DF$Area[is.na(Subset_BackGround_Reduced_DF$Area)] <- 1
Subset_BackGround_Reduced_DF$Background[is.na(Subset_BackGround_Reduced_DF$Background)] <- 1

#unstack the updated area values for each run
Subset_BackGround_Reduced_DF$Replicate.Name <- as.character(Subset_BackGround_Reduced_DF$Replicate.Name)
Unstack_Background_Reduced_DF <- unstack(Subset_BackGround_Reduced_DF, form = Area ~ Replicate.Name)

# rebuilt a dataframe containing transition protein and peptide infor
n <- length(unique(Subset_BackGround_Reduced_DF$Replicate.Name))
Subset_BackGround_Reduced_DF$Replicate.Name <- c(1:n)
Subset_BackGround_Reduced_DF$Replicate.Name <- as.numeric(Subset_BackGround_Reduced_DF$Replicate.Name)
IDs_DF <- Subset_BackGround_Reduced_DF[Subset_BackGround_Reduced_DF$Replicate.Name == "1",]
Unstack_Background_Reduced_DF$Transition <- IDs_DF$Transition
Unstack_Background_Reduced_DF$Protein <- IDs_DF$Protein
Unstack_Background_Reduced_DF$Peptide <- IDs_DF$Peptide

#estimate the variation for the case of background reduction model
Unstack_Background_Reduced_DF$mean <- apply(Unstack_Background_Reduced_DF[, c(1,2,n-1,n)], 1, mean)
Unstack_Background_Reduced_DF$STDev <- apply(Unstack_Background_Reduced_DF[, c( 1, 2, n-1, n)], 1, sd) 
Unstack_Background_Reduced_DF$variation <- ((Unstack_Background_Reduced_DF$STDev)/(Unstack_Background_Reduced_DF$mean))*100
Unstack_Background_Reduced_DF$variation[is.na(Unstack_Background_Reduced_DF$variation)] <- 100

# plot the variation distribution for the backgroun reduction model
test3 <- density(Unstack_Background_Reduced_DF$variation)
lines(test3, col = "green")


# re import the stacked data for the MS2 normalzation
Unstack_MS2_DF <- Unstack_DF

# import the retention time values
Unstack_MS2_DF$Retention.Time <- IDs_DF$Retention.Time

#order the data based on the retention time
Unstack_MS2_DF <- Unstack_MS2_DF[order(Unstack_MS2_DF$Retention.Time),]

# divide control 1 to other 3 controls
Unstack_MS2_DF$temp1 <- Unstack_MS2_DF[, 1]/Unstack_MS2_DF[, 1]
Unstack_MS2_DF$temp2 <- Unstack_MS2_DF[, 1]/Unstack_MS2_DF[, 2]
Unstack_MS2_DF$temp3 <- Unstack_MS2_DF[, 1]/Unstack_MS2_DF[, n-1]
Unstack_MS2_DF$temp4 <- Unstack_MS2_DF[, 1]/Unstack_MS2_DF[, n]


# average the normalized coefficients in a range of 20 values
i<- 1 
while (i <= as.numeric(nrow(Unstack_MS2_DF))) {
  Unstack_MS2_DF[i, "temp5"] <- mean(Unstack_MS2_DF[i-10:i+10,"temp1"])
  i <- i +1
}

i<- 1 

while (i <= as.numeric(nrow(Unstack_MS2_DF))) {
  Unstack_MS2_DF[i, "temp6"] <- mean(Unstack_MS2_DF[i-10:i+10,"temp2"])
  i <- i +1
}

i<- 1 

while (i <= as.numeric(nrow(Unstack_MS2_DF))) {
  Unstack_MS2_DF[i, "temp7"] <- mean(Unstack_MS2_DF[i-10:i+10,"temp3"])
  i <- i +1
}

i<- 1 

while (i <= as.numeric(nrow(Unstack_MS2_DF))) {
  Unstack_MS2_DF[i, "temp8"] <- mean(Unstack_MS2_DF[i-10:i+10,"temp4"])
  i <- i +1
}


# multiply the original raw control values with the coefficents
Unstack_MS2_DF$temp9 <- Unstack_MS2_DF$temp5*Unstack_MS2_DF[, 1]
Unstack_MS2_DF$temp10 <- Unstack_MS2_DF$temp5*Unstack_MS2_DF[, 2]
Unstack_MS2_DF$temp11 <- Unstack_MS2_DF$temp5*Unstack_MS2_DF[, 3]
Unstack_MS2_DF$temp12 <- Unstack_MS2_DF$temp5*Unstack_MS2_DF[, 4]

# estimate the variation after MS2 normalization
Unstack_MS2_DF$mean <- apply(Unstack_MS2_DF[, c("temp9", "temp10", "temp11", "temp12")], 1, mean)
Unstack_MS2_DF$STDev <- apply(Unstack_MS2_DF[, c( "temp9", "temp10", "temp11", "temp12")], 1, sd) 
Unstack_MS2_DF$variation <- ((Unstack_MS2_DF$STDev)/(Unstack_MS2_DF$mean))*100
Unstack_MS2_DF$variation[is.na(Unstack_MS2_DF$variation)] <- 100

# add the variation model to the plot
test4 <- density(Unstack_MS2_DF$variation)
lines(test4, col = "blue")

# a hybrid model based on MS2 normalization and TIC normalization
Unstack_DF_TIC_MS2 <- Unstack_DF_TIC

# import the retention times and sort based on retention time for MS2 normalization
Unstack_DF_TIC_MS2$Retention.Time <- IDs_DF$Retention.Time
Unstack_DF_TIC_MS2 <- Unstack_DF_TIC_MS2[order(Unstack_DF_TIC_MS2$Retention.Time),]


# TIC normalization
Unstack_DF_TIC_MS2$temp1 <- Unstack_DF_TIC_MS2[, 1]/Unstack_DF_TIC_MS2[, 1]
Unstack_DF_TIC_MS2$temp2 <- Unstack_DF_TIC_MS2[, 1]/Unstack_DF_TIC_MS2[, 2]
Unstack_DF_TIC_MS2$temp3 <- Unstack_DF_TIC_MS2[, 1]/Unstack_DF_TIC_MS2[, n-1]
Unstack_DF_TIC_MS2$temp4 <- Unstack_DF_TIC_MS2[, 1]/Unstack_DF_TIC_MS2[, n]

#MS2 normalization coefficients

i<- 1 
while (i <= as.numeric(nrow(Unstack_DF_TIC_MS2))) {
  Unstack_DF_TIC_MS2[i, "temp5"] <- mean(Unstack_DF_TIC_MS2[i-10:i+10,"temp1"])
  i <- i +1
}

i<- 1 
while (i <= as.numeric(nrow(Unstack_DF_TIC_MS2))) {
  Unstack_DF_TIC_MS2[i, "temp6"] <- mean(Unstack_DF_TIC_MS2[i-10:i+10,"temp2"])
  i <- i +1
}

i<- 1 
while (i <= as.numeric(nrow(Unstack_DF_TIC_MS2))) {
  Unstack_DF_TIC_MS2[i, "temp7"] <- mean(Unstack_DF_TIC_MS2[i-10:i+10,"temp3"])
  i <- i +1
}

i<- 1 
while (i <= as.numeric(nrow(Unstack_DF_TIC_MS2))) {
  Unstack_DF_TIC_MS2[i, "temp8"] <- mean(Unstack_DF_TIC_MS2[i-10:i+10,"temp4"])
  i <- i +1
}

#multiplication of the MS2 normalization coefficients
Unstack_DF_TIC_MS2$temp9 <- Unstack_DF_TIC_MS2$temp5*Unstack_DF_TIC_MS2[, 1]
Unstack_DF_TIC_MS2$temp10 <- Unstack_DF_TIC_MS2$temp5*Unstack_DF_TIC_MS2[, 2]
Unstack_DF_TIC_MS2$temp11 <- Unstack_DF_TIC_MS2$temp5*Unstack_DF_TIC_MS2[, 3]
Unstack_DF_TIC_MS2$temp12 <- Unstack_DF_TIC_MS2$temp5*Unstack_DF_TIC_MS2[, 4]

#Estimation of the variations
Unstack_DF_TIC_MS2$mean <- apply(Unstack_DF_TIC_MS2[, c("temp9", "temp10", "temp11", "temp12")], 1, mean)
Unstack_DF_TIC_MS2$STDev <- apply(Unstack_DF_TIC_MS2[, c( "temp9", "temp10", "temp11", "temp12")], 1, sd) 
Unstack_DF_TIC_MS2$variation <- ((Unstack_DF_TIC_MS2$STDev)/(Unstack_DF_TIC_MS2$mean))*100
Unstack_DF_TIC_MS2$variation[is.na(Unstack_DF_TIC_MS2$variation)] <- 100

# adding the line of TIC and MS2 normalization to the plot
test5 <- density(Unstack_DF_TIC_MS2$variation)
lines(test5, col = "orange")

# legend for each model with the same color that comes at the end prior to dev off
legend("topright", legend= c("Raw data", "TIC norm", "Background reduced", "MS2 norm", "MS2 TIC norm"), col= c("black", "red", "green", "blue", "orange"), lty=1:2, cex=0.8)

dev.off()


# report the variations for every fragment based on every model

Fragments_variations_DF <- Unstack_DF[, c("Protein", "Peptide", "Transition")]
Fragments_variations_DF$raw_variation <- round(Unstack_DF$variation, digits = 1)
Fragments_variations_DF$TIC_variation <- round(Unstack_DF_TIC$variation, digits = 1)
Fragments_variations_DF$MS2_variation <- round(Unstack_MS2_DF$variation, digits = 1)
Fragments_variations_DF$TIC_MS2_variation <- round(Unstack_DF_TIC_MS2$variation, digits =1)
Fragments_variations_DF$Background_variation <- round(Unstack_Background_Reduced_DF$variation, digits = 1)

# number of columns, we'll use this later
cols <- length(Fragments_variations_DF[1, ])-3 


# exporting data.frame to excel is easy with xlsx package
require(xlsx)
sheetname <- "mysheet"
write.xlsx(Fragments_variations_DF, "D:/Quan10273/Fragment_Report.xlsx", sheetName=sheetname)
file <- "D:/Quan10273/Fragment_Report.xlsx"

# we want to highlight cells based on their value
# load workbook
wb <- loadWorkbook(file)              

# Red
fo2 <- Fill(foregroundColor="#FF0000")    
a2 <-  Alignment(h="ALIGN_CENTER")
cs2 <- CellStyle(wb, fill=fo2, alignment = a2) 

#Yellow
fo4 <- Fill(foregroundColor="#FFFF00")    
a4 <-  Alignment(h="ALIGN_CENTER")
cs4 <- CellStyle(wb, fill=fo4, alignment = a4) 

#Green
fo5 <- Fill(foregroundColor="#00FF00")
a5 <-  Alignment(h="ALIGN_CENTER")
cs5 <- CellStyle(wb, fill=fo5, alignment = a5)         


# get all sheets
sheets <- getSheets(wb)  

# get specific sheet
sheet <- sheets[[sheetname]]

# get rows: 1st row is "not" headers
rows <- getRows(sheet, rowIndex=1:(nrow(Fragments_variations_DF)+1))     

# get cells: data begins from first column and goes on
cells <- getCells(rows, colIndex = 1:cols+4)          


# extract the cell values
values <- lapply(cells, getCellValue) 


# find cells meeting red conditional criteria
highlightred <- NULL
for (i in names(values)) {
  x <- as.numeric(values[i])
  if (x > 50 && !is.na(x)) {
    highlightred <- c(highlightred, i)
  }    
}


# find cells meeting yellow conditional criteria 
highlightyellow <- NULL
for (i in names(values)) {
  x <- as.numeric(values[i])
  if (50 > x && x > 20 && !is.na(x)) {
    highlightyellow <- c(highlightyellow, i)
  }    
}

# find cells meeting green conditional criteria 
highlightgreen <- NULL
for (i in names(values)) {
  x <- as.numeric(values[i])
  if (20 >= x && !is.na(x)) {
    highlightgreen <- c(highlightgreen, i)
  }    
}

# find cells meeting conditional criteria less than -200%
highlightblue <- NULL
for (i in names(values)) {
  x <- as.numeric(values[i])
  if (x < -200 && !is.na(x)) {
    highlightblue <- c(highlightblue, i)
  }    
}




lapply(names(cells[highlightred]),
       function(ii) setCellStyle(cells[[ii]], cs2))

lapply(names(cells[highlightyellow]),
       function(ii) setCellStyle(cells[[ii]], cs4))

lapply(names(cells[highlightgreen]),
       function(ii) setCellStyle(cells[[ii]], cs5))


saveWorkbook(wb, file)


# to generate a new folder Raw_Data_Fragments and set it as the working directory   



# Unique number of proteins that are represented in the data set 
proteins <- as.character(unique(Unstack_DF$Protein))
proteins <- proteins[!is.na(proteins)]
n <- length(proteins)


# export the fragmentation results from the raw data
dir.create("Raw_Data_Fragments")
setwd("Raw_Data_Fragments")
for (m in 1:n){
  temp1 <-Unstack_DF[ Unstack_DF$Protein == proteins[m], ]
  temp1 <- as.data.frame(temp1)
  write.csv(temp1, file = ((paste(proteins[m], "Raw_Result.CSV", sep = "_"))))
}

# export the fragmentation results from the TIC data
setwd(dir = "D:/Quan10273")
dir.create("D:/Quan10273/TIC_Data_Fragments")
setwd("TIC_Data_Fragments")
for (m in 1:n){
  temp1 <-Unstack_DF_TIC[ Unstack_DF_TIC$Protein == proteins[m], ]
  temp1 <- as.data.frame(temp1)
  write.csv(temp1, file = ((paste(proteins[m], "TIC_Result.CSV", sep = "_"))))
}

# export the fragmentation results from the background reduction data
setwd(dir = "D:/Quan10273")
dir.create("D:/Quan10273/Background_Reduction_Fragments")
setwd("Background_Reduction_Fragments")
for (m in 1:n){
  temp1 <-Unstack_Background_Reduced_DF[ Unstack_Background_Reduced_DF$Protein == proteins[m], ]
  temp1 <- as.data.frame(temp1)
  write.csv(temp1, file = ((paste(proteins[m], "Background_Reduced_Result.CSV", sep = "_"))))
}

# export the fragmentation results from the MS2 normalzied data
setwd(dir = "D:/Quan10273")
dir.create("D:/Quan10273/MS2_Normalized_Fragments")
setwd("MS2_Normalized_Fragments")
for (m in 1:n){
  temp1 <-Unstack_MS2_DF[ Unstack_MS2_DF$Protein == proteins[m], ]
  temp1 <- as.data.frame(temp1)
  write.csv(temp1, file = ((paste(proteins[m], "MS2_Normalzied_Result.CSV", sep = "_"))))
}

# export the fragmentation results from the MS2 and TIC normalzied data
setwd(dir = "D:/Quan10273")
dir.create("D:/Quan10273/MS2_TIC_Normalized_Fragments")
setwd("MS2_TIC_Normalized_Fragments")
for (m in 1:n){
  temp1 <-Unstack_DF_TIC_MS2[ Unstack_DF_TIC_MS2$Protein == proteins[m], ]
  temp1 <- as.data.frame(temp1)
  write.csv(temp1, file = ((paste(proteins[m], "MS2_TIC_Normalzied_Result.CSV", sep = "_"))))
}
setwd(dir = "D:/Quan10273")

# report the mean variation value for each of the models gridExtra is a 
# package for writing tables in pdf format

#library(gridExtra)
#pdf("D:/Quan10273/test.pdf")
#a <- mean(Unstack_DF_TIC$variation)
#b <- mean(Unstack_DF$variation)
#c <- mean(Unstack_Background_Reduced_DF$variation)
#d <- mean(Unstack_MS2_DF$variation)
#e <- mean(Unstack_DF_TIC_MS2$variation)
#variations <- c(a, b, c, d, e)
#grid.table(variations)
#dev.off()




# plot controls to each other


#pdf(paste(samples[n],".pdf",sep="")) 
#plot(DF3$time,log2(DF3[, n]), xlab = "Time / Minute", ylab = "Log (concentartion)", main = c("Half life time (in minutes) for", as.character(samples[n]), "is", round(Half_Time, 0)))

#dev.off() 
#pdf(paste(samples[n],"percentage.pdf",sep="_"))
#require(gridExtra)
#signal <- as.data.frame(log2(DF3[, n]))
#i <- length(unique(DF3$time))
#intercept1 <- rep(intercept, times = i)
#result <- round(2^(as.numeric(signal$`log2(DF3[, n])`)-as.numeric(intercept1))*100, 3)
#percentage <- data.frame("Time / minutes" = DF3$time, "Percentage" = result)
#grid.table(percentage)
#dev.off() 

#Mod1 <- lm(Unstack_DF[,1] ~ Unstack_DF[, 2])
#require(broom)
#tidy_Mod1 <- tidy(Mod1)
#slope <- tidy_Mod1[2, 2]
#intercept <- tidy_Mod1[1, 2]
#x11(20,20,12)