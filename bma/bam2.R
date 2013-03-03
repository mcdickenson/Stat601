# Exploring ensembleBMA

install.packages('ensembleBMA')
library(ensembleBMA)


# Example - p22
data(ensBMAtest)

ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo") 
obs <- paste("PCP24","obs", sep = ".")
ens <- paste("PCP24", ensMemNames, sep = ".")

prcpTestData <- ensembleData(forecasts = ensBMAtest[,ens], 
	dates = ensBMAtest[,"vdate"], 
	observations = ensBMAtest[,obs], 
	station = ensBMAtest[,"station"], 
	forecastHour = 48, 
	initializationTime = "00")

prcpTestFit <- ensembleBMAgamma0( prcpTestData, trainingDays = 30)

# My Example

curve(dgamma(x,2,scale=1/5))