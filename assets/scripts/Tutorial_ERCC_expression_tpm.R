#!/usr/bin/env Rscript

#Jason Walker, jason.walker[AT]wustl.edu
#Malachi Griffith, mgriffit[AT]wustl.edu
#Obi Griffith, obigriffith[AT]wustl.edu
#The McDonnell Genome Institute, Washington University School of Medicine

#R tutorial for Informatics for RNA-sequence Analysis workshops

library(ggplot2)

args <- commandArgs(TRUE)
filename <- args[1]

data = read.delim(filename)

data$logTPM  = log2(data$TPM + 1)
data$logConc= log2(data$Concentration)

tpm_model <- lm(logTPM ~ logConc, data=data)
tpm_r_squared = summary(tpm_model)[['r.squared']]
tpm_slope = coef(tpm_model)["logConc"]

p = ggplot(data, aes(x=logConc, y=logTPM))
p = p + geom_point(aes(shape=Label,color=Label))
p = p + geom_smooth(method=lm) + annotate('text', 5, -3, label=paste("R^2 =", tpm_r_squared, sep=' ')) + annotate('text', 5, -4, label=paste("Slope =", tpm_slope, sep=' '))
p = p + xlab("Expected concentration (log2 scale)") + ylab("Observed TPM (log2 scale)")

pdf('Tutorial_ERCC_expression_tpm.pdf')
print(p)
dev.off()

