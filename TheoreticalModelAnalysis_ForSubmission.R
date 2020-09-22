#
###########################################################################################################################################
###########################################################################################################################################
###   Script Details:                                                                                                                   ###
###   Plotting model results for Singh, A. A.F., Agrawal. 2020. Sex-Specific Variance in Fitness and the Efficacy of Selection          ###
###   Script written by Amardeep Singh (amardeep.singh[at]utoronto.ca)                                                                  ###
###   This script will plot the result of a model used to assess the relationship between condition dependency of fitness in males      ###
###     and the efficacy of selection                                                                                                   ###
###########################################################################################################################################
###########################################################################################################################################

## Theoretical Model Plotting 
# NOTE: YOU WILL NEED TO SPECIFY THE PATH TO DATA FILES BELOW BY REPLACING "/FILE/PATH/TO/" WITH YOUR ACTUAL FILE PATHS 

# Here I evaluate and plot the results of Equation 6 in the paper for a set of values of variance in fitness among females (sigma) and a set of 
#       alpha values (strength of selection in male relative to females)
# Plots generated with this script: Figure S1 and Figure S2
# NOTE: Plots will not look exactly like those in the paper (the data is identical, but I processed images to change fonts, labels etc.)

rm(list=ls())

##  Loading packages
require(RColorBrewer)

# Note on variable meanings
#v = variance in female fitness
#a = alpha
#b = beta
#f = f (or phi) (relative contribution of condition versus stochastic processes)

# Range of f values over which to evaluate 
range = c(0,1/3,2/3)
x.range = c(0.5,3.0) # Range of values of beta

# Defining Plotting Settings 

# Plotting parameters 
par(mfrow=c(2,2))
par(mai=c(0,0,0,0), family = "Times")
par(cex = 1)
par(mar = c(2,10,0,0) + 0.4)
par(oma = c(4,4,0.5,10))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
par(mai = c(0.5,0.1,0.1,0.1))

# Colour palette settings and line type
palette = c("#046C9A", "#FF0000", "#F98400", "#4DAF4A")
line.type = rep(1,length(range))
x=rep(1:1, each=21)

# Panel A: Female Variance in fitness = 2
VarF = 2
plot(0,0, axes = FALSE, type="n", ylim=c(0.6,2), las = 2, xlim=x.range, xlab = "", ylab = "", main = "")
for (i in 1:length(range)) {
  curve(((2 + VarF)*(1+x)/(4 + VarF*(2 - range[i] + range[i]*x^2))), xlim=c(0.5,(x.range[2]) + 1.5), add=TRUE, lty = line.type[i], lwd=5, col=palette[i])
}
axis(2, at = seq(from = 0.6,to = 2,by = 0.2),cex.axis = 1.5, las = 2)
abline(h=1, lwd = 2)
box(col = "black")

# Panel B: Female Variance in fitness = 5
VarF = 5
plot(0,0, axes = FALSE, type="n", ylim=c(0.6,2), las = 1, xlim=x.range,xlab = "", ylab = "", main = "") #I used this plot to set the axis limits (mostly for the y axis)
for (i in 1:length(range)) {
  curve(((2 + VarF)*(1+x)/(4 + VarF*(2 - range[i] + range[i]*x^2))), xlim=c(0.5,(x.range[2]) + 1.5), add=TRUE, lty = line.type[i], lwd=5, col=palette[i])
}
abline(h=1, lwd = 2)
box(col = "black")

# Panel C: Female Variance in fitness = 10
par(mai=c(0.1,0.1,0.5,0.1))
#Plot for a variance in female fitness value of 10
VarF = 10
plot(0,0, axes = FALSE, type="n", ylim=c(0.6,2), las = 1,  xlim=x.range,xlab = "", ylab = "", main = "") #I used this plot to set the axis limits (mostly for the y axis)
for (i in 1:length(range)) {
  curve(((2 + VarF)*(1+x)/(4 + VarF*(2 - range[i] + range[i]*x^2))), xlim=c(0.5,(x.range[2]) + 1.5), add=TRUE, lty = line.type[i], lwd=5, col=palette[i])
}
axis(1, at = seq(from = 0.5,to = 3,by = 0.5), cex.axis = 1.5, las = 1)
axis(2, at = seq(from = 0.6,to = 2,by = 0.2), cex.axis = 1.5, las = 2)
abline(h=1, lwd = 2)
box(col = "black")

# Panel D: Female Variance in fitness = 1000
VarF = 1000
plot(0,0, axes = FALSE, type="n", ylim=c(0.6,2), las = 1, xlim=x.range,xlab = "", ylab = "", main = "") #I used this plot to set the axis limits (mostly for the y axis)
for (i in 1:length(range)) {
  curve(((2 + VarF)*(1+x)/(4 + VarF*(2 - range[i] + range[i]*x^2))), xlim=c(0.5,(x.range[2]) + 1.5), add=TRUE, lty = line.type[i], lwd=5, col=palette[i])
}
axis(1, at = seq(from = 0.5,to = 3,by = 0.5), cex.axis = 1.5, las = 1)
abline(h=1, lwd = 2)
box(col = "black")

# Adding a legend to plot 
#par(fig = c(0, 1, 0, 1), oma = c(1, 1, 1, 5), mar = c(0, 0, 0, 0), new = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right",legend = c("0", "1/3", "2/3"), col = palette, lwd = 5, lty = 1, xpd = TRUE, horiz = FALSE, seg.len=1, bty = 'n')
# xpd = TRUE makes the legend plot to the figure

# Saving final plot 
pdf("FILE/PATH/TO/DATA/Figure 1.pdf", width = 15, height = 15)
lsmeans.plot.by.stock
dev.off()

#
