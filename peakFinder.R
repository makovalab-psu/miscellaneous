par(mfrow=c(1,2))

setwd("/Users/alice/Desktop/projects/histos")
library("fda")

files <- list.files(pattern=".kmer_dist")
for (input in files){
  h = read.table(input)
  fitDat <- data.frame(x=h$V1, y=log(h$V2))
  smoothingSplineHist = smooth.spline(fitDat$y, spar=0.35)  

  plot(fitDat$x, fitDat$y,xlim=c(0,200),main=input, xlab="Number of times a distinct k-mer occurs",       ylab="Number of k-mers")

  bspline_basis <- create.bspline.basis(c(1,max(h$V1)), norder=4, 
                                     breaks=seq(1,200,length=20))
  data_functional <- Data2fd(fitDat$x, fitDat$y, bspline_basis, 2, 10e-20)
  deriv=eval.fd(seq(1,max(h$V1)), data_functional, 1)
  lines(data_functional,col="green")

  #abline(v=min(which(smoothingSplineDiff$y>cummin(smoothingSplineDiff$y))),col="red")
  #abline(v=min(which(smoothingSplineHist$y>cummin(smoothingSplineHist$y))),col="blue")

  if (min(which(deriv>0))==Inf) {
    abline(v=min(which(abs(deriv)>cummin(abs(deriv)))),col="violet")
  } else {
    abline(v=min(which(deriv>0)),col="brown") #looking for first positive value
  }


  #RIGHT PLOT
  fitDat <- data.frame(smoothingSplineHist$x,smoothingSplineHist$y)
  shifted_fitDat <- fitDat[-1,] #remove first observation from the histogram
  fitDat <- fitDat[-nrow(fitDat),] #remove last observation from the histogram

  differences_between_neighbors<-fitDat$smoothingSplineHist.y-shifted_fitDat$smoothingSplineHist.y #store differences between neighboring points in the histogram, this will reveal the "drop" in the differences or in other words "bump" in the coverage
  
  plot(differences_between_neighbors,xlim=c(0,200))
  smoothingSplineDiff = smooth.spline(differences_between_neighbors, spar=0.35)
  
  #abline(v=min(which(smoothingSplineDiff$y>cummin(smoothingSplineDiff$y))),col="red")
  #abline(v=min(which(smoothingSplineHist$y>cummin(smoothingSplineHist$y))),col="blue")

}