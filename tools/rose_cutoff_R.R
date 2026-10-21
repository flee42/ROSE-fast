numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}
calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
 	inputVector <- sort(inputVector)
	inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
	slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
	xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
	
	if(drawPlot){  #if TRUE, draw the plot
		plot(1:length(inputVector), inputVector,type="l",...)
		b <- y_cutoff-(slope* xPt)
		abline(v= xPt,h= y_cutoff,lty=2,col=8)
		points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
		abline(coef=c(b,slope),col=2)
		title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
		axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
	}
	return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}
#!/usr/bin/env Rscript
# rose_cutoff_R.R -- compute ROSE cutoff with the AUTHORS own R function.
# usage: Rscript tools/rose_cutoff_R.R <dir with *_signal.tsv>
#
# QA ONLY. This is a copy of ROSE's calculate_cutoff() used to cross-check
# that rose_call_super_R_native.py really reproduces ROSE_callSuper.R.
# Production runs always go through the upstream ROSE_callSuper.R.
args <- commandArgs(trailingOnly=TRUE)
d <- if (length(args) >= 1) args[1] else "."
files <- list.files(d, pattern="_signal", full.names=TRUE)
cat("set\tn\txstar\tcutoff\toverMedian\tfrac_super\n")
for (f in files) {
  tb <- read.delim(f, stringsAsFactors=FALSE)
  v  <- as.numeric(tb$signal)
  r  <- calculate_cutoff(v, drawPlot=FALSE)
  sv <- sort(v)
  xstar <- max(which(sv <= r$absolute))
  nm <- basename(f)
  nm <- substr(nm, 1, nchar(nm) - 11)
  cat(sprintf("%s\t%d\t%d\t%.6f\t%.4f\t%.6f\n",
      nm, length(v), xstar, r$absolute, r$overMedian, mean(v >= r$absolute)))
}
