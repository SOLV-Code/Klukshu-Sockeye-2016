

# function to create smoothed average (7 day moving average)

mv.avg <- function(x,len){
	f_len <- rep(1/len,len) # weights for weighted avg
	ma_x <- filter(x, f_len, sides=2)
	ma_x
	}





# function to prepare data for cluster analysis
# (convert each daily count into 1 record of a return day for each fish)
counts.reorg <- function(X,cols){
	out.list <- list(Label="Klukshu Weir Counts Reorganized") # start storage list
	for(col in cols){	
		X.sub <- X[,c("Day_Num",as.character(col))]
		out.vec <- NA # start storage vector
		for(day in X.sub[,1]){
			num.reps <- round( X.sub[X.sub[,1]==day,2])
			#print(col);print(day);print(num.reps)
			if(num.reps>0 & !is.na(num.reps)){out.vec <-c(out.vec,rep(day,num.reps))}
			}

		out.vec <- out.vec[-1] # remove first entry
		#print(out.vec)
		eval(parse(text=paste("out.list <- c(out.list,list(",paste("Yr",col,"=",sep=""), "out.vec))")))
	} # end looping through cols
	out.list
} # end function





########################################
# CLUSTER ANALYSIS VERSION 1
# function to do cluster analysis using the {mclust} package

mclust.fun <- function(dat, G.set=2,plot.layout="panels",x.range=c(1,152),x.markers="default"){



clust.fit <- Mclust(dat,G=G.set)

out.mat <-matrix(NA,nrow=4,ncol=G.set,dimnames=list(c("Mean","StDev","Pop1p90","PercPop1atp90"),paste("Pop",1:G.set,sep="")))

out.mat["Mean",] <- clust.fit$parameters$mean
out.mat["StDev",] <- sqrt(clust.fit$parameters$variance$sigmasq)
out.mat["Pop1p90","Pop1"] <- qnorm(.9,out.mat["Mean","Pop1"],out.mat["StDev","Pop1"])

print("Fits-----------")
print(out.mat)


if(G.set==2){
	# calculate proportion of p1 on day when 90% of Pop has passed by

	#TBI

} # end if G.set=2






#x.vals <- x.range[1]:x.range[2]
#plot(x.vals,dnorm(x.vals,,),



if(plot.layout=="panels"){layout(matrix(1:6,ncol=3,byrow=TRUE))}


# Plot 1
print("starting Plot 1")
#hist(dat,breaks=30,axes=FALSE,xlab="",main="Histogram of Data")
plot(density(dat), xlim=c(1,152),main="Observed Density",axes=FALSE,xlab="")

	if(x.markers=="default"){
			axis(side=1, at=c(1,31,62,93,123,152),labels=FALSE,lwd.ticks=2)
			axis(side=1, at=c(15,45,76,107,137),labels=c("Jun","Jul", "Aug", "Sep", "Oct"))
		}


# Plot 2
print("starting Plot 2")
plot(clust.fit, what="density", axes=FALSE,xlab="",main="")
title(main="Model Density")
	if(x.markers=="default"){
			axis(side=1, at=c(1,31,62,93,123,152),labels=FALSE,lwd.ticks=2)
			axis(side=1, at=c(15,45,76,107,137),labels=c("Jun","Jul", "Aug", "Sep", "Oct"))
		}

# Plot 3
print("starting Plot 3")
plot(clust.fit, what="uncertainty",axes=FALSE,xlab="",main="")
title(main="Model Uncertainty")
	if(x.markers=="default"){
			axis(side=1, at=c(1,31,62,93,123,152),labels=FALSE,lwd.ticks=2)
			axis(side=1, at=c(15,45,76,107,137),labels=c("Jun","Jul", "Aug", "Sep", "Oct"))
		}


# Plot 4 - MIXTURE MODEL
print("starting Plot 4")

x.vals <- x.range[1]:x.range[2]
y.range <- c(0,max(dnorm(out.mat[1,],out.mat[1,],out.mat[2,])))

plot(x.vals,dnorm(x.vals,clust.fit$parameters$mean[1],sqrt(clust.fit$parameters$variance$sigmasq[1])),
	type="l",axes=FALSE,xlab="",ylab="Density",main="Mixture Model",col="darkblue",ylim=y.range)
for(i in 2:G.set){lines(x.vals,dnorm(x.vals,clust.fit$parameters$mean[i],sqrt(clust.fit$parameters$variance$sigmasq[i])))}
	if(x.markers=="default"){
			axis(side=1, at=c(1,31,62,93,123,152),labels=FALSE,lwd.ticks=2)
			axis(side=1, at=c(15,45,76,107,137),labels=c("Jun","Jul", "Aug", "Sep", "Oct"))
		}
		abline(v=out.mat["Pop1p90","Pop1"],col="red",lwd=2)
		
# Plot 5 - SCALED MIXTURE MODEL
print("starting Plot 5")


x.vals <- x.range[1]:x.range[2]

pop.counts <- table(clust.fit$classification) 
#print("classes");print(pop.counts)

if(length(pop.counts)==G.set){

peak.vals <- dnorm(out.mat[1,],out.mat[1,],out.mat[2,]) * pop.counts
# this multiplies each density function with the total number of obs in each classification
# e.g. density(Pop1) * count(Pop1)

y.range <- c(0,max(peak.vals))


plot(x.vals,dnorm(x.vals,clust.fit$parameters$mean[1],sqrt(clust.fit$parameters$variance$sigmasq[1]))*pop.counts[1],
	type="l",axes=FALSE,xlab="",ylab="Density",main="Scaled Mixture Model",col="darkblue",ylim=y.range)
for(i in 2:G.set){lines(x.vals,dnorm(x.vals,clust.fit$parameters$mean[i],sqrt(clust.fit$parameters$variance$sigmasq[i]))*pop.counts[i])}
	if(x.markers=="default"){
			axis(side=1, at=c(1,31,62,93,123,152),labels=FALSE,lwd.ticks=2)
			axis(side=1, at=c(15,45,76,107,137),labels=c("Jun","Jul", "Aug", "Sep", "Oct"))
		}
		abline(v=out.mat["Pop1p90","Pop1"],col="red",lwd=2)
		
}
	
if(length(pop.counts)!=G.set){plot(1:10,1:10,bty="n", type="n",axes=FALSE,xlab="",ylab="")}


	
		
		
print(summary(clust.fit, parameters = TRUE))
#print(clust.fit$classification)

out.mat


}  # end mclust.fun



########################################
# CLUSTER ANALYSIS VERSION 1
# function to do cluster analysis using the kmeans function

kmeans.fun <- function(dat, G.set=2,plot.layout="panels",x.range=c(1,152),x.markers="default"){

kmeans.fit <- kmeans(dat,centers=G.set)
out.mat <-matrix(NA,nrow=4,ncol=G.set,dimnames=list(c("Centers","TBD","TBD","TBD"),paste("Pop",1:G.set,sep="")))

out.mat["Centers",] <- sort(kmeans.fit$centers)

#print(out.mat)

if(plot.layout=="panels"){layout(matrix(1:4,ncol=2,byrow=TRUE))}


x.vals <- x.range[1]:x.range[2]
plot(x.vals, rep(1,length(x.vals)),type="n",axes=FALSE,xlab="",ylab="",main="Peaks (kmean)")
abline(v=out.mat["Centers",],col="darkblue",lwd=2)
	if(x.markers=="default"){
			axis(side=1, at=c(1,31,62,93,123,152),labels=FALSE,lwd.ticks=2)
			axis(side=1, at=c(15,45,76,107,137),labels=c("Jun","Jul", "Aug", "Sep", "Oct"))
		}




out.mat

} # end kmeans.fun









## TIMING PLOT
# function for basic timing plots

plot.timing <- function(x,y, mv.avg=7,x.lab="Date", y.lab="Counts",x.markers="default"){
	# x = timestep
	# y = obs
	# mv.avg = period for moving average
	
	ylims <- pretty(y)
	ylims <- ylims[ylims<max(y)]
	
	plot(x,y,bty="n", axes=FALSE,type="l",col="dodgerblue",xlab=x.lab, ylab=y.lab)
	lines(x,mv.avg(y,mv.avg),col="red",lwd=3)
	axis(side=2,at=)
	if(x.markers=="default"){
			axis(side=1, at=c(1,31,62,93,123,152),labels=FALSE,lwd.ticks=2)
			axis(side=1, at=c(15,45,76,107,137),labels=c("Jun","Jul", "Aug", "Sep", "Oct"))
		}		
	}


	
plot.timing.sm <- function(x,y, mv.avg=7,x.lab="Date", y.lab="7d Avg",x.markers="default"){
	# x = timestep
	# y = obs
	# mv.avg = period for moving average
	
	plot.vals <- mv.avg(y,mv.avg)
	ylims <- pretty(plot.vals)
	ylims <- ylims[ylims<max(plot.vals)]
	
	plot(x,plot.vals,bty="n", axes=FALSE,type="l",col="red",lwd=3,xlab=x.lab, ylab=y.lab)
	axis(side=2,at=)
	if(x.markers=="default"){
			axis(side=1, at=c(1,31,62,93,123,152),labels=FALSE,lwd.ticks=2)
			axis(side=1, at=c(15,45,76,107,137),labels=c("Jun","Jul", "Aug", "Sep", "Oct"))
		}		
	}











