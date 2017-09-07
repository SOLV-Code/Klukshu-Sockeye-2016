
################################################################################
# compare two sets of size samples


TwoSample.comp.test <- function(sample1,sample2,sample.labels=c("Sample 1","Sample 2"),file.label=	"Sample_Comparison.csv", boot.size=100,boot.prop=0.9){

library(moments)

sample.mat.labels <- c("n","Mean","Median","SD", "Skewness","S","Kurtosis","K")  #,"S or K > 1.96?")
sample.mat <- matrix(NA,nrow=2,ncol=length(sample.mat.labels),dimnames=list(sample.labels,			
					sample.mat.labels))

sample.mat[,"n"] <- c(length(na.omit(sample1)),length(na.omit(sample2)))
sample.mat[,"Mean"]  <- c(mean(sample1,na.rm=TRUE),mean(sample2,na.rm=TRUE))
sample.mat[,"Median"]  <- c(median(sample1,na.rm=TRUE),median(sample2,na.rm=TRUE))
sample.mat[,"SD"]  <- c(sd(sample1,na.rm=TRUE),sd(sample2,na.rm=TRUE))
sample.mat[,"Skewness"]  <- c(skewness(sample1,na.rm=TRUE),skewness(sample2,na.rm=TRUE))
sample.mat[,"Kurtosis"]  <- c(kurtosis(sample1,na.rm=TRUE),kurtosis(sample2,na.rm=TRUE))
sample.mat[,"S"] <- sample.mat[,"Skewness"] /sample.mat[,"SD"] 
sample.mat[,"K"] <- sample.mat[,"Kurtosis"] /sample.mat[,"SD"] 
#sample.mat[,"S or K > 1.96?"] <- sample.mat[,"S"] >1.96 | sample.mat[,"K"] > 1.96

t.test.out <- t.test(sample1,sample2,na.action="na.omit")
conf.int.out <- t.test.out$conf.int
attributes(conf.int.out ) <- NULL


p.values.boot <- rep(NA,boot.size)
for(i in 1:boot.size){
	sub1 <- sample(na.omit(sample1),length(na.omit(sample1))*boot.prop,replace=FALSE)
	sub2 <- sample(na.omit(sample2),length(na.omit(sample2))*boot.prop,replace=FALSE)
	p.values.boot[i] <- t.test(sub1,sub2,na.action="na.omit")$p.value
	
	
}
			

out.list <- list( p.value= t.test.out$p.value,
				 prop.sig.p.values = sum(p.values.boot<=0.05)/boot.size,
				conf.int.sample.diff= conf.int.out,
				sample.summary = sample.mat  #, bootstrapped.p.values = p.values.boot
				)
			
#write.csv(sample.mat, file=file.label,row.names=FALSE)


return(out.list)

}



sample.summary <- function(X){			
#X is a vector of numbers			
library(moments)				

out.labels <- c("n","Mean","Median","SD", "Skewness","S","Kurtosis","K")  #,"S or K > 1.96?")
out.vec <- rep(NA,length(out.labels))
names(out.vec) <- out.labels



out.vec["n"] <- length(na.omit(X))
out.vec["Mean"]  <- mean(X,na.rm=TRUE)
out.vec["Median"]  <- median(X,na.rm=TRUE)
out.vec["SD"]  <- sd(X,na.rm=TRUE)
out.vec["Skewness"]  <- skewness(X,na.rm=TRUE)
out.vec["Kurtosis"]  <- kurtosis(X,na.rm=TRUE)
out.vec["S"] <- out.vec["Skewness"] /out.vec["SD"] 
out.vec["K"] <- out.vec["Kurtosis"] /out.vec["SD"] 

out.vec <- round(out.vec,3)

return(out.vec)

}


FourSample.comp.test <- function(sample.list, file.label="OUTPUT/Sample_Comparison", boot.size.in=1000, boot.prop.in=0.9){

library(moments)		
		
tmp <- sample.summary(sample.list[[1]])


sample.summary.mat <- matrix(NA, nrow=4,ncol=length(tmp),dimnames=list(names(sample.list),names(tmp)))

sample.summary.mat[1,] <- sample.summary(sample.list[[1]])
sample.summary.mat[2,] <- sample.summary(sample.list[[2]])
sample.summary.mat[3,] <- sample.summary(sample.list[[3]])
sample.summary.mat[4,] <- sample.summary(sample.list[[4]])


boots.mat<- matrix(NA,ncol=4,nrow=4,dimnames=list(names(sample.list),names(sample.list)))

sets.list <- list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))

for(i in 1:length(sets.list)){
		boots.mat[sets.list[[i]][1],sets.list[[i]][2]] <- TwoSample.comp.test(sample.list[[sets.list[[i]][1]]], sample.list[[sets.list[[i]][2]]],
					sample.labels = names(sample.list[c(sets.list[[i]][1],sets.list[[i]][2])]),boot.size=boot.size.in,boot.prop=boot.prop.in)[["prop.sig.p.values"]]
	}
					

boots.mat <- round(boots.mat*100,0)
boots.mat <- boots.mat[-4,-1]


write.csv(boots.mat,file=paste(file.label,"_P_Value_Bootstrap.csv",sep=""))

write.csv(sample.summary.mat,file=paste(file.label,"_SampleSummary.csv",sep=""))			

out.list <- list(sample.summary = sample.summary.mat,boots.mat=boots.mat)

return(out.list)


}


##############################################################################
#  Chi-Squared Test

chi.boot <- function(X,boot.size=1000, boot.prop=0.9,continuity.correct=FALSE){
# X is a matrix with 2 columns and two levels in each column
# bootstrapping is done on the first column
	X <- X[order(X[,1]),] # SORT baseD on col 1
	
	library(moments)
	
	prop.test.out <- prop.test(table(X ),correct=continuity.correct)

	x1.levels <- unique(X[,1])
	
	x1.L1.mat <- X[ X[,1]== x1.levels[1],]
	L1.n <- dim(x1.L1.mat)[1]
	x1.L2.mat <- X[ X[,1]== x1.levels[2],]
	L2.n <- dim(x1.L2.mat)[1]

	p.values.boot <- rep(NA,boot.size)
	for(i in 1:boot.size){
		L1.sub <- x1.L1.mat[sample(L1.n,L1.n*boot.prop),]
		L2.sub <- x1.L2.mat[sample(L2.n,L2.n*boot.prop),]
		x.sub <- rbind(L1.sub, L2.sub)
		p.values.boot[i] <- prop.test(table(x.sub),correct=continuity.correct)$p.value	
	}

out.list <- list( est= prop.test.out$estimate,
				p.value= prop.test.out$p.value,
				 prop.sig.p.values = sum(p.values.boot<=0.05)/boot.size,
				conf.int.sample.diff= prop.test.out$conf.int,
				sample.summary = table(X ) #, bootstrapped.p.values = p.values.boot
				)

return(out.list)

}





##############################################################################
#  ALTERNATIVE WEIGHTED MEANS



alt.means.fun <- function(X,Y,x.col="Perc River"){
# X is one of the weekly composition matrices created earlier (% fem  by week
# Y is a weights table created earlier (weights by week, based on run size, tag number etc.)

wk.labels <- paste("W",28:41,sep="")
early.labels  <- paste("W",28:33,sep="")
late.labels  <- paste("W",35:41,sep="")
early.weeks <- as.character(28:33)
late.weeks <- as.character(35:41)

means.out <- list(
	WtMean_Run = weighted.mean(X[,x.col], Y[wk.labels ,"wt_Run"],na.rm=TRUE),
	WtMean_TagRatio = weighted.mean(X[,x.col], Y[wk.labels ,"wt_TF"]*Y[wk.labels ,"Tag_Fates"],na.rm=TRUE),
	WtMean_n = weighted.mean(X[,x.col], Y[wk.labels ,"Tag_Fates"],na.rm=TRUE),
	RawMean = mean(X[,x.col],na.rm=TRUE)
	)
return(means.out)

}



















