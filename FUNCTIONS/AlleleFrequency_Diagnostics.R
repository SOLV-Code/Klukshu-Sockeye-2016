# THIS SCRIPT HAS PROTOTYPES AND TESTS FOR R FUNCTIONS THAT SUMMARIZE ALLELE NUMBERS AND FREQUENCIES
# allele.boots() is a subroutine that takes a vector of alleles and does a bootstrap estimate of num.alleles/100samples
# allele.count() is a subroutine that creates a summary of alleles for a locus (2 cols in object created by genepop.read)
# allele.diag() is a wrapper function that applies allele.count to all the loci for each sample group,  and user-specified aggregations



###########################################################################

num.uniques<- function(x){ return(length(unique(na.omit(x))))}



###########################################################################


allele.boots <- function(x,boot.size=1000,n=100){
# x is a vector of alleles

boot.mat <- matrix(sample(x,size=boot.size*n,replace=TRUE), ncol=n,nrow=boot.size)

boot.avg <- apply(boot.mat,MARGIN=1,FUN=num.uniques)


return(c(mean(boot.avg),range(boot.avg)))

}


###########################################################################


allele.count <- function(X,invalid.label="00"){
# X is a matrix with 2 columns, subset from a genepop.read() output


allele.vec.all <- unname(unlist(X))
allele.vec.valid <-  allele.vec.all[allele.vec.all != invalid.label]

out.mat <- matrix(NA,ncol=1, nrow=10,dimnames=list(
				c("n.samples","n.valid","n.invalid","lowest","highest",
					"num.alleles","num.alleles.p100.boots.avg",
					"num.alleles.p100.boots.min", "num.alleles.p100.boots.max","prop.top5")
				,"Value"))

out.mat["n.samples","Value"] <- dim(X)[1]				
out.mat["n.valid","Value"] <- length(allele.vec.valid)
out.mat["n.invalid","Value"] <- length(allele.vec.all) - out.mat["n.valid","Value"]

out.mat[c("lowest","highest"),"Value"] <- range(as.numeric(allele.vec.valid))
out.mat["num.alleles","Value"] <- length(unique(allele.vec.valid ))

out.mat[c("num.alleles.p100.boots.avg","num.alleles.p100.boots.min", "num.alleles.p100.boots.max"),
					"Value"] <- allele.boots(allele.vec.valid,boot.size=1000,n=100)


out.mat[,"Value"] <- round(out.mat[,"Value"],3)
					
freq.vec <- sort(table(allele.vec.valid ),decreasing=TRUE)


out.mat["prop.top5","Value"] <- round(sum(freq.vec[1:min(length(freq.vec),5)])/ out.mat["n.valid","Value"],3)


return(list(summary.table=out.mat,freq.vec = freq.vec))

}

###########################################################################

allele.diag <- function(X,report.groups=NULL,loci.list=NULL){
# X is a matrix component from  object create by genepop.read(), or some modified version of it
# report.groups is a list, with each list element a vector of sample groups to be combined and the element name used as a label

grp.list <- unique(X[,"GROUP"])

# do summary for 1st locus across all groups
all.summary <- allele.count(X[, grep(loci.list[1] , dimnames(X)[[2]])])


# use it as a spec for the output storage
out.summary.array  <- array(NA,dim=c(dim(all.summary$summary.table)[1],length(loci.list), length(grp.list)+length(report.groups)+1),
							dimnames=list(dimnames(all.summary$summary.table)[[1]],loci.list,c("All",names(report.groups),grp.list) ))

out.freq.list <- list(TEMP=NA)

for(locus.use in loci.list){

locus.idx <-  grep(locus.use , dimnames(X)[[2]])
 
	counts.tmp <- allele.count(X[,locus.idx])
	out.summary.array[ ,locus.use,"All"] <- counts.tmp$summary.table
	
	# create a template to store allele freq for this locus
	out.freq.list[[locus.use]] <- matrix(0,ncol=length(counts.tmp$freq.vec),nrow= dim(out.summary.array)[3],dimnames=list(
								c("All",names(report.groups),grp.list),sort(names(counts.tmp$freq.vec)) )   )
	
	# add "all" freq vector
	out.freq.list[[locus.use]]["All", names(counts.tmp$freq.vec)] <- counts.tmp$freq.vec
	 
	 
	for(rep.grp.use in names(report.groups)){
		rep.grp.idx <- X[,"GROUP"] %in% report.groups[[rep.grp.use]]
		counts.tmp <- allele.count(X[rep.grp.idx ,locus.idx])
		out.summary.array[ ,locus.use,rep.grp.use] <- counts.tmp$summary.table 
		out.freq.list[[locus.use]][rep.grp.use, names(counts.tmp$freq.vec)] <- counts.tmp$freq.vec
		}	
	 
	for(grp.use in grp.list){
		grp.idx <- X[,"GROUP"] == grp.use
		counts.tmp <- allele.count(X[grp.idx ,locus.idx])
		out.summary.array[ ,locus.use,grp.use] <- counts.tmp$summary.table 
		out.freq.list[[locus.use]][grp.use, names(counts.tmp$freq.vec)] <- counts.tmp$freq.vec
	
		} # end looping through groups
	
} # end looping through loci

out.freq.list <- out.freq.list[-1] # drop the placeholder element 
return(list(Table=out.summary.array,Freq=out.freq.list ))

}


