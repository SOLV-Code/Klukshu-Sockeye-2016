# THIS SCRIPT HAS PROTOTYPES AND TESTS FOR GENEPOP-RELATED R FUNCTIONS
# - genepop.read(): reads in a tab-delimited text file that was manually modified from a GENEPOP file
# - genepop.plot.basicbar(): uses the data.frame created by genepop.read() to create a basic frequency plot (hist)
# genepop.plot.basicscatter(): uses the data.frame created by genepop.read() to create a basic scatter plot 



# mini function
genotype.merge <- function(x){x.out <- paste(x,sep="",collapse="");return(x.out)}


###################################################################
# genepop.read()
###################################################################


genepop.read <- function(X,ncode=2,label.comp=c("SAMPLE","DAY","SAMPLE_ID_DAY")){
# - X is the path for a tab-delimited file created by loading the genepop file into excel, then cleaning, then saving as tab-delimited
# - ncode 	is an integer indicating the number of characters used to code an allele. (i.e if 27/58)
# - for now this is very specific to outputs from PBS lab (assumes diploidy etc)
# TO DO:
# - fix allele range calc

raw.table <-  read.table(X,sep="\t", header=TRUE, colClasses = "character")
loci.idx <- grep("loc_",dimnames(raw.table)[[2]])
loci.list <- sort(dimnames(raw.table)[[2]][loci.idx]) # assuming that current order of columns is not meaningful in any way -> is it?
loci.labels <- paste(sort(rep(loci.list,2)),1:2,sep="_")

loci.split.list <- paste(sort(rep(loci.list,2)),1:2,sep="_")

labels.list <-  dimnames(raw.table)[[2]][ !( (1:dim(raw.table)[2]) %in% loci.idx)]
num.records <- dim(raw.table)[1]




out.table <- data.frame(matrix(NA, nrow=num.records,ncol=1+length(c(labels.list, loci.list,loci.list)),
				dimnames=list(1:num.records,c("GROUP",labels.list,loci.split.list)) ))

out.table[,labels.list]<-raw.table[,labels.list]

for(locus in loci.list){
		out.table[,paste(locus,1,sep="_")] <- substr(raw.table[,locus],start=1,stop=ncode)
		out.table[,paste(locus,2,sep="_")] <- substr(raw.table[,locus],start=ncode+1,stop=ncode+ncode)
		}


out.table <- cbind(ID=NA,ID_Split=NA,out.table)
out.table[,"ID"] <- paste(out.table[,label.comp[1]],out.table[,label.comp[2]],out.table[,label.comp[3]],sep="_")
out.table[,"ID_Split"] <- paste(out.table[,label.comp[1]],out.table[,label.comp[2]],out.table[,label.comp[3]],sep=" ")
out.table[,"GROUP"] <- out.table[,"SAMPLE"]
out.table <- cbind( Full.Genotype =  paste("GT",apply(out.table[,loci.labels],MARGIN=1,FUN=genotype.merge)),out.table)



allele.range <- range(as.numeric(unlist(out.table[,loci.split.list])),na.rm=TRUE)


out.list <- list(mat=out.table, loci.list=loci.list, loci.labels=loci.labels, allele.range=allele.range)


return(out.list)

}



###################################################################
# allele.coding.trim()
###################################################################

allele.coding.trim <- function(X){
# X is the matrix component of an object created by genepop.read()
# for now only handles 3 digit coding -> 2 digit coding
loci.idx <- grep("loc",dimnames(X)[[2]] )
for(idx in loci.idx){
	X[,idx] <- substr(X[,idx],start=2,stop=3)
	}
return(X)
}










###################################################################
# genepop.write()
###################################################################
# X.full is the matrix component of an object created by genepop.read
# e.g. klukshu.sox$mat
# loci.list is the klukshu.sox$loci.list
# label.use = if 'GROUP' then all records in a POP have identical labels, 
#				if "ID" then have unique labels (single text string with under scores)
#				if "ID_Split" then have unique labels with blank spaces (like original genepop files from PBS)

# --------------------------------------
genepop.write <- function(X.full, loci.list, file.label = "test.txt",
					file.description = "Short description of genepop file goes here",
					label.use = "GROUP" ){


group.list <- unique(X.full[,"GROUP"])


sink(file.label)

cat(file.description,"\n")



for(i in 1:length(loci.list)){
	cat(loci.list[i],"\n")
}




for(grp in group.list){

X <- X.full[X.full[,"GROUP"]==grp,]

cat("POP","\n")



for(i in 1:dim(X)[1]){
	

	cat(paste(X[i,label.use],    # grp
			",",		
			paste(X[i,c(paste(loci.list[1],1,sep="_"),paste(loci.list[1],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[2],1,sep="_"),paste(loci.list[2],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[3],1,sep="_"),paste(loci.list[3],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[4],1,sep="_"),paste(loci.list[4],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[5],1,sep="_"),paste(loci.list[5],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[6],1,sep="_"),paste(loci.list[6],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[7],1,sep="_"),paste(loci.list[7],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[8],1,sep="_"),paste(loci.list[8],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[9],1,sep="_"),paste(loci.list[9],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[10],1,sep="_"),paste(loci.list[10],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[11],1,sep="_"),paste(loci.list[11],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[12],1,sep="_"),paste(loci.list[12],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[13],1,sep="_"),paste(loci.list[13],2,sep="_"))],collapse=""),
			paste(X[i,c(paste(loci.list[14],1,sep="_"),paste(loci.list[14],2,sep="_"))],collapse=""),
				sep="  ",collapse="  "),"\n")

	
	} # end looping through records for that group

} # end looping through groups

sink()

} # end of genepop.write()
#------------------------------------------




###################################################################
# genepop.incomplete.filter()
###################################################################

# for now, remove any records with more than n missing alleles 
# also create before and after summary
# Future extension: check and remove at loci level


genepop.incomplete.filter<- function(X, incomplete.code="00", incomplete.limit=6){
# X is an object created by genepop.read

before.summary <- genepop.summary(X)

loci.idx <- grep("loc_",dimnames(X$mat)[[2]])

X$mat <- X$mat[rowSums(X$mat[,loci.idx]==incomplete.code) < incomplete.limit,]


after.summary <- genepop.summary(X)


track.mat <- cbind(before.summary[,1:2],NA,NA,NA)
names(track.mat)<- c("Group","Rec.Before","Rec.After","Num.Dropped","Perc.Dropped")

match.idx <- dimnames(track.mat)[[1]] %in% dimnames(after.summary )[[1]]

track.mat[match.idx,"Rec.After"] <- after.summary[,"Samples"]
track.mat[match.idx,"Num.Dropped"] <- track.mat[match.idx,"Rec.Before"] -track.mat[match.idx,"Rec.After"]
track.mat[match.idx,"Perc.Dropped"] <- round(track.mat[match.idx,"Num.Dropped"] / track.mat[match.idx,"Rec.Before"]*100,0)



return(list(genepop.obj=X,incomplete.tracker= track.mat))

}




###################################################################
# genepop.summary()
###################################################################

# for now just count number of samples in each group, and if there's a YEAR
# column, list the sample years.


genepop.summary <- function(X,group.label="GROUP"){
# X is an object created by genepop.read
grp.list <- unique(X$mat[,group.label])
out.mat <- data.frame(table(X$mat[,group.label]),     stringsAsFactors=FALSE)
dimnames(out.mat)[[2]] <- c("Group","Samples")
if("YEAR" %in% dimnames(X$mat)[[2]]){

	out.mat <- cbind(out.mat,SampleYears=NA)

	for(grp in grp.list){
	out.mat[out.mat[,"Group"]==grp,"SampleYears"] <-   paste(unique( X$mat[X$mat[,group.label]==grp,"YEAR"]),collapse="/")
	}	
	
}
out.mat
} # end genepop.summary




genepop.subsample <- function(X,drop.prop=0.1, group.label="GROUP",trace=TRUE){
# X is an object originally created by genepop.read (but may be subsequently modified)
# drop.prop is the proportion of sample to be dropped (will be rounded to the nearest numer of records
# capped at a max 80%

grp.list <- unique(X$mat[,group.label])

for(grp in grp.list){

grp.idx <- which(X$mat[,group.label]==grp) # gives locations of TRUE values in the vector
drop.num <- round(length(grp.idx) * 0.1)
drop.idx <- sample(grp.idx,drop.num)
if(trace){
print("--------------")
print(grp)
print(paste("Records =",length(grp.idx)))
print(paste("Dropping =",drop.num))
print(paste("Retained =",length(grp.idx)-length(drop.idx)))
}
X$mat <- X$mat[-drop.idx,,drop=FALSE]
}



return(X)

}




########################################################################################################
# PLOT A PANEL WITH HISTOGRAM OF ALLELE VALUES FOR 1 LOCUS AND 1 GROUP 
# across the 2 chromosomes, if give number or locus names (e.g.  5 or "loc_1b")
# for one of the chromosomes, if give specific label (e.g. "loc_1b_1"

genepop.plot.basicbar <- function(X,locus.sub=1,group.sub=1,plot.title=NULL,x.lim=NULL,breaks.use=30){
# X is a list object created by genepop.read, possibly after modifying the GROUP variable for thr X$mat object
# locus.sub = either a number or a label
# group.sub = either a number or a label
# break.use = either a number, or a vector as per ?hist()

if(is.numeric(group.sub)){group.plot <- unique(X$mat[,"GROUP"])[group.sub]}
if(!is.numeric(group.sub)){group.plot <- group.sub}

if(is.numeric(locus.sub)){locus.plot <- X$loci.list[group.sub]}
if(!is.numeric(locus.sub)){locus.plot <- locus.sub}

if(is.null(plot.title)){plot.title <- paste(group.plot,locus.plot,sep="")}



var.list <- dimnames(X$mat)[[2]]


sub.data <- as.numeric(unlist(X$mat[grep(group.plot,X$mat[,"GROUP"])   ,  grep(locus.plot,var.list)]))

if(is.null(x.lim)){x.lim <- range(sub.data)}

hist(sub.data, main=plot.title,xlab="Allele #",col="dodgerblue",xlim=x.lim,breaks=breaks.use) # simple plot, 




} # end genepop.plot.basicbar




#---------------------------------------------------------------------------
######################################################################################
# PLOT A PANEL WITH SCATTERPLOT OF ALLELE VALUES FOR 1 LOCUS AND 1 GROUP 
#  (e.g. plot "loc_1b_1" vs "loc_1b_2" , by using input "loc_1b")

genepop.plot.basicscatter <- function(X,locus.sub=1,group.sub=1,plot.title=NULL,x.lim=NULL,y.lim=NULL,
		add=FALSE, pch.use=1,col.use="dodgerblue",cex.use=1.3,jitter.do=TRUE){
# X is a list object created by genepop.read, possibly after modifying the GROUP variable for the X$mat object
# locus.sub = either a number or a label
# group.sub = either a number or a label
# add = if FALSE, create a new plot, if TRUE, add to existing plot
# jitter.do shuffles the points a bit to show density

if(is.numeric(group.sub)){group.plot <- unique(X$mat[,"GROUP"])[group.sub]}
if(!is.numeric(group.sub)){group.plot <- group.sub}

if(is.numeric(locus.sub)){locus.plot <- X$loci.list[group.sub]}
if(!is.numeric(locus.sub)){locus.plot <- locus.sub}

if(is.null(plot.title)){plot.title <- paste(group.plot,locus.plot,sep="")}

var.list <- dimnames(X$mat)[[2]]

print(group.plot)
print(locus.plot)

sub.data.x <- as.numeric(unlist(X$mat[grep(group.plot,X$mat[,"GROUP"]) ,  grep(paste(locus.plot,"_1",sep=""),var.list)]))
sub.data.y <- as.numeric(unlist(X$mat[grep(group.plot,X$mat[,"GROUP"]) ,  grep(paste(locus.plot,"_2",sep=""),var.list)]))

if(jitter.do){
	sub.data.x <- jitter(sub.data.x)
	sub.data.y <- jitter(sub.data.y)
	}

if(is.null(x.lim)){x.lim <- range(sub.data.x)}
if(is.null(y.lim)){y.lim <- range(sub.data.y)}

plot(sub.data.x,sub.data.y , main=plot.title,xlab="Lower Allele #",ylab="Higher Allele #",
		col=col.use,xlim=x.lim,ylim=y.lim,pch=pch.use,cex=cex.use) 


} # end genepop.plot.scatterplot

