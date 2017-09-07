
#function to fit and plot phylogenetic trees

treefit.fun <- function(file.in,dist.method=2,tree.package="ape",tree.method = "nj", plot.type="phylogram"){
# Note: applies functions from {ape} , {adegenet} , {phangorn}
# for genetic distance method see:  dist.genpop() in {adegenet)
# for plot type see plot.phylo() in {ape}
# Note: taps into different tree fitting packages - for now choose between "ape" and "phangorn". phangorn has a bootstrap approach built in
# tree.method: 
#      - if using {ape}, then either "nj" or "bionj" for now - both of these are functions in {ape}, see documentation there
#      - if using {phangorn}, then either "NJ", "UNJ", "upgma", "wpgma" for now - both of these are functions in {ape}, see documentation there


 require("adegenet")

test.obj <- read.genepop(file=file.in,ncode = 2L, quiet = FALSE)# read a genepop file (creates object of class genind)
test.obj.genpop <- genind2genpop(test.obj)# convert to object of class "genpop"
dist.mat <- dist.genpop(test.obj.genpop ,method=dist.method)   ## calculate genetic distances(method=1 gives Nei's distance, method= 2 gives Edward's distance, method=3 gives Cavalli Sforza)


if(tree.package=="ape"){
	require("ape")
	if(tree.method != "nj" & tree.method != "bionj"){warning("Specified tree.method not available for the specified tree.package!");stop()}
	# fit a tree (creates object of class phylo)
	if(tree.method=="nj"){       # fit neighbour-joining tree
		tree.fit <- nj(dist.mat) 
		plot(tree.fit,type=plot.type,cex=1.3)
		#tiplabels(tree.fit$tip,font=2,family="mono",frame="n",cex=1.5)

		
		
		} 


	if(tree.method=="bionj"){           # fit improved neighbour-joining tree
		tree.fit <- bionj(dist.mat) 
		plot(tree.fit,type=plot.type,,cex=1.3)
		#tiplabels(tree.fit$tip,font=2,family="mono",frame="n",cex=1.5)


		} 





	}


if(tree.package=="phangorn"){
	require("phangorn")
	if(!(tree.method  %in% c("NJ","UNJ","upgma","wpgma"))){warning("Specified tree.method not available for the specified tree.package!");stop()}
	
	
	# fit a tree (creates object of class phylo)
	if(tree.method=="NJ"){tree.fit <- NJ(dist.mat) } # fit neighbour-joining tree
	if(tree.method=="UNJ"){tree.fit <- bionj(dist.mat) } # fit unweighted neighbour-joining tree
	if(tree.method=="upgma"){tree.fit <- upgma(dist.mat) } # fit unweighted clustering
	if(tree.method=="wpgma"){tree.fit <- wpgma(dist.mat) } # fit weighted clustering ?
	plot(tree.fit,type=plot.type)  #,tip.label=FALSE)
	#tiplabels(tree.fit$tip,font=2,family="mono",frame="n",cex=1.5)
	
	# parsimony tests not yet implemented
	#pars.out <- parsimony(tree.fit,dist.mat)
	#print(names(pars.out))
	#print(summary(pars.out))

	}

return(list(tree=tree.fit, distances=dist.mat))


}
