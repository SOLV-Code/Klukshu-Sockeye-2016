
boxplot.comparison.v1 <- function(X, var.plot="ForkLength_mm",var.label="Fork Length (mm)" ,
					groupings.set="Set1",x.lim=NULL){

if(is.null(x.lim)){ x.ticks <- pretty(range(X[,var.plot],na.rm=TRUE))}
if(!is.null(x.lim)){ x.ticks <- pretty(range(x.lim))}
numobs.x <- 1 + max(x.ticks)*1.04

plot(1:10,1:10,bty="n",xlim=range(x.ticks), type="n",ylim=c(0,10),axes=FALSE,xlab=var.label, ylab="")

text(numobs.x,10.5, label="# Obs",xpd=NA,font=2,adj=0)

# all samples
boxplot(X[,var.plot],horizontal=TRUE,add=TRUE,at=9.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
text(numobs.x,9.5, label=length(na.omit(X[,var.plot])),xpd=NA,adj=0)

if(groupings.set %in% c("Set1","Set2")){

male.idx <- X[,"Gender"]=="M"
female.idx <- X[,"Gender"]=="F"

boxplot(X[female.idx,var.plot],horizontal=TRUE,add=TRUE,at=8,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[male.idx,var.plot],horizontal=TRUE,add=TRUE,at=7.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
text(numobs.x,8, label=length(na.omit(X[female.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,7.5, label=length(na.omit(X[male.idx,var.plot])),xpd=NA,adj=0)

}

if(groupings.set=="Set1"){

river.idx <- X[,"DNA_Group"]=="River"
lake.idx <- X[,"DNA_Group"]=="Lake"

boxplot(X[river.idx & female.idx,var.plot],horizontal=TRUE,add=TRUE,at=6,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[lake.idx & female.idx,var.plot],horizontal=TRUE,add=TRUE,at=5.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
text(numobs.x,6, label=length(na.omit(X[river.idx & female.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,5.5, label=length(na.omit(X[lake.idx &female.idx,var.plot])),xpd=NA,adj=0)


#segments(median(X[river.idx & female.idx,var.plot],na.rm=TRUE),6, 
#		median(X[lake.idx & female.idx,var.plot],na.rm=TRUE),5.5,
#		col="red", lwd=2)

boxplot(X[river.idx & male.idx,var.plot],horizontal=TRUE,add=TRUE,at=4.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[lake.idx & male.idx,var.plot],horizontal=TRUE,add=TRUE,at=4,bty=NULL,border="darkblue",col="lightgray",lty=1 )
text(numobs.x,4.5, label=length(na.omit(X[river.idx & male.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,4, label=length(na.omit(X[lake.idx & male.idx,var.plot])),xpd=NA,adj=0)

#segments(median(X[river.idx & male.idx,var.plot],na.rm=TRUE),4.5, 
#		median(X[lake.idx & male.idx,var.plot],na.rm=TRUE),4,
#		col="red", lwd=2)

}


if(groupings.set=="Set2"){

rivertag.idx <- X[,"Tag_Fate"]=="River"
laketag.idx <- X[,"Tag_Fate"]=="Lake"

boxplot(X[rivertag.idx,var.plot],horizontal=TRUE,add=TRUE,at=6,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[laketag.idx,var.plot],horizontal=TRUE,add=TRUE,at=5.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
text(numobs.x,6, label=length(na.omit(X[rivertag.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,5.5, label=length(na.omit(X[laketag.idx,var.plot])),xpd=NA,adj=0)


rivergsi.idx <- X[,"RB_Match1"] %in% c("AdSpn_KlukshuRiver2016","Neskataheen")
lakegsi.idx <- X[,"RB_Match1"]=="AdSpn_KlukshuLake2016"

boxplot(X[rivergsi.idx,var.plot],horizontal=TRUE,add=TRUE,at=4.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[lakegsi.idx ,var.plot],horizontal=TRUE,add=TRUE,at=4,bty=NULL,border="darkblue",col="lightgray",lty=1 )
text(numobs.x,4.5, label=length(na.omit(X[rivergsi.idx ,var.plot])),xpd=NA,adj=0)
text(numobs.x,4, label=length(na.omit(X[lakegsi.idx,var.plot])),xpd=NA,adj=0)

}


if(groupings.set %in% c("Set1","Set2")){

weir.idx <- X[,"Location"] == "Weir"
early.idx <- X[,"StatWeek"] < 34
mix.idx <- X[,"StatWeek"] == 34
late.idx <- X[,"StatWeek"] > 34

boxplot(X[male.idx & weir.idx & early.idx ,var.plot],horizontal=TRUE,add=TRUE,at=3,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[male.idx & weir.idx & mix.idx ,var.plot],horizontal=TRUE,add=TRUE,at=2.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[male.idx & weir.idx & late.idx ,var.plot],horizontal=TRUE,add=TRUE,at=2,bty=NULL,border="darkblue",col="lightgray",lty=1 )
text(numobs.x,3, label=length(na.omit(X[male.idx & weir.idx & early.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,2.5, label=length(na.omit(X[male.idx & weir.idx & mix.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,2, label=length(na.omit(X[male.idx & weir.idx & late.idx,var.plot])),xpd=NA,adj=0)

boxplot(X[female.idx & weir.idx & early.idx ,var.plot],horizontal=TRUE,add=TRUE,at=1,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[female.idx & weir.idx & mix.idx ,var.plot],horizontal=TRUE,add=TRUE,at=.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[female.idx & weir.idx & late.idx ,var.plot],horizontal=TRUE,add=TRUE,at=0,bty=NULL,border="darkblue",col="lightgray",lty=1 )
text(numobs.x,1, label=length(na.omit(X[female.idx & weir.idx & early.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,.5, label=length(na.omit(X[female.idx & weir.idx & mix.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,0, label=length(na.omit(X[female.idx & weir.idx & late.idx,var.plot])),xpd=NA,adj=0)

}




if(groupings.set=="Set1"){
axis(2,at=c(9.5,8,7.5,6,5.5,4.5,4,3,2.5,2,1,0.5,0), las=2,
	labels=c("All Samples","All Females", "All Males", "Females-River Spn", "Females-Lake Spn", "Males-River Spn", "Males-Lake Spn", 
			"M-Weir-W28toW33","M-Weir-W34","M-Weir-W35toW41","F-Weir-W28toW33","F-Weir-W34","F-Weir-W35toW41" ))
}
			
if(groupings.set=="Set2"){
axis(2,at=c(9.5,8,7.5,6,5.5,4.5,4,3,2.5,2,1,0.5,0), las=2,
	labels=c("All Samples","All Females", "All Males", "Tag-River", "Tag-Lake", "GSI-River", "GSI-Lake", 
			"M-Weir-W28toW33","M-Weir-W34","M-Weir-W35toW41","F-Weir-W28toW33","F-Weir-W34","F-Weir-W35toW41" ))
}			
			
			
			
			
			
			
			
			
			
			
if(groupings.set=="Set3"){



male.idx <- X[,"Gender"]=="M"
female.idx <- X[,"Gender"]=="F"

rivertag.idx <- X[,"Tag_Fate"]=="River"
laketag.idx <- X[,"Tag_Fate"]=="Lake"

rivergsi.idx <- X[,"RB_Match1"] %in% c("AdSpn_KlukshuRiver2016","Neskataheen")
lakegsi.idx <- X[,"RB_Match1"]=="AdSpn_KlukshuLake2016"


boxplot(X[male.idx  & rivertag.idx,var.plot],horizontal=TRUE,add=TRUE,at=8,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[male.idx  & rivergsi.idx,var.plot],horizontal=TRUE,add=TRUE,at=7.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[male.idx  & laketag.idx,var.plot],horizontal=TRUE,add=TRUE,at=7,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[male.idx  & lakegsi.idx,var.plot],horizontal=TRUE,add=TRUE,at=6.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )

text(numobs.x,8, label=length(na.omit(X[male.idx  & rivertag.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,7.5, label=length(na.omit(X[male.idx  & rivergsi.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,7, label=length(na.omit(X[male.idx  & laketag.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,6.5, label=length(na.omit(X[male.idx  & lakegsi.idx,var.plot])),xpd=NA,adj=0)

boxplot(X[female.idx  & rivertag.idx,var.plot],horizontal=TRUE,add=TRUE,at=5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[female.idx  & rivergsi.idx,var.plot],horizontal=TRUE,add=TRUE,at=4.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[female.idx  & laketag.idx,var.plot],horizontal=TRUE,add=TRUE,at=4,bty=NULL,border="darkblue",col="lightgray",lty=1 )
boxplot(X[female.idx  & lakegsi.idx,var.plot],horizontal=TRUE,add=TRUE,at=3.5,bty=NULL,border="darkblue",col="lightgray",lty=1 )

text(numobs.x,5, label=length(na.omit(X[female.idx  & rivertag.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,4.5, label=length(na.omit(X[female.idx  & rivergsi.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,4, label=length(na.omit(X[female.idx  & laketag.idx,var.plot])),xpd=NA,adj=0)
text(numobs.x,3.5, label=length(na.omit(X[female.idx  & lakegsi.idx,var.plot])),xpd=NA,adj=0)

axis(2,at=c(9.5,8,7.5,7,6.5, 5,4.5,4,3.5), las=2,
	labels=c("All Samples","M-Tag-River", "M-GSI-River","M-Tag-Lake","M-GSI-Lake","F-Tag-River", "F-GSI-River","F-Tag-Lake","F-GSI-Lake"))
		
			
			
			
			
}

			
}
