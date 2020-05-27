mtz<- uniq.ph[c(2,3,6,14,16,17,21,31,32,39,40)]
unique(sk.ph[sk.ph$phylum %in%  mtz, "phylum"])

sk.f.mtz.agg<- aggregate(pass[pass$phylum %in%  mtz, "acc"], by=list(pass[pass$phylum %in%  mtz, "phylum"]), length)
macro.agg<-aggregate(j.tax[j.tax$phylum %in%  mtz, "number"], by=list(j.tax[j.tax$phylum %in%  mtz, "phylum"]), sum)
sk.mtz.agg<- aggregate(sk.ph[sk.ph$phylum %in%  mtz, "acc"], by=list(sk.ph[sk.ph$phylum %in%  mtz, "phylum"]), length)


agg.plot<- function(taxa, aggroup="acc", bygroup="phylum", fun.type=length, inv=FALSE, mt=""){
###agg.plot<- function(taxa, aggroup="acc", bygroup="phylum", fun.type=length, inv=FALSE, mt=""){
#    taxa<- taxa[!duplicated(taxa$aggroup),]
    if (inv==FALSE){
        x<- aggregate(taxa[taxa$phylum %in%  mtz, aggroup], by=list(taxa[taxa$phylum %in%  mtz, bygroup]), fun.type)
    }else{
        x<- aggregate(taxa[!taxa$phylum %in%  mtz, aggroup], by=list(taxa[!taxa$phylum %in%  mtz, bygroup]), fun.type)
    }
    
    colnames(x)<- c("phylum","freq")
    
    x<- x[order(x$freq, decreasing=TRUE),]
          
    if(nrow(x)%%2==0){
    y<-rep(c("white","steelblue3"), round((nrow(x)/2),0))
    }else{
        y<- rep(c("white","steelblue3"), round((nrow(x)+1)/2,0))
        y <- y[-length(y)]
    }
    
    x <- cbind("col"= y,x)
    x<- cbind("pct"=(x$freq/sum(x$freq))*100, x)
    
    sm.pct<- x[x$pct<3,]
    sm.pct<- sm.pct[order(sm.pct$freq),]
    l.pct<- x[x$pct>3,]
    
####plotting the ratiocolumn
plot(x$pct, ylim=c(0,110),
     xlim=c(0,180),
     type="n", axes=F, ylab="", xlab="",
     xaxs="i", yaxs="i"
     )
    
    axis(2, seq(0,100,
                by=5),cex.axis=.9, las=1, col="darkgrey")
    
    rect(0,0,70,(x$pct[1]), col=as.character(x$col[1]), border="black")
###    text(x=0.5, y=x$pct[1]-1, labels=x$phylum[1], cex=0.85, adj=-0.1, srt=90, col="white", las=2)
    text(x=35, y=x$pct[1]/2, labels=x$phylum[1], cex=0.8, col="black")

    for(i in 1:(length(x$pct)-1)){
        rect(0,
             sum(x$pct[1:i]),
             70,
             (sum(x$pct[1:(i+1)])),
             col=as.character(x$col[i+1]),
             ##             col=rep(c("white","steelblue3"), round((nrow(x)/2)-1,0)+1),
             ##   border=as.character(x$col[i+1])
             border="black"
             )    
#        text(x=0.5,y=sum(x$pct[1:(i+1)])-1, labels=x$phylum[i+1], cex=0.85, adj=-0.1, srt=90, col="white")
        
        if(x$pct[i]>3){
        #    if(x$pct[i] %in% sm.pct){
            ##            text(x=40,y=sum(x$pct[1:(i+1)])-0.7, labels=x$phylum[i+1], cex=0.85, col="black")
                                        #text(x=40,y=sum(x$pct[1:(i+1)])-1, labels=x$phylum[i+1], cex=0.80, col="black")
                text(x=35,y=(sum(x$pct[1:i])+sum(x$pct[1:(i+1)]))/2, labels=x$phylum[i+1], cex=0.75, col="black")
                ##  text(x=sum(x$pct[1:(i+1)])-0.5, y=0, labels=x[x$pct> 3, "phylum"][i+1], cex=0.9, adj=-0.1, srt=90)
         #   }       
        }
}    
    base<- (100+sum(l.pct$pct))/2
    
    polygon(c(75,83),c(sum(l.pct$pct),(sum(l.pct$pct))),border="darkgrey")
    polygon(c(75,83),c(100,100),border="darkgrey")
    polygon(c(83,83),c(sum(l.pct$pct),100),border="darkgrey")
    polygon(c(83,88),c(base,base),border="darkgrey")
    
    
    for (i in nrow(sm.pct):1){
        text(x=130,y=base-(i*2.2-2), labels=paste(sm.pct$phylum[i],",",round(sm.pct$pct[i],1),"%"), cex= 0.70)
    }
###And finally for the main title    
    text(x=65,y=105,labels=paste(mt,"\n","n=",sum(x$freq)), cex=0.85)

}   
