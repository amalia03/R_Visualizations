agg.plot<- function(dataset, phyla, aggroup="acc", bygroup="phylum", fun.type=length, inv=FALSE, mt=""){
     
   ###Create a dataset with the aggregate function that contains the number of alignments per phylum 
    x<- aggregate(dataset[dataset$phylum %in%  phyla, aggroup], by=list(dataset[dataset$phylum %in% phyla, bygroup]), fun.type)
    ###Name the table columns
    colnames(x)<- c("phylum","freq")
    ###Order them by number of alignments per phylum.
    x<- x[order(x$freq, decreasing=TRUE),]
    ###Make the color of each bracket to interchange between grey and white.
    if(nrow(x)%%2==0){
    y<-rep(c("white","grey"), round((nrow(x)/2),0))
    }else{
        y<- rep(c("white","grey"), round((nrow(x)+1)/2,0))
        y <- y[-length(y)]
    }
    #Bind the color column to table x that has the alignment number information.
    x <- cbind("col"= y,x)
   
    #Create a column that contains the percentages of alignment number per phylum.
    x<- cbind("pct"=(x$freq/sum(x$freq))*100, x)
    
    #Separate phyla by small (not readable) and large percentages in two subsets called "sm.pct" and "l.pct".
    sm.pct<- x[x$pct<3,]
    sm.pct<- sm.pct[order(sm.pct$freq),]
    l.pct<- x[x$pct>3,]
    
    ##Making the plotting space
    plot(x$pct, ylim=c(0,110),
     xlim=c(0,160),
     type="n", axes=F, ylab="", xlab="",
     xaxs="i", yaxs="i"
     )
    
    ### Add an axis containing the 0-100 scale.
    axis(2, seq(0,100,
                by=5),cex.axis=.9, las=1, col="darkgrey")
    ###Add the first (most abundant phylum).
    rect(0,0,70,(x$pct[1]), col=as.character(x$col[1]), border="black")
    ###Add the text in the first rectangle. 
    text(x=35, y=x$pct[1]/2, labels=x$phylum[1], cex=0.83, col="black")

    ###Add the rest of the phyla who's height is equal to their percentage, sorted by most alignments to the bottom and less alignments to the top. 
    
    for(i in 1:(length(x$pct)-1)){
        rect(0,
             sum(x$pct[1:i]),
             70,
             (sum(x$pct[1:(i+1)])),
             col=as.character(x$col[i+1]),
             border="black"
             )            
        
        ###Add text to the larger percentage entries: 
        if(x$phylum[i] %in% l.pct$phylum){
                text(x=35,y=(sum(x$pct[1:i])+sum(x$pct[1:(i+1)]))/2, labels=x$phylum[i+1], cex=0.83, col="black")
        }
}    
    base<- (100+sum(l.pct$pct))/2
    
    #Create polygonal shapes that will represent the bracket and the topology where the text will be placed for the smaller percentages.
    #The first polygon creates a vertical line the height of all the small percentages.
    polygon(c(75,83),c(sum(l.pct$pct),(sum(l.pct$pct))),border="darkgrey")
    #The second polygon creates a horizontal line top of the plot.
    polygon(c(75,83),c(100,100),border="darkgrey")
    ##The third polygon creates a line on bottom of the first small percentage entry.
    polygon(c(83,83),c(sum(l.pct$pct),100),border="darkgrey")
    polygon(c(83,88),c(base,base),border="darkgrey")
    
    ###Add the text on the right hand size of all the phyla represented by smaller percentages in the correct order. 
    for (i in nrow(sm.pct):1){
        text(x=120,y=base-(i*3-3.4), labels=paste(sm.pct$phylum[i],",",round(sm.pct$pct[i],1),"%"), cex= 0.83)
    }
###And finally for the main title, it will have the main title as it has been filled in the main command, and the number of entries.     
    text(x=72,y=105,labels=paste(mt,"\n","n=",sum(x$freq)), cex=0.88)

}   
