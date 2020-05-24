###This Rscript is a specialised script that will plot the families per given phylum along, coloredby type of dataset ( whether it has been contaminatied by a sourc or not)
###The input file needs to have the mentioned elements: 
###sk.ph -> dataset
###"phylum"= Phylum, "n"= number of phyla per family (
sk.ph[!is.na(sk.ph[,"Phylum"]),"phylum"] <- sk.ph[!is.na(sk.ph[,"Phylum"]),"Phylum"]
##Last part was not useful but it was good to make sure.. like lmj also found out as well below
## table(sk.ph[sk.ph$lump==T,"n"])
#  0 
#351 

br.col <- c("goldenrod2","cyan4", "tomato2", "dodgerblue")

color.picker.raw <- function(x){                                                                                     
    if(x == -1 ){return(br.col[1])}                                                                            
            else if(x == 0 ){return(br.col[2])}
            else if(x == 1 ){return(br.col[3])}
            else {return(br.col[4])}                                                                                               
}

sk.ph$cols <- sapply(sk.ph$n, color.picker.raw)                                                                   

sk.ph$category <- cut(sk.ph$n, c(-Inf, -1, 0, 1,Inf),
                      c("Failed mp+lf", "Failed lf", "Passed mp+lf","Failed mp"))  
ticks <- function(x){
    if(x <= 10)
    {return(1)}
    else if(x>50)
    {return(10)}
    else{return(5)}
}
nrow(sk.ph)
## That makes a lot of sense, but if we want to make plots for the family composition of each
## family then we need to do something a little bit trickier..
##############
phylo.bar <- function(p, leg.h=6, leg=TRUE){
    ##where p stands for phylum and leg.h stands for legend height
    phyl.all<- sk.ph[grep(p, sk.ph$phylum),]
    #phyl.all<- sk.ph[grep(p, sk.ph$phylum),]
    # phyl.all<- sk.bl.2[sk.bl.2$phylum == p,]
    ##First see how many species are there per family
    fam.ph
    fam.ph<- aggregate(phyl.all$acc, by=list(phyl.all$family), length)
    colnames(fam.ph)<- c("family", "Freq")
    phyl.all<- merge(fam.ph, phyl.all, by=c("family"))
    phyl<- phyl.all[!duplicated(phyl.all$family),]
    phyl <- phyl[order(phyl$Freq, decreasing=FALSE),]
    
    plot(phyl$Freq, xlim=c(0,eval(max(phyl$Freq)+(0.1*max(phyl$Freq)))),
         ylim=c(0,nrow(phyl)+ 1) ,# yaxs="i",
         type="n", axes=F, ylab="", xlab="Number of alignments", main=unique(phyl$p))
    axis(1, seq(0,round(eval((max(phyl$Freq)+(0.1*max(phyl$Freq))))),
                by=eval(ticks(max(phyl$Freq)))),
         cex.axis=0.8, las=1)
    
    axis(2, 1:nrow(phyl), labels=phyl$family, las=2, cex.axis=1)
        ##So the square would be as long as all the families, as wide as the maximum score of alignments
    y<- 1:nrow(phyl)
                                        #for(i in 1:nrow(phyl)
    for(j in 1:(nrow(phyl[order(phyl$n, phyl$Freq, decreasing=FALSE),]))){
        fam<- phyl[j,"family"]
        sk.hom <- phyl.all[phyl.all$family== fam,]
        sk.hom$n[sk.hom$n >=2] <- 2
        fam.freq=as.data.frame(table(sk.hom$n))
        colnames(fam.freq) <- c("id","freq")
        phyl[phyl$family==fam,"Freq"]
        ##You dont need a percent dumbass..
        colindex <- data.frame("colors"=c("goldenrod2","cyan4","tomato2", "dodgerblue"),"id"=c(-1,0,1,2))
###if lmj asks me to putt a better order
### co<-c("B", "A" ,"C")           
### df$b<- factor( as.character(df$b), levels=co )
### dd <- df[order(df$b),]
    
        fam.freq<- merge(colindex, fam.freq, by=c("id"))
        fam.freq$colors<- as.character(fam.freq$colors)
        
        rect(0,j-0.4,(fam.freq$freq[1]),j+0.4, col=fam.freq$colors[1], border=fam.freq$colors[1])
        
        for(i in 1:(length(fam.freq$freq)-1)){
            rect(sum(fam.freq$freq[1:i]),
                 j-0.4,
                 (sum(fam.freq$freq[1:(i+1)])),
                 j+0.4,
                 col=(fam.freq$colors[i+1]),
                 border=(fam.freq$colors[i+1]))
        }
    }       
    box()
    if (leg==TRUE){    
    ##To create the legend, it has to be modular..
        xl=eval(max(phyl$Freq-(max(phyl$Freq)/2.3)))
        yb =0
        xr=eval(max(phyl$Freq)+(0.1*max(phyl$Freq)))
####I made the legend height a variable so that if I need to transform specific graphs, it doesnt squish or
####expand significantly
        yt=eval(nrow(phyl)/leg.h)
        
        rect(xleft=xl,
             ybottom= yb ,
             xright=xr,
             ytop=yt,
             col="white"
             )
        
        text(x=eval((xl+xr)/2), y=(yt+(0.25*yt)), labels="Alignment types", cex=1)
        
    ##Get a unique color for each group
        ##Edit, changed c input from phyl to phyl.all as the former would ignore more than one color on the same family
        c<- phyl.all[!duplicated(phyl.all$cols),]
        n<- length(unique(phyl$cols))
        ##Give dimensions of the content area of the legend
        xlplus<- xl+(0.07*(xr-xl))
        xrplus<- xr-(0.07*(xr-xl))
        ybplus<- yb+(0.07*(yt-yb))
        ytplus<- yt-(0.15*(yt-yb))
        xhplus<- ((xlplus+xrplus)/2)
        
        ##And fill the content area with the legend content
        for(i in 1:nrow(c)){
            y<- (ytplus-ybplus)/(nrow(c)+1)
        rect(xlplus*0.98, y*i, xhplus*0.97, y*(i+1), col = c$cols[order(c$n)[i]], border="white")
            
            text(x=eval((xrplus+xhplus)/1.98), y=eval((y*i+y*(i+1))/2),
                 labels=c$category[order(c$n)[i]], cex=0.78)
        }
    }
}
