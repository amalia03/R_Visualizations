##12.04.2020

j.tax <- read.table(file="jelly_taxa.csv", stringsAsFactors=FALSE, header=TRUE, na.strings="N/A")
##28.04.2020
##Adding the stra
plot.phylum<- function(phyl, type, title){
j.tax<- j.tax[j.tax$phylum ==phyl,]
    
#    j.tax<- j.tax[j.taxphylum =="Annelida",]
    j.tax<- j.tax[!j.tax$family=="Uknown",]
    j.tax[is.na(j.tax$number),"number"] <- 0
    

    j.ph<- aggregate(j.tax[,type], by=list(j.tax$family), FUN=sum)
    #j.ph<- aggregate(j.tax[,"number"], by=list(j.tax$family), FUN=sum)
    colnames(j.ph)<- c("family", "count")
    j.ph
    
    j.ph$r<- (j.ph$count/sum(j.ph$count)*100)
    ##barplot(j.ph$r, names.arg=j.ph$family)
    j.ph$group<- paste("Morpho",type)
    j.ph$colors<- "dodgerblue"
  
    ##sk.ph.n<-sk.ph
    sk.ph.n<- sk.ph[sk.ph$n==1 & sk.ph$family!="Hominidae" & sk.ph$family!="Monodo
psidaceae" & sk.ph$lump == FALSE,]
    sk.ph.n<- sk.ph.n[sk.ph.n$phylum == phyl,]
   # sk.ph.n<- sk.ph.n[sk.ph.n$phylum == "Annelida",]

    sk.ph.n<- aggregate(sk.ph.n$acc, by=list(sk.ph.n$family), length)
   
    colnames(sk.ph.n)<- c("family", "count")
    sk.ph.n$r<- (sk.ph.n$count/sum(sk.ph.n$count)*100)
    ##barplot(sk.ph.n$r, names.arg=sk.ph.n$family)
    sk.ph.n$group<- "Number of contigs"
    sk.ph.n$colors<- "tomato2"
    
    phyl.u<- unique(c(sk.ph.n$family, j.tax$family))
    
    for(i in 1:length(phyl.u)){
        if((phyl.u[i] %in% sk.ph.n$family)==FALSE){
            empty.fam <- c(phyl.u[i],0,0,"Number of contigs","tomato2")
            sk.ph.n   <-  rbind(sk.ph.n, empty.fam)
        }   
    }

    for(i in 1:length(phyl.u)){
        if((phyl.u[i] %in% j.ph$family)==FALSE){
            empty.fam <- c(phyl.u[i],0,0,paste("Morpho",type),"dodgerblue")
            #empty.fam <- c(phyl.u[i],0,0,paste("Morpho","number"),"dodgerblue")
            j.ph  <-  rbind(j.ph, empty.fam)
        }   
    }
    
sk.ph.n$r<- as.numeric(sk.ph.n$r)
    j.ph$r<- as.numeric(j.ph$r)
    
    ph.n<-rbind(sk.ph.n,j.ph)

    phyl<- ph.n
    phyl<- phyl[order(phyl$family, decreasing=TRUE),]

#####
    ##for the actual plotting section: 
    par(mar=c(5, 8, 2, 2))
    plot(phyl$r, xlim=c(0,1.1*max(phyl$r)),
         ylim=c(0,nrow(phyl) +1 ), xaxs="i",yaxs="i",
         type="n", axes=F, ylab="", xlab="Family ratio", main=title)
    
    axis(1, seq(0,100, 5),
         cex.axis=0.8, las=1)
    
    axis(2, seq(1, nrow(phyl),2), labels=unique(phyl$family), las=2, cex.axis=1)
    
    
    seq(1,nrow(phyl),4)
    ##For zebra plotting
    for(i in seq(2,nrow(phyl),4)){
                                        #     for(i in 1:nrow(phyl)){
        rect(0,
             i-2,
             max(phyl$r),
             i, col="lightgray", border="lightgray"
             )
    }
    
    for(j in 1:eval(nrow(phyl))){
        rect(0,
        (j-0.8),
        phyl$r[j],
        (j-0.2),
        col=phyl$colors[j],
###         border=phyl$colors[j]
        border="black"
        )
    }
    
    
    ##ifyou want to add a line for the base of the barplots    
    rect(0,0.2,0,nrow(phyl)-0.2, border="black")
    
###to make the differnces between the lines greater I am using abline. 
    ##for(i in seq(1,nrow(phyl),2)){
    ##   abline(h=i-1, lty=4, col="darkgray")
    ##}
    ##abline(h=nrow(phyl), lty=4, col="lightgray")

    box()
    
    ##To create the legend, it has to be modular..
    xl=eval(max(phyl$r-(max(phyl$r))/2.3))
    yb =0
    xr=eval(max(phyl$r)+(0.1*max(phyl$r)))
####I made the legend height a variable so that if I need to transform specific graphs, it doesnt squish or
####expand significantly
    yt=eval(nrow(phyl)/6)
    
    c(xl,yb,xr,yt)
    rect(xleft=xl,
         ybottom= yb ,
         xright=xr,
         ytop=yt,
         col="white"
         )
    
text(x=eval((xl+xr)/2), y=(yt+(0.25*yt)), labels="Data type", cex=1)
    
    ##Get a unique color for each group
    ##Edit, changed c input from phyl to phyl.all as the former would ignore more than one color on the same family
    c<- phyl[!duplicated(phyl$colors),]
    
    ##Give dimensions of the content area of the legend
    xlplus<- xl+(0.07*(xr-xl))
    xrplus<- xr-(0.07*(xr-xl))
    ybplus<- yb+(0.07*(yt-yb))
    ytplus<- yt-(0.15*(yt-yb))
    xhplus<- ((xlplus+xrplus)/2)
    
    ##And fill the content area with the legend content
    for(i in 1:nrow(c)){
        y<- (ytplus-ybplus)/(nrow(c)+1)
        rect(xlplus*0.98, y*i, xhplus*0.97, y*(i+1), col = c$colors[order(c$r)[i]], border="white")
        text(x=eval((xrplus+xhplus)/1.98), y=eval((y*i+y*(i+1))/2),
             labels=c$group[order(c$r)[i]], cex=0.8)
    }
}

#pdf("j_farm_biased_rna.pdf", width=8.37, height=11.67)
pdf("rm_rna.pdf", width=8.37, height=11.67)
par(mfrow=c(2,1))
sapply("Annelida", function(x) plot.phylum(x, "number", paste("Annelida (Alignments vs Specimen Number)")))
sapply("Annelida", function(x) plot.phylum(x, "biomass", paste("Annelida (Alignments vs Specimen Biomass)")))
sapply("Mollusca", function(x) plot.phylum(x, "number", paste("Mollusca (Alignments vs Specimen Number)")))
sapply("Mollusca", function(x) plot.phylum(x, "biomass", paste("Mollusca (Alignments vs Specimen Biomass)")))
sapply("Echinodermata", function(x) plot.phylum(x, "number", paste("Echinodermata (Alignments vs Specimen Number)")))
sapply("Echinodermata", function(x) plot.phylum(x, "biomass", paste("Echinodermata (Alignments vs Specimen Biomass)")))
#l.2p.m(d3.w = 36, d3.h=45, p1.h=1.5, p2.h=1.5, s=1)
#sapply("Annelida", function(x) plot.phylum(x))
dev.off()
