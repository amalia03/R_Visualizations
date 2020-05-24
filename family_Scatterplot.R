###Let's now try to design some scatterplots, I will be taking the first parts of the function above and plotting alignments against the biomass and count

j.tax <- read.table(file="jelly_taxa.csv", stringsAsFactors=FALSE, header=TRUE, na.strings="N/A")

j.sp<- read.table("jfarm_special_phyla.csv", stringsAsFactors=FALSE, header=TRUE, na.strings="N/A")
j.tax <- rbind(j.sp,j.tax)

j.tax<- j.tax[!j.tax$family=="Uknown",]
j.tax[is.na(j.tax$number),"number"] <- 0


##Establish the color palette
pal2<- colorRampPalette(c("royalblue","violetred1","red1"))
col.pal<- data.frame("pct"= seq(0, 100 ,1),"col.index"=as.character(pal2(101)))
col.pal$col<- as.character(col.pal$col.index)
sk.ph.n<- sk.ph[sk.ph$n==1 & sk.ph$family!="Hominidae" & sk.ph$family!="Monodopsidaceae" & sk.ph$lump == FALSE,]

scatter.fam<- function(phyl){
    jtax<- j.tax[j.tax$phylum == phyl,]
#jtax<- j.tax[j.tax$phylum == "Annelida",]

j.n<- aggregate(jtax[,"number"], by=list(jtax$family), FUN=sum)
    colnames(j.n)<- c("family", "m_ct")
  
    j.b<- aggregate(jtax[,"biomass"], by=list(jtax$family), FUN=sum)
    colnames(j.b)<- c("family", "m_b")
    j.ph<- merge(j.b, j.n, by="family")
    j.ph$m_r<- ((j.ph$m_b*1000)/j.ph$m_ct)+1

sk.ph.n<- sk.ph.n[sk.ph.n$phylum ==phyl,]
#    sk.ph.n<- sk.ph.n[sk.ph.n$phylum =="Annelida",]
    sk.ph.n<- aggregate(sk.ph.n$acc, by=list(sk.ph.n$family), length)
    
    colnames(sk.ph.n)<- c("family", "eRNA")
    
    phyl.u<- unique(c(sk.ph.n$family, j.ph$family))
    
    for(i in 1:length(phyl.u)){
        if((phyl.u[i] %in% sk.ph.n$family)==FALSE){
            empty.fam <- c(phyl.u[i],0)
            sk.ph.n   <-  rbind(sk.ph.n, empty.fam)
        }   
    }
    
    for(i in 1:length(phyl.u)){
        if((phyl.u[i] %in% j.ph$family)==FALSE){
            empty.fam <- c(phyl.u[i],0,0,0)
            j.ph  <-  rbind(j.ph, empty.fam)
        }   
    }
    
    ph.n<-merge(sk.ph.n,j.ph, by="family")
    empty.fam <- c(phyl.u[i],0)
    ph.n$m_r<- as.numeric(ph.n$m_r) +1
    ph.n$eRNA<- as.numeric(ph.n$eRNA) +1
    ##To fix the values on the m_ct column:
#    ph.n <- ph.n[!ph.n$family=="Amphinomidae",]
        
####for the damn color..
    col.r<- function(x, perc.subj="ph"){
        if(perc.subj=="ph"){
            if(nrow(sk.ph[sk.ph$family==x,])>0){
                sk.hom <- sk.ph[sk.ph$family==x,]
                sk.hom$n[sk.hom$n >=2] <- 2
                fam.freq<- as.data.frame(table(sk.hom$n))
                colnames(fam.freq) <- c("type","freq")
                fam.freq$r <- (fam.freq$freq/sum(fam.freq$freq))*100
                
                if (fam.freq[1,"type"]==1){
                    col.b<-as.character(col.pal[col.pal$pct==round(fam.freq[1,"r"],0), "col"])
                    
                }else{
                    col.b<- as.character(col.pal[1,2])
                }        
                return(col.b)
            }
            else{
                return("darkgrey")
            }                                             
    }
        if(perc.subj=="bm"){
            fam.r<- as.numeric(ph.n[ph.n$family==x,"m_b"])/sum(as.numeric(ph.n$m_b))*100
            col.b<-as.character(col.pal[col.pal$pct==round(fam.r,0), "col"])
            return(col.b)
        }
        if(perc.subj=="ct"){
            fam.r<- as.numeric(ph.n[ph.n$family==x,"m_ct"])/sum(as.numeric(ph.n$m_ct))*100
            col.b<-as.character(col.pal[col.pal$pct==round(fam.r,0), "col"])
        }
    }

ph.n$col<-sapply(ph.n[,"family"], function(x) col.r(x, perc.subj="ct"))

    ph.n$log.r<- round(log(ph.n$m_r,10),1)
    ph.n$log.erna<- round(log(ph.n$eRNA,2),2)

    ##And now the plotting part
par(mar=c(5, 5, 3, 1))
    plot(ph.n$log.erna, round(ph.n$log.r,1),
         main= phyl,

         xlim=c(min(ph.n$log.erna)-((range(ph.n$log.erna)[2]-range(ph.n$log.erna)[1])*0.02),
                max(ph.n$log.erna)*1.1),
         ylim=c(min(log(ph.n$m_r, 10))-((range(log(ph.n$m_r),10)[2]-range(log(ph.n$m_r, 10))[1])*0.02),
         max(log(ph.n$m_r, 10)*1.1)),

         ##For forams:
      ##   ylim=c(min(ph.n$m_r)-((range(ph.n$m_r)[2]-range(ph.n$m_r)[1])*0.02),
       ##         max(ph.n$m_r)*1.1),
         xaxs="i",
         yaxs="i",
         xlab="Number of aligning genes (log2)",,
         ##For forams
         ylab="Ratio of of biomass (mg) to counts of morphological samples (log10)", pch=19,
         ##ylab="Percent of morphological biomass to counts ratio (mg, log10 +1)",pch=19,
         col=ph.n$col, cex=1.5
         )
    mx<- max(log(ph.n$eRNA,2))
    
###This function will correct for the position of the text so that it will always be visible within the graph
    pos.t<- function(x){
        if(x <= mx-mx/5)
        {return(4)}
        else
        {return(2)}
    }
    
###so here (21.4.2020) I want to do something for the overlapping text.
x.range<- range(ph.n$log.erna)[2]- range(ph.n$log.erna)[1]
y.mv<- as.data.frame(table(ph.n$log.r))

###but (8.05.2020) i need to fix the x coordinates too 
#Take1:
#x.mv<- unique(ph.n$log.erna); names(x.mv)<- x.mv
#x.dist<- outer(x.mv,x.mv, '-')
#x.dist<- abs(x.dist)
#Take2: 
#ph.y<- ph.y[order(ph.y$log.erna),]
#x.mv<- cbind(ph.y,dist=c(ph.y$log.erna[1],diff(ph.y$log.erna)))
                                        #a<- x.mv[x.mv$dist!=0 & x.mv$dist> x.range/6.5,]
                                        #b<- x.mv[!(a[1,]%in%x.mv[1,]),]
    
    max.log<- max(ph.n$log.r)
    
    ##if(nrow((y.mv[y.mv$Freq >1,])) > 0){
    n.s<- ph.n[ph.n$log.r %in% y.mv[y.mv$Freq>1,"Var1"],]
    n.s<- unique(n.s$log.r)
    
    for(j in 1:length(n.s)){ 
            high.y<- ph.n[ph.n$log.r==n.s[j],]
            high.y<- high.y[order(high.y$log.erna),]
            high.y<- cbind(high.y,dist=c(0,diff(high.y$log.erna)))
                                        # high.y<- cbind(high.y,dist=c(high.y$log.erna[1],diff(high.y$log.erna)))
            
            b <- high.y[high.y$dist==0 | high.y$dist < x.range/6.5,]
            
            if(nrow(b)>1){
                b$new.y<- seq(b$log.r[1], b$log.r[1]+max.log*0.1,
                          by=((b$log.r[1] + max.log*0.1)-b$log.r[1])/(nrow(b)-1))

            for( i in 1:nrow(b)){
                polygon(c(b$log.erna[i],b$log.erna[i]+0.25), c(b$log.r[i], b$new.y[i]), border="grey")
            }
            text(b$log.erna+0.2, b$new.y, labels=b$family, pos=sapply(b$log.erna, pos.t), cex=0.8) 
            }else{
                text(b$log.erna, b$log.r, labels=b$family, pos=sapply(b$log.erna, pos.t ), cex=0.8)
           }
            if(nrow(high.y[high.y$dist!=0 & high.y$dist > x.range/6.5,]) !=0)
            {
            dist.y<- high.y[high.y$dist!=0 & high.y$dist > x.range/6.5,]
            text(dist.y$log.erna, dist.y$log.r, labels=dist.y$family, pos=sapply(dist.y$log.erna, pos.t ), cex=0.8)
            }
        }
    
        c<- ph.n[!ph.n$log.r %in% y.mv[y.mv$Freq>1,"Var1"],]
        
    text(c$log.erna, c$log.r, labels=c$family, pos=sapply(c$log.erna, pos.t ), cex=0.78)

###For forams:
     ## text(ph.n$eRNA, ph.n$m_r, labels=ph.n$family, pos=sapply(ph.n$eRNA, pos.t ), cex=0.8)

#}
    
    
    ##For the spectrum legend (20.04.2029)
       
    plot(0:1, 0:1, ylim=c(0,150), 
         type="n", axes=F, ylab="", xlab="", xaxs="i")#, main="Percent of family passing second-phase filtering")
    
    rect(0.1,       
         seq(5, nrow(col.pal), by=10),
         0.5,
         seq(15, nrow(col.pal)+5, by=10),
         col=col.pal[seq(1, nrow(col.pal)-10, by=10),"col"],
         border=col.pal[seq(1, nrow(col.pal)-10, by=10),"col"]
         )
   
    
    box(col="grey")
    
    seq(0, nrow(col.pal), by= nrow(col.pal)/2)
    seq(5,nrow(col.pal)+5,by=50)
    text(x=0.8, y=seq(5, nrow(col.pal)+5, by=50), labels=paste(seq(0, 100, by=50), "%"), cex= 0.8)
#    text(x=0.5, y=135, labels="Genes passing filtering steps", cex=0.8)
    text(x=0.5, y=135, labels="Relative counts (%)", cex=0.8)
    
}


pdf("morpho_vs_erna_mctcol.pdf", width=8.37, height=11.67)


layout(matrix(c(1,0,0,
                1,2,0,
                0,0,0), ncol=3, byrow=T), widths=c(75,30,20), heights=c(90,50,165))
sapply("Annelida", function(x) scatter.fam(x))
dev.off()

layout(matrix(c(1,0,1,2,0,0), ncol=2, byrow=T), widths=c(70,25), heights=c(80,60,160))
sapply("Mollusca", function(x) scatter.fam(x))
layout(matrix(c(1,0,1,2,0,0), ncol=2, byrow=T), widths=c(70,25), heights=c(80,60,160))
sapply("Echinodermata", function(x) scatter.fam(x))
#sapply("Arthropoda", function(x) scatter.fam(x))

dev.off()
