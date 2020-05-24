###This script will plot barplots of two different groups that share the same factors, in this case families.  
###The script is a bit specialized for the type of datasets we have, one is a molecular dataset that has been finalized after multiple filtrings, 
###The other dataset has the morphological metrics such as biomass and number.

plot.phylum<- function(phyl, type, title){##where phyl is the phylum group in question, type is the type of numbers from the morphological dataset that are chosen ( eg biomass/count ratio), and title is the title of the plot    
###j.tax is the morphological dataset. Important entries include site (sampling site), "taxon" = the scientific name, "number", the number of individual finds, "biomass" and different taxonomic levels (from genus to phylum) 
###Keep the entries that match to the same phylum
    j.tax<- j.tax[j.tax$phylum ==phyl,]
###Remove all sequences that match to the string "uknown"
    j.tax<- j.tax[!j.tax$family=="Uknown",]
###Replace all entries that have N/A values in their numbers with 0      
    j.tax[is.na(j.tax$number),"number"] <- 0
    
###Create a sum per group table where the group is the one you designated in the command (total biomass/number ect per family)     
    j.ph<- aggregate(j.tax[,type], by=list(j.tax$family), FUN=sum)
    colnames(j.ph)<- c("family", "count")

###Create a percent ratio of the the frequency of the summed group (inexplicably named count for now)
    j.ph$r<- (j.ph$count/sum(j.ph$count)*100)
###Create a color column and a dataset column that only has the value "Moprho" that will serve as a label for the entries later
    j.ph$group<- paste("Morpho",type)
    j.ph$colors<- "dodgerblue"

###For the molecular dataset, keep only the values that dont have any contaminants ( dont have any human, nannochloropsis or lumpfish positives while they also do not match more than one family per accession)
    sk.ph.n<- sk.ph[sk.ph$n==1 & sk.ph$family!="Hominidae" & sk.ph$family!="Monodopsidaceae" & sk.ph$lump == FALSE,]
###Keep only the entries that match the phylum of interest
    sk.ph.n<- sk.ph.n[sk.ph.n$phylum == phyl,]

###Create a table That aggregaytes the acc column per family(how many alignments we have per family)
    sk.ph.n<- aggregate(sk.ph.n$acc, by=list(sk.ph.n$family), length)
    colnames(sk.ph.n)<- c("family", "count")
###Create a percent ratio of the the frequency of the summed group 
    sk.ph.n$r<- (sk.ph.n$count/sum(sk.ph.n$count)*100)
###Create a color column and a dataset column that only has the value "Number of contigs" that will serve as a label for the entries later
    sk.ph.n$group<- "Number of contigs"
    sk.ph.n$colors<- "tomato2"
    
###For the next section, we need to create entries for the absent values for each family per dataset so that both dataset have an equal number of families
###First we create a variable with all the unique family values 
    phyl.u<- unique(c(sk.ph.n$family, j.tax$family))

###For each dataset we ask whether the entries in the phyl.u variable match with each family name in the dataset. If it returns false, it creates an entry for that family that has 0 values in the numeric categories and rbinds them to the rest of the dataset
    
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

####Just to be sure, we ask for the ratio to be numeric 
sk.ph.n$r<- as.numeric(sk.ph.n$r)
j.ph$r<- as.numeric(j.ph$r)
###We compile the morphological and molecular datasets together in one called "ph.n"    
    ph.n<-rbind(sk.ph.n,j.ph)
###We make a copy of ph.n so that we can modify it but keep the unchanged one for later
    phyl<- ph.n
###Order the entries alphabetically by the family groups 
    phyl<- phyl[order(phyl$family, decreasing=TRUE),]
    par(mar=c(5, 8, 2, 2))
###Here we start the plotting section:
    plot(phyl$r, We plot the phylun ratios
         xlim=c(0,1.1*max(phyl$r)), ##" set the x limits which is the maximum ratio value
         ylim=c(0,nrow(phyl) +1 ),  ###set the y limit where it is the number of families +1
         xaxs="i",yaxs="i", ###These two values will removre the default offset of the axes, making them look better 
         type="n", ###Ask not to plot anything yet
         axes=F, ##Ask not to add any axes yet
         ylab="", xlab="Family ratio", 
         main=title ###title being the one that will be add in the command later
        )
####Set the axes values     
    axis(1, seq(0,100, 5),
         cex.axis=0.8, las=1)
####for the y axes, I am asking it to add a value in every second row so that it does not have to repeat the family name every second time    
    axis(2, seq(1, nrow(phyl),2), labels=unique(phyl$family), las=2, cex.axis=1)

###To make the barplots more clear, I first created background interchanging bars (zebra bars), by drawing a grey bar every second to fourth row  
    for(i in seq(2,nrow(phyl),4)){
        rect(0,
             i-2,
             max(phyl$r),
             i, col="lightgray", border="lightgray"
             )
    }
###Then I added the actual rectangles    
    for(j in 1:(nrow(phyl)){
        rect(0,
        (j-0.8),###this is to create some distance between each bar, the sum should be =1
        phyl$r[j],
        (j-0.2),
        col=phyl$colors[j],
        border="black"
        )
    }
    
    
###if you want to add a line for the base of the barplots    
    rect(0,0.2,0,nrow(phyl)-0.2, border="black")
    box()

###Next, creating the legend, it will be plotted based on some quite arbitrary modular coordinates that have been working quite fine so far for the majority of the plots 
        
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

###Here are a two examples on how to call this plot function
sapply("Annelida", function(x) plot.phylum(x, "number", paste("Annelida (Alignments vs Specimen Number)")))
sapply("Mollusca", function(x) plot.phylum(x, "biomass", paste("Mollusca (Alignments vs Specimen Biomass)")))
