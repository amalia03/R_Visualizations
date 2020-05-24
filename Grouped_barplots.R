###This Rscript is a specialised script that will plot the families per given phylum along, coloredby type of dataset ( whether it has been contaminatied by a sourc or not)
###The input file needs to have the mentioned elements: 
###sk.ph -> dataset
###"phylum"= Phylum, "n"= number of phyla per family (found at a previous analysis), "family"= taxonomic family
###"acc"= accession id (signifying the different genes)

###Before the plot, I first start with picking colors and creating some useful functions beforehand.
###Since in our case there are four categories of data ( lumpfish contaminated, multiple phyla contaminated, both or neither), I needed to pick four colors so I poicked four I liked manually 
br.col <- c("goldenrod2","cyan4", "tomato2", "dodgerblue")

###Here a function will return the color from the array that fits with the numbered category of the variable n. N in this case has to have a numeric label for each group (in this case -1 for both failing", 0 for lumpfish, 1 for the sequences that pass the filters and >1 for multiple phyla). Obviously those categories can be changed depending on the situation
color.picker.raw <- function(x){                                                                                     
    if(x == -1 ){return(br.col[1])}                                                                            
            else if(x == 0 ){return(br.col[2])}
            else if(x == 1 ){return(br.col[3])}
            else {return(br.col[4])}                                                                                               
}
###Use sapply to create a new dataset variable that includes that has the color for each category
sk.ph$cols <- sapply(sk.ph$n, color.picker.raw)                                                                   

###Create a label for each category by using cut.
sk.ph$category <- cut(sk.ph$n, c(-Inf, -1, 0, 1,Inf),
                      c("Failed mp+lf", "Failed lf", "Passed mp+lf","Failed mp"))  
    
###This is a useful function that will be used for the xaxis ticks so that there will always be a reasonable and legible amount of ticks, aka when values are 1-10 , the increments increase by 1, if above 50 they increase by 10 and anything in between , increase by 5.
    
ticks <- function(x){
    if(x <= 10)
    {return(1)}
    else if(x>50)
    {return(10)}
    else{return(5)}
}
    
###Here we start with the plot function
phylo.bar <- function(p, leg.h=6, leg=TRUE){
    ##where p stands for phylum, leg.h stands for legend height and leg stands for a legend being used
###Here we pick the phylum of interest. I would usually use the command below but I think grep would pick any string that contains the elements in general, or it didnt work otherwise.
    # phyl.all<- sk.bl.2[sk.bl.2$phylum == p,]
    phyl.all<- sk.ph[grep(p, sk.ph$phylum),]
    
    ##First see how many species are there per family by aggregated the acc column per family(how many alignments we have per family)
    fam.ph<- aggregate(phyl.all$acc, by=list(phyl.all$family), length)
    colnames(fam.ph)<- c("family", "Freq")
    ###We merge the intial dataset so that we now have the Frequency column 
    phyl.all<- merge(fam.ph, phyl.all, by=c("family"))
    ###We only need one value per family and since all categories should be the same as per family, we only need a single representative per family 
    phyl<- phyl.all[!duplicated(phyl.all$family),]
    ###Then I order per frequency so that when they are plotted, they are sorted by number of alignments
    phyl <- phyl[order(phyl$Freq, decreasing=FALSE),]
    
    ###Set the plotting stage
    plot(phyl$Freq, xlim=c(0,eval(max(phyl$Freq)+(0.1*max(phyl$Freq)))),
         ylim=c(0,nrow(phyl)+ 1) , ##Axes limits are the max number of frequencis for x-axis and number of family groups (or selected phylum rows) for xaxis
         type="n", ##to not plot anything 
         axes=F, ###to leave the axes for the next command
         ylab="", xlab="Number of alignments", 
         main=unique(phyl$p)) ###Use the phylum column to name the plot 
    ###for the xaxis I use the sequence function to say that I want to count from 0 to the maximum frequency of any family + 0.1 of that number and label on what the tick function will return based on the magnitude of the frequency value
    axis(1, seq(0,round(eval((max(phyl$Freq)+(0.1*max(phyl$Freq))))),
                by=eval(ticks(max(phyl$Freq)))),
         cex.axis=0.8, las=1)
    ###yaxis is more self explained, counts from 1 to the number of families, labeled by family 
    axis(2, 1:nrow(phyl), labels=phyl$family, las=2, cex.axis=1)
        ##So the square would be as long as all the families, as wide as the maximum score of alignments

    ###Here I process the groups so that they each become an individual rectangle that stacks as a whole bar per family
    ###Starting with the largest category ( so that the order of categories ar always consistent) and then the frequency:
  ###Loop starts by taking the each family by entry number and making a category rectangle
    for(j in 1:(nrow(phyl[order(phyl$n, phyl$Freq, decreasing=FALSE),]))){
        fam<- phyl[j,"family"]
        sk.hom <- phyl.all[phyl.all$family== fam,]
        sk.hom$n[sk.hom$n >=2] <- 2 ##This could be done before but regardless, for the multiple phyla category, the variable is changed from everything >1 to 2
        fam.freq=as.data.frame(table(sk.hom$n)) ##make a frequency table for each category
        colnames(fam.freq) <- c("id","freq")
        phyl[phyl$family==fam,"Freq"] ###
        ##You dont need a percent dumbass..
        ###Now in reference to each category (id), we create a dataframe with where each category has a designated color
        colindex <- data.frame("colors"=c("goldenrod2","cyan4","tomato2", "dodgerblue"),"id"=c(-1,0,1,2))

    ###Merge the new colindex dataframe to the frequency table for each category that we made before. 
        fam.freq<- merge(colindex, fam.freq, by=c("id"))
    ###Sometimes the colors do not work because they have been taken in as levels, so in case that has happened, we change the \m into characters
        fam.freq$colors<- as.character(fam.freq$colors)
   ###Barplotting time. We start with the first rectangle that is separate because the scripting specifics are different
   ###So what it does is it makes the first category value and uses the frequency for the variable length of the rectangle
        rect(0,j-0.4,(fam.freq$freq[1]),j+0.4, col=fam.freq$colors[1], border=fam.freq$colors[1])
   ###The next part can be a loop, since we are using vectos based on the previous and the current category. ie: in the command above we created rectangles with lenghth being frequency of category 1. Here the same goes but the frequency is category 2 - category 1 , next one is category 3 - category 2 ect. 
            for(i in 1:(length(fam.freq$freq)-1)){ ##-1 is used to make the vector 1 less since we taken out the first value in the command before
            rect(sum(fam.freq$freq[1:i]), ##the offset should start at the sum of all previous categories
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
