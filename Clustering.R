
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (!require("optparse")) install.packages("optparse")


option_list = list(
  make_option(c("-R", "--R"), type="character", default=NULL, 
              help="Ranked helices", metavar="Path To file"),
  make_option(c("-C", "--C"), type="character", default=NULL, 
              help="Contributing pairs", metavar="Path To file"),
  make_option(c("-P", "--P"), type="character", default=NULL, 
              help="Maximal helices", metavar="Path To file"),
   make_option(c("-N", "--N"), type="integer", default=NULL, 
              help="Number of MH compounds to plot in arc diagram") ,
  make_option(c("-F", "--F"), type="character", default=NULL, 
        help="Fasta File Name")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("tidygraph")) install.packages("tidygraph")
if (!require("ggraph")) install.packages("ggraph")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("cowplot")) install.packages("cowplot")
if (!require("egg")) install.packages("egg")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("dplyr")) install.packages("dplyr")
if (!require("magrittr")) install.packages("magrittr")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("purrr")) install.packages("purrr")
if (!require("ggpubr")) install.packages("ggpubr")

contributingpairs<- read.csv(file = opt$C)
Rankedhelices<- read.csv(file = opt$R)
FASTA<-opt$F
target<-as.character(Rankedhelices$Helix)

M<-mean(Rankedhelices$averageRLtdiff)
sigma<-sd(Rankedhelices$averageRLtdiff)

Rankedhelices <-Rankedhelices  %>%
  mutate(Statt= case_when(Rankedhelices$averageRLtdiff >=M+2*sigma  ~ '>= M+2sigma',
                          Rankedhelices$averageRLtdiff >=M+1*sigma & Rankedhelices$averageRLtdiff <M+2*sigma ~ '>=M+sigma & <M+2sigma',
                          Rankedhelices$averageRLtdiff >=M & Rankedhelices$averageRLtdiff <M+1*sigma ~ '>=M & <M+sigma',
                          Rankedhelices$averageRLtdiff >=M-sigma & Rankedhelices$averageRLtdiff <M~ ' >=M-sigma& <M',
                          Rankedhelices$averageRLtdiff <M-sigma ~ ' <M-sigma& <M' ))

########################################"clusteing based on the mean
set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(Rankedhelices$averageRLtdiff, k, nstart = 10 )$tot.withinss
}
kmeans(Rankedhelices$averageRLtdiff, 1, nstart = 10 )
# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)
N<-0.1*max(wss_values)
thresh<-1+sum(  wss_values >=N )
thresh
L<-plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares (Mean)")

k <-kmeans(Rankedhelices$averageRLtdiff, centers=thresh) #Create thresh clusters, Remove columns 1 and 2
#length(k$cluster)
#length(target)

Numberofhelicesinoptimalcluster<-table( k$cluster)[names(table (k$cluster)) == k$cluster[1]]
k$clusters
Numberofhelicesinoptimalcluster
k3 <-kmeans(Rankedhelices$averageRLtdiff, centers=3) #Create3 clusters
k3$cluster
Numberofhelicesinoptimalcluster<-table( k3$cluster)[names(table (k3$cluster)) == k3$cluster[1]]
Numberofhelicesinoptimalcluster
Newdata<-data.frame(target,Rankedhelices$averageRLtdiff,k$cluster,Rankedhelices$Statt)
write.csv(Newdata,'RedMaxHoutput/Kmeans_elbow_mean.csv')
Newdata3<-data.frame(target,k3$cluster,Rankedhelices$Statt)

table(Newdata$k.cluster)[names(table(Newdata$k.cluster))==Newdata$k.cluster[1]]
table(Newdata3$k3.cluster)[names(table(Newdata3$k3.cluster))==Newdata3$k3.cluster[1]]
contributingpairsX<-merge(contributingpairs, Newdata, by.x="Helix", by.y="target")
contributingpairs3<-merge(contributingpairs, Newdata3, by.x="Helix", by.y="target")
contributingpairsX$Helix <- factor(contributingpairsX$Helix, levels=unique(target))
contributingpairs3$Helix <- factor(contributingpairs3$Helix, levels=unique(target))


p1<-ggplot(contributingpairsX,aes(x=1,y=Rforpairs,color=as.factor(k.cluster)))+#as.factor(Helix))) +
  ggtitle("Kmeans on Mean, nbr clusters (elbow)")+
  geom_violin(colour = "grey50")+
  geom_boxplot(width=0.1,colour = "grey50")+
  stat_summary(fun.y=mean, geom="point", size=2, color="black")+
  geom_jitter( size = 1.5)+ylim(-1,3)+
  
  facet_wrap(~as.factor(contributingpairsX$Helix),nrow = 1)

p1<-p1+ theme(legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x=element_blank())

p2<-ggplot(contributingpairs3,aes(x=1,y=Rforpairs,color=as.factor(k3.cluster)))+#as.factor(Helix))) +
  ggtitle("Kmeans on Mean, nbr clusters=3")+
  geom_violin(colour = "grey50")+
  geom_boxplot(width=0.1,colour = "grey50")+
  stat_summary(fun.y=mean, geom="point", size=2, color="black")+
  geom_jitter( size = 1.5)+ylim(-1,3)+
  
  facet_wrap(~as.factor(contributingpairs3$Helix),nrow = 1)

p2<-p2+ theme(legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x=element_blank())

g<-ggplot(contributingpairsX,aes(x=1,y=Rforpairs,color=as.factor(Rankedhelices.Statt)))+#as.factor(Helix))) +
  ggtitle(">=Mean+2sigma,..")+
  geom_violin(colour = "grey50")+
  geom_boxplot(width=0.1,colour = "grey50")+
  stat_summary(fun.y=mean, geom="point", size=2, color="black")+
  geom_jitter( size = 1.5)+ylim(-1,3)+
  
  facet_wrap(~as.factor(contributingpairsX$Helix),nrow = 1)

g<-g+ theme(legend.title= element_blank(),legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x=element_blank())
library(ggpubr)

multi.page <- ggarrange(p1, p2, g,
                        nrow =9, ncol = 1)
ggexport(multi.page, filename =paste( "RedMaxPlots/",FASTA,"_Clusters_Helices.pdf"), width =30.0, height = 15)



ggexport(p1, filename = paste("RedMaxPlots/",FASTA,"_Elbow_Clusters_Helices.pdf"), width =30.0, height = 15)



##############" TODO  bring back the mypalette to the first Spectral values
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
sc <- scale_colour_gradientn(name ="Average R",colours = myPalette(100), limits=c(-1, 1))


p1<-ggplot(contributingpairsX,aes(x=1,y=Rforpairs,colour=Rankedhelices.averageRLtdiff))+#as.factor(Helix))) +
  sc+
  geom_violin(colour = "grey50")+
  geom_boxplot(width=0.1,colour = "grey50")+
  stat_summary(fun.y=mean, geom="point", size=2, color="black")+
  geom_jitter( size = 3)+ylim(-1,3)+
  facet_wrap(~as.factor(contributingpairsX$Helix),nrow = 1)

p1<-p1+ theme(legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x=element_blank(),axis.text=element_text(size=25,face="bold"),legend.title = element_text( size=30, face="bold"),legend.text = element_text( size=25, 
                                     ),legend.key.size = unit(2, 'cm'))


ggexport(p1, filename = paste("RedMaxPlots/",FASTA,"_Color_by_average_Helices.pdf"), width =30.0, height = 15)


myPalette <- colorRampPalette(brewer.pal(11, "RdBu"))
sc <- scale_colour_gradientn(name ="Average R",colours = myPalette(100), limits=c(-1, 1))


p1<-ggplot(contributingpairsX,aes(x=1,y=Rforpairs,colour=Rankedhelices.averageRLtdiff))+#as.factor(Helix))) +
  sc+
  geom_violin(colour = "grey50")+
  geom_boxplot(width=0.1,colour = "grey50")+
  stat_summary(fun.y=mean, geom="point", size=2, color="black")+
  geom_jitter( size = 3)+ylim(-1,3)+
  facet_wrap(~as.factor(contributingpairsX$Helix),nrow = 1)

p1<-p1+ theme(legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x=element_blank(),axis.text=element_text(size=25,face="bold"),legend.title = element_text( size=30, face="bold"),legend.text = element_text( size=25, 
                                     ),legend.key.size = unit(2, 'cm'))


ggexport(p1, filename = paste("RedMaxPlots/",FASTA,"_Color_by_average_Helices_RdBupalette.pdf"), width =30.0, height = 15)
####################Plotting arc diagrams
BPs<- read.csv(file = opt$P)

################## Get all contributing BPs for a given helix
library(tidyverse)
df_total = data.frame()
for (row in 1:nrow(BPs)) {
  Helix <- BPs[row, "Helix"]
  opening  <- BPs[row, "opening"]
closing <- BPs[row, "closing"]
L<- BPs[row, "length"]-1
df<-bind_cols(Helix,seq(opening, opening+L, by=1 ),seq(closing, closing-L, by=-1 ))
df_total <- rbind(df_total,df)
}
colnames(df_total)=c("Helix", "i","j")
############################ split helix compounds
df_rankedHelices = data.frame()
for (row in 1:nrow(Rankedhelices)) {
  Helix <- Rankedhelices[row, "Helix"]
 Rlt<-Rankedhelices[row, "averageRLtdiff"]
  df <-bind_cols(str_split(Helix,"-"), Rlt,row)
  df_rankedHelices<- rbind( df_rankedHelices,df)
}
colnames( df_rankedHelices)=c("Helix", "Rd","order")
#head( df_rankedHelices)

############################## Combine  df_total with  df_rankedHelices based on values in ranked helices
WeightedBasepairs<-merge(df_rankedHelices, df_total, by.x="Helix", by.y="Helix")
#Numberofcompoundsforarcplot=12
#head(WeightedBasepairs[WeightedBasepairs$order<=Numberofcompoundsforarcplot,])

arc_theme <-  theme(axis.line=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                    panel.grid.major.x = element_line( size=.1, color="grey" ),
                     plot.background=element_blank(),
                     plot.title = element_text(family="Helvetica Neue Light", size=20, face="plain"),
                     plot.subtitle =  element_text(family="Helvetica Neue Light", size=12.9, face="plain"),
                     legend.text =  element_text(family="Helvetica Neue Light", face="plain"),
                     legend.title = element_text(family="Helvetica Neue Light", face="plain"),
                     axis.text = element_text(family="Helvetica Neue Light", face="plain", size=6),
                     legend.key = element_blank())

Numberofcompoundsforarcplot=opt$N

SelectedHelices<-df_rankedHelices[df_rankedHelices$order<=Numberofcompoundsforarcplot,]


LabelsHelices<-merge(SelectedHelices, BPs, by.x="Helix", by.y="Helix")

myPalette <- colorRampPalette(brewer.pal(11, "RdBu"))

Parcplots<-ggraph(WeightedBasepairs[WeightedBasepairs$order<=Numberofcompoundsforarcplot & ( WeightedBasepairs$Rd >=0.1 |WeightedBasepairs$Rd <=-0.1) ,], layout = 'linear') + 
  geom_edge_arc(aes( x=WeightedBasepairs[WeightedBasepairs$order<=Numberofcompoundsforarcplot & ( WeightedBasepairs$Rd >=0.1 |WeightedBasepairs$Rd <=-0.1),]$i, xend=WeightedBasepairs[WeightedBasepairs$order<=Numberofcompoundsforarcplot & ( WeightedBasepairs$Rd >=0.1|WeightedBasepairs$Rd <=-0.1),]$j,colour=WeightedBasepairs[WeightedBasepairs$order<=Numberofcompoundsforarcplot & ( WeightedBasepairs$Rd >=0.1 |WeightedBasepairs$Rd <=-0.1),]$Rd), 
                fold = TRUE)+#, strength = 1)+#abs(WeightedBasepairs[WeightedBasepairs$order<=Numberofcompoundsforarcplot & ( WeightedBasepairs$Rd >=0.1 |WeightedBasepairs$Rd <=-0.1),]$Rd))+
 # geom_text(data=WeightedBasepairs[WeightedBasepairs$order<=Numberofcompoundsforarcplot,],
 #   label=WeightedBasepairs[WeightedBasepairs$order<=Numberofcompoundsforarcplot,]$Helix, 
  #  x = WeightedBasepairs[WeightedBasepairs$order<=Numberofcompoundsforarcplot,]$i, y = 20, 
 #   check_overlap = T
 # )+
 # scale_x_continuous(expand = c(0, 0),limits = c(1,120))  +
  scale_edge_colour_gradientn(name ="Rd average",colours = myPalette(100), limits=c(-1, 1))+
  theme_minimal()+
  scale_x_continuous(breaks = seq(0,max(WeightedBasepairs[WeightedBasepairs$order<=Numberofcompoundsforarcplot &( WeightedBasepairs$Rd >=0.1 |WeightedBasepairs$Rd <=-0.1),]$j)+10, by = 10))+
  arc_theme
 


if(Numberofcompoundsforarcplot<10){
 Parcplots<-  Parcplots+
        geom_label(data=LabelsHelices,
              label=LabelsHelices$Helix, 
             x = LabelsHelices$opening+LabelsHelices$length/2, y = 1.5*log2(LabelsHelices$length), 
         hjust = rep(-1:1, length.out=length(LabelsHelices$opening)),
            vjust = rep(-0.6:0.6, length.out=length(LabelsHelices$opening)),
            #angle = c( 45, -45),
            #check_overlap = T,
              label.padding = unit(0.55, "lines"), # Rectangle size around label
             label.size = 0.35,
  #            #color = "black",
  #           size=4,
            fill="grey",alpha=0.1
  )
  }
  
ggexport(Parcplots, filename = paste("RedMaxPlots/",FASTA,"_Arc_diagram_Color_by_average_Helices.eps"), width =30, height = 5)


