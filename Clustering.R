
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (!require("optparse")) install.packages("optparse")


option_list = list(
  make_option(c("-R", "--R"), type="character", default=NULL, 
              help="Ranked helices", metavar="Path To file"),
  make_option(c("-C", "--C"), type="character", default=NULL, 
              help="Contributing pairs", metavar="Path To file")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

if (!require("ggplot2")) install.packages("ggplot2")
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
ggexport(multi.page, filename = "RedMaxPlots/Clusters_Helices.pdf", width =30.0, height = 15)



ggexport(p1, filename = "RedMaxPlots/Elbow_Clusters_Helices.pdf", width =30.0, height = 15)




myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
sc <- scale_colour_gradientn(name ="Average R",colours = myPalette(100), limits=c(-1, 1))
p1<-ggplot(contributingpairsX,aes(x=1,y=Rforpairs,colour=Rankedhelices.averageRLtdiff))+#as.factor(Helix))) +
  sc+
  ggtitle("Kmeans on Mean, nbr clusters (elbow)")+
  geom_violin(colour = "grey50")+
  geom_boxplot(width=0.1,colour = "grey50")+
  stat_summary(fun.y=mean, geom="point", size=2, color="black")+
  geom_jitter( size = 3)+ylim(-1,3)+
  
  facet_wrap(~as.factor(contributingpairsX$Helix),nrow = 1)

p1<-p1+ theme(legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x=element_blank())

ggexport(p1, filename = "RedMaxPlots/Color_by_average_Helices.pdf", width =30.0, height = 15)

