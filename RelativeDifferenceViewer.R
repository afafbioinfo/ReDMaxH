#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (!require("optparse")) install.packages("optparse")

#library("optparse")

option_list = list(
  make_option(c("-R", "--RD1"), type="character", default=NULL, 
              help="Raw data Experiment 1", metavar="Path To file"),
  make_option(c("-S", "--RD2"), type="character", default=NULL, 
              help="Raw data Experiment 2", metavar="Path To file")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

#df = read.table(opt$file, header=TRUE)
#num_vars = which(sapply(df, class)=="numeric")
#df_out = df[ ,num_vars]
#write.table(df_out, file=opt$out, row.names=FALSE)

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("cowplot")) install.packages("cowplot")
if (!require("egg")) install.packages("egg")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("dplyr")) install.packages("dplyr")

RawC1<- read.csv(file =opt$RD1)
RawC2 <- read.csv(file =opt$RD2)

#MutationrateC1 <- read.csv(file =opt$MD1)
#MutationrateC2 <- read.csv(file = opt$MD2)

RawC1[RawC1== 0] <- NaN
RawC1[RawC1==Inf] <-NaN
RawC1 <- na.omit( RawC1 )

RawC2[RawC2== 0] <- NaN
RawC2[RawC2==Inf] <-NaN
RawC2 <- na.omit( RawC2 )

#######################Set the range ####################################""
#start<-20
#end <-227
#shift=14
#shift=37
#RawC1trimmed<- RawC1[ RawC1$i>start & RawC1$j>start &  RawC1$i<end & RawC1$j<end   ,]#&  RawC1$Rltdiff >2
#RawC1trimmed$i<- RawC1trimmed$i-start
#RawC1trimmed$j<- RawC1trimmed$j-start
#RawC1trimmed[RawC1trimmed== 0] <- NaN
#RawC1trimmed[RawC1trimmed==Inf] <-NaN
#RawC1trimmed <- na.omit( RawC1trimmed )

#RawC2trimmed<- RawC2[ RawC2$i>start & RawC2$j>start &  RawC2$i<end & RawC2$j<end ,]
#RawC2trimmed$i<- RawC2trimmed$i-start
#RawC2trimmed$j<- RawC2trimmed$j-start
#RawC2trimmed[RawC2trimmed== 0] <- NaN
#RawC2trimmed[RawC2trimmed==Inf] <-NaN
#RawC2trimmed <- na.omit( RawC2trimmed )

######################"""" Plot heatmap with muttaion rates:

flip<-RawC2 %>% 
  rename(
    i = j,
    j= i
  )
flip <-flip[c(2,1,3,4,5)]
Combine<-rbind(flip, RawC1)

###########palettes
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
sc <- scale_colour_gradientn(name ="R",colours = myPalette(100), limits=c(-1, 1))
sc1<- scale_colour_gradientn(name ="Raw count",colours = myPalette(100), limits=c(0, 10000))

# plot the heatmap
gg_hm = ggplot(data = Combine, aes(x=i, y=j)) + 
  geom_point(aes(colour = MM),size =1.5) +
  sc1+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  coord_equal()+ 
  theme(legend.position = "right")

gg_hm <-gg_hm+ggtitle("Raw mutation count. Tri-Sup: Experiment1, Tri-Inf: Experiment2 ")
#gg_hm
 ggsave("RedMaxPlots/RawMutationCount.png")
gg_hm1 = ggplot(data = Combine, aes(x=i, y=j)) + 
  geom_point(aes(colour = Rltdiff),size =1.5) +
  sc+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  coord_equal()+ 
  theme(legend.position = "bottom")

gg_hm1 <-gg_hm1+ggtitle(" R.diff. Tri-Sup: Experiment1, Tri-Inf: Experiment2 ")
#gg_hm1
ggsave("RedMaxPlots/RelativeDifference.png")

