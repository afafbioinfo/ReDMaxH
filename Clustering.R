RawC1<- read.csv(file ='/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication/5SFREQUENCY_COUNTex_Rawdata.csv')
RawC2 <- read.csv(file ='/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication/5SincelFREQUENCY_COUNTex_Rawdata.csv')

helixC1<- read.csv(file ='/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication/5SFREQUENCY_COUNTex_Norm_helices.csv')
helixC2<- read.csv(file ='/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication/5SincelFREQUENCY_COUNTex_Norm_helices.csv')

MutationrateC1 <- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication/5SFREQUENCY_COUNTex_MutationIndex.csv')
MutationrateC2 <- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication/5SincelFREQUENCY_COUNTex_MutationIndex.csv')

####################################"""By helix"
head(RawC1)
head(helixC1)
RawC1<-merge(RawC1,helixC1, by=c("i","j","Rltdiff")) 
RawC2<-merge(RawC2,helixC2, by=c("i","j","Rltdiff")) 

###########################################################
RawC1[RawC1$Rltdiff>=3, "Rltdiff"] <- 3
RawC1[RawC1==Inf] <-NaN
#RawC1 <- na.omit( RawC1 )

RawC2[RawC2$Rltdiff>=3, "Rltdiff"] <- 3
RawC2[RawC2==Inf] <-NaN
RawC2 <- na.omit( RawC2 )

#start<-16
##RawC1$i<- RawC1$i-start
#RawC1$j<- RawC1$j-start

RawC2$i<- RawC2$i-start
RawC2$j<- RawC2$j-start

MutationrateC2$position<-MutationrateC2$position-start
MutationrateC1$position<-MutationrateC1$position-start



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
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)   
flip<-RawC2 %>% 
  rename(
    i = j,
    j= i
  )
flip <-flip[c(2,1,3,4,5,6,7,8,9)]
Combine<-rbind(flip, RawC1)

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

###########palettes
library(ggplot2)
library(RColorBrewer)
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
sc <- scale_colour_gradientn(name ="Rltdiff",colours = myPalette(100), limits=c(-1, 1))

#sc<-  scale_colour_gradient2(
#  low = "deeppink3",
#  mid = "lightcyan",
 # high = "blue",
#  midpoint = 0.1,
#  space = "Lab",
#  na.value = "grey50",
#  guide = "edge_colourbar"
#)
sc1<- scale_colour_gradientn(name ="Raw count",colours = myPalette(100), limits=c(0, 10000), labels =fancy_scientific)
library(ggplot2)
# plot the heatmap
setEPS()
postscript("MMcount.eps")

gg_hm = ggplot(data = RawC1, aes(x=i, y=j)) + 
  geom_point(aes(colour = MM),size =1.5) +
  geom_rect(aes(xmin=1, xmax=15, ymin=1,ymax=120), fill="grey", alpha=0.1)+
  geom_rect( aes(ymin=102, ymax=120, xmin=1,xmax=120), fill="grey", alpha=0.1)+
  sc1+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  coord_equal()+  theme(legend.position = "left")+
  coord_cartesian(xlim=c(1,120), ylim = c(1,120))+  scale_y_continuous(breaks = c(1,20,40,60,80,100,120))+  scale_x_continuous(breaks = c(1,20,40,60,80,100,120))+
  theme(legend.position = "right")
# scale_x_continuous(limits = c(1, 115)) + 

#+ theme(legend.position = "none")+ 
# scale_z_continuous(labels=fancy_scientific) 
#theme(legend.position = "right")

#gg_hm +ggtitle("Raw count  ")

gg_hm
dev.off()
library(cowplot)
library(grid)
library(gridExtra)
setEPS()
postscript("Rcountlegend.eps")
legend <- cowplot::get_legend(gg_hm1)
grid.newpage()
grid.draw(legend)
dev.off()

flip<-RawC1 %>% 
  rename(
    i = j,
    j= i
  )
flip <-flip[c(2,1,3,4,5,6,7,8,9)]
Combine<-rbind(flip, RawC2)
setEPS()
postscript("Rtriagular.eps")
gg_hm1 = ggplot(data = RawC1, aes(x=i, y=j)) + 
  geom_point(aes(colour = Rltdiff),size =1.5) +
#  geom_rect(aes(ymin=1, ymax=15, xmin=1,xmax=120), fill="grey", alpha=0.1)+
#geom_rect( aes(xmin=102, xmax=120, ymin=1,ymax=120), fill="grey", alpha=0.1)+
  sc+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  coord_equal()+  theme(legend.position = "right")+
  coord_cartesian(xlim=c(1,120), ylim = c(1,120))+  scale_y_continuous(breaks = c(1,20,40,60,80,100,120))+  scale_x_continuous(breaks = c(1,20,40,60,80,100,120))+
  theme(legend.position = "none")
gg_hm1
dev.off()



library(patchwork)
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pembedded<-blankPlot + inset_element(gg_hm, left=-0.08, bottom = 0.05, right = 1, top = 1) + inset_element(gg_hm1, left=-0.08, bottom = 0.05, right = 1, top = 1)


mycols <- c( "#52854C","#C3D7A4", "#FFDB6D", "#C4961A",
             "#D16103","red" )
count.data <- data.frame(
  class = c("4", "5", "6", "7", "8",  "10"),
  n = c(23,6,1,0,2,1),prop=  c(23/33,6/33,1/33,0/33,2/33,1/33))
#after merging
count.data <- data.frame(
  class = c("4", "5", "6", "7", "8",  "10"),
  n = c(23-7,6+1,1,0,2+1,1),prop=  c(16/28,7/28,1/28,0/28,3/28,1/28))
setEPS()
postscript("28helicescircle.eps")
ggplot(count.data, aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  # geom_text(aes( label=n), color = "black")+ #label = format(round(prop*100, 0), nsmall =0))"%"
  scale_fill_manual(values = mycols) +
  theme_void()+
  xlim(0.5, 2.5)
dev.off()
Mutations5Scellfree<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/5SFREQUENCY_COUNTex_Norm_helices.csv')

v<-S5helices[S5helices$length>4,]$Helix
#v=[C5  C6  C7  C10 C14 C19 C22 C26 C27 C31]
v<-S5helices[S5helices$length==4,]$Helix
library(dplyr)
Mutations5Scellfree <-Mutations5Scellfree[Mutations5Scellfree$Helix=="C26" |Mutations5Scellfree$Helix=="C63",]
Mutations5Scellfree <-Mutations5Scellfree[Mutations5Scellfree$Helix=="C64" |Mutations5Scellfree$Helix=="C32",]#|Mutations5Scellfree$Helix=="C40"|Mutations5Scellfree$Helix=="C43",]



#Mutations5Scellfree <-Mutations5Scellfree[Mutations5Scellfree$Helix=="C5" | Mutations5Scellfree$Helix=="C14"| Mutations5Scellfree$Helix=="C19"|  Mutations5Scellfree$Helix=="C31",]
#Mutations5Scellfree <-Mutations5Scellfree[Mutations5Scellfree$Helix!="C5" & Mutations5Scellfree$Helix!="C6"& Mutations5Scellfree$Helix!="C7"& Mutations5Scellfree$Helix!="C10"& Mutations5Scellfree$Helix!="C14"& Mutations5Scellfree$Helix!="C19"& Mutations5Scellfree$Helix!="C22"& Mutations5Scellfree$Helix!="C26"& Mutations5Scellfree$Helix!="C27"& Mutations5Scellfree$Helix!="C31",]


#write.csv(Mutations5Scellfree ,'/home/user/Documents/Alain/longhditributions.csv')
dataB <- Mutations5Scellfree[, c("i", "j", "Rltdiff","Helix")]
dataB[dataB$Rltdiff>=3, "Rltdiff"] <- 3
setEPS()

postscript("/home/user/Documents/Alain/figure2/64-32Rltddiffallhelices.eps")
library(ggraph)#
ggraph(dataB, layout = 'linear') + 
  geom_edge_arc(aes( x=dataB$i, xend=dataB$j,colour=Rltdiff), alpha=2, 
                fold = TRUE)+
  scale_x_continuous(expand = c(0, 0),limits = c(1,120))  +
  scale_edge_colour_gradientn(name ="Rltdiff",colours = myPalette(100), limits=c(-1, 1))+
 theme_minimal()+
  coord_cartesian(xlim=c(1,120))+   scale_x_continuous(breaks = c(1,20,40,60,80,100,120))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank())+
  theme(legend.position = "none")
#p+scale_x_reverse()
dev.off()
#################demi-plan
Mutations5Scellfree<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/5SFREQUENCY_COUNTex_Norm_helices.csv')

Mutations5Scellfree <-Mutations5Scellfree[Mutations5Scellfree$Helix=="C26"|Mutations5Scellfree$Helix=="C63"|Mutations5Scellfree$Helix=="C32"|Mutations5Scellfree$Helix=="C43"|Mutations5Scellfree$Helix=="C64" |Mutations5Scellfree$Helix=="C40",]
#Mutations5Scellfree <-Mutations5Scellfree[Mutations5Scellfree$Helix=="C5" | Mutations5Scellfree$Helix=="C14"| Mutations5Scellfree$Helix=="C19"|  Mutations5Scellfree$Helix=="C31",]
#Mutations5Scellfree <-Mutations5Scellfree[Mutations5Scellfree$Helix!="C5" & Mutations5Scellfree$Helix!="C6"& Mutations5Scellfree$Helix!="C7"& Mutations5Scellfree$Helix!="C10"& Mutations5Scellfree$Helix!="C14"& Mutations5Scellfree$Helix!="C19"& Mutations5Scellfree$Helix!="C22"& Mutations5Scellfree$Helix!="C26"& Mutations5Scellfree$Helix!="C27"& Mutations5Scellfree$Helix!="C31",]


#write.csv(Mutations5Scellfree ,'/home/user/Documents/Alain/longhditributions.csv')
dataB <- Mutations5Scellfree[, c("i", "j", "Rltdiff","Helix")]
dataB[dataB$Rltdiff>=3, "Rltdiff"] <- 3

library(ggraph)#
ggraph(dataB, layout = 'linear') + 
  geom_edge_arc(aes( x=dataB$i, xend=dataB$j,colour=Rltdiff), alpha=2, 
                fold = TRUE)+
  scale_x_continuous(expand = c(0, 0),limits = c(1,120))  +
  scale_edge_colour_gradientn(name ="Rltdiff",colours = myPalette(100), limits=c(-1, 1))+
  theme_minimal()+
coord_cartesian(xlim=c(1,120))+   scale_x_continuous(breaks = c(1,20,40,60,80,100,120))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank())+
theme(legend.position = "none")

###########################""jitter plots:
dataB <- Mutations5Scellfree[, c("i", "j","ijdist", "Rltdiff","Helix")]
dataB[dataB$Rltdiff>=3, "Rltdiff"] <- 3
p<-ggplot(mapping = aes(x=1,y=C64C63 )) +
  geom_boxplot(colour = "grey50")+
  geom_jitter( size = 1.5)+ylim(-1,3)
p<-p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x=element_blank())

  # xlim(5, 20)+
  # xlim(45, 180)+
  # xlim(5,32)+#xlim(5,32)+
  # ylim(-1.1,3.1)+
  facet_wrap(~Helix,nrow = 1)+
  # scale_color_manual(values=c( "red", "blue","pink","brown2", "darkgreen", "black","darkmagenta","chartreuse","blue"))
  #scale_color_manual(values=c("deeppink", "brown2","darkolivegreen4","purple"))
  scale_color_manual(values=c("cyan","#44AA99","cyan", "red", "green","orange","deeppink","brown2","cyan", "deeppink","brown2","darkgreen","darkolivegreen4","purple", "black","darkmagenta","chartreuse"))

C40C43 <-c(0.06401216023179604, 0.007440845523367854, -0.05178085664068233, 0.05935159319326254, 0.023140238141638726, -0.014430455102838993, 0.03534819304346376, -0.014672127983128287, 0.02304801540009173, 0.11100046544483311, 0.06991620621239983, -0.05353748449274265, 0.06716649707417836, -0.03839824030477869, -0.10062888935953593, 0.011337157986587151, -0.03249578591955345, -0.20977143364594888, -0.050123729844174104, 0.09522630941981307, 0.008266945926825182, -0.19120357758206116, 0.1421355520163104, 0.16270995734268104, -0.377323282366379, 0.08811143076582877, 0.14566395658453607, -0.1886437211258854, 0.11052506117522992, -0.20861920036741882, -0.37882414417564686, 0.001945065297736542, 0.0033556623539415533, -0.052860934921876664)
C64C63 <-c(-0.030170544925281224, 0.1364793436139271, -0.11908020330567691, -0.09554253264348327, 1.9543964470455668, -0.3384169010289339, 0.13528515540443484, 0.24267602216943157, -0.40888495461873364, 0.32794746491007043, 0.46070993959162754, 0.46122148372184113, 0.35347563184161107, 0.5797537687605101, 0.04310584956247029, -0.09714449962058225, 1.214563077407716, 0.027585368568469423, -0.09822521632866611, 1.1400197110704475, -0.09645561204983498, 0.5643089185825565, 1.4049247108710567, 1.9549190170414812, 3, 1.7286128831738234, 0.17771979559182197, 2.0789399757673617, 0.27003110501282895, 0.848152908072343, 0.7121564665005605, 0.022548741844794034, 0.020607431804457538, -0.10982697727573933, -0.2545214734765022, 2.5262055216780435, -0.16144394185979025, -0.035127951577609805, -0.08153540917334691, 0.17930052019456091, -0.020139786343241352, -0.016178883127994582, 0.27746671235714315)
C26C32 <-c(0.21633112454421333, 0.14725541288824942, -0.05662870051020537, 0.08930778391384894, 0.1891184352508092, -0.09392389090262307, 0.2682154237358921, 0.47703352173652985, -0.07269669874501854, 3, 0.10524523440072397, 0.18547021111675208, -0.02353274175235899, 1.2599685311603652, 0.04822963904388631, 0.34161116823939386, 0.14930039899865466, 0.15528408022272402, 0.5844358880376115, 0.2519013361586909, -0.4457366662738452, -0.16420707173060464, 0.04328753302693621, -0.07943569676051317, 0.07931530315471377, 0.10697606556761517, -0.07503477900688285, -0.03821230086919272, 0.05993951962898181, 0.030108087675607513, -0.08933358964426091, 0.04687053997813832, -0.021537922200446306, 0.3444094146926248, 0.20588950109504484, -0.0168338897669552, -0.2927628483539176, -0.15311839940802946, -0.17038639795535793, 0.32083986850644486, 0.0019685661014240694, -0.1520051976520514, 0.07124085289303336, -0.23637312514936124, -0.009037284521788534, -0.2154987291107246, -0.014845810516453095, -0.34613879824914817, -0.34235056929670804, -0.1016187920241955, 0.03266124367003676, -0.0379667423886662, -0.23487117610591526)
C21C24 <-c(-0.12638101140121807, -1.0, -0.05204196560498434, -0.09626249703313872, 0.0709125714701659, -0.18940479993562628, -0.11042751694646463, 0.00726287837032177, -0.31656070224497107, -0.2307454815100855, -0.16016615758243866, -1.0, -1.0, -1.0, -0.5133488844147525, 0.12864787639752948, -0.05882939069427924, -0.05364385162791299, -1.0, -0.29968491160191, -0.26121019262912554, 0.5129667266156068, 0.10880818515197807, 0.4379832746273899, 0.4578336686592514, 0.6410152438624929, 0.04505745487156291, 0.009320840498865012, -0.020523154351014667, 0.15769522391129656, -0.06258071892051822, -0.32998771255367343)
C51C52<-c(-0.026800743748199102, 0.06143781619783842, 0.2638120455179165, -0.22637137386687595, -0.11251678451158974, 0.16285541165079287, 0.3448661028020803, 0.11032425579585989, -0.202603672236863, 0.20275884364997449, 0.021536080257503056, 0.09769910839919703, 0.04065443581198678, 0.016374516757634818, -0.04941672542098832, 0.004623650619335018, -0.03698279648337918, 0.04542971088281707, 0.029894380679541355, 0.006709929607947772, 0.1259699234656209, -0.13024251384970847, 0.011410810975428844, -0.1649990258899369, -0.05988959070358722, -0.0077795105884423486, 0.2589272069617304, 0.05920341373526867, 0.02012209142441832, 0.04875346310567188, -0.003247853249949209, -0.019200657794131867, -0.0036123823687246345, 0.0929339537345263)

R<-c(0.06401216023179604, 0.007440845523367854, -0.05178085664068233, 0.05935159319326254, 0.023140238141638726, -0.014430455102838993, 0.03534819304346376, -0.014672127983128287, 0.02304801540009173, 0.11100046544483311, 0.06991620621239983, -0.05353748449274265, 0.06716649707417836, -0.03839824030477869, -0.10062888935953593, 0.011337157986587151, -0.03249578591955345, -0.20977143364594888, -0.050123729844174104, 0.09522630941981307, 0.008266945926825182, -0.19120357758206116, 0.1421355520163104, 0.16270995734268104, -0.377323282366379, 0.08811143076582877, 0.14566395658453607, -0.1886437211258854, 0.11052506117522992, -0.20861920036741882, -0.37882414417564686, 0.001945065297736542, 0.0033556623539415533, -0.052860934921876664,-0.030170544925281224, 0.1364793436139271, -0.11908020330567691, -0.09554253264348327, 1.9543964470455668, -0.3384169010289339, 0.13528515540443484, 0.24267602216943157, -0.40888495461873364, 0.32794746491007043, 0.46070993959162754, 0.46122148372184113, 0.35347563184161107, 0.5797537687605101, 0.04310584956247029, -0.09714449962058225, 1.214563077407716, 0.027585368568469423, -0.09822521632866611, 1.1400197110704475, -0.09645561204983498, 0.5643089185825565, 1.4049247108710567, 1.9549190170414812, 3, 1.7286128831738234, 0.17771979559182197, 2.0789399757673617, 0.27003110501282895, 0.848152908072343, 0.7121564665005605, 0.022548741844794034, 0.020607431804457538, -0.10982697727573933, -0.2545214734765022, 2.5262055216780435, -0.16144394185979025, -0.035127951577609805, -0.08153540917334691, 0.17930052019456091, -0.020139786343241352, -0.016178883127994582, 0.27746671235714315,0.21633112454421333, 0.14725541288824942, -0.05662870051020537, 0.08930778391384894, 0.1891184352508092, -0.09392389090262307, 0.2682154237358921, 0.47703352173652985, -0.07269669874501854, 3, 0.10524523440072397, 0.18547021111675208, -0.02353274175235899, 1.2599685311603652, 0.04822963904388631, 0.34161116823939386, 0.14930039899865466, 0.15528408022272402, 0.5844358880376115, 0.2519013361586909, -0.4457366662738452, -0.16420707173060464, 0.04328753302693621, -0.07943569676051317, 0.07931530315471377, 0.10697606556761517, -0.07503477900688285, -0.03821230086919272, 0.05993951962898181, 0.030108087675607513, -0.08933358964426091, 0.04687053997813832, -0.021537922200446306, 0.3444094146926248, 0.20588950109504484, -0.0168338897669552, -0.2927628483539176, -0.15311839940802946, -0.17038639795535793, 0.32083986850644486, 0.0019685661014240694, -0.1520051976520514, 0.07124085289303336, -0.23637312514936124, -0.009037284521788534, -0.2154987291107246, -0.014845810516453095, -0.34613879824914817, -0.34235056929670804, -0.1016187920241955, 0.03266124367003676, -0.0379667423886662, -0.23487117610591526,-0.12638101140121807, -1.0, -0.05204196560498434, -0.09626249703313872, 0.0709125714701659, -0.18940479993562628, -0.11042751694646463, 0.00726287837032177, -0.31656070224497107, -0.2307454815100855, -0.16016615758243866, -1.0, -1.0, -1.0, -0.5133488844147525, 0.12864787639752948, -0.05882939069427924, -0.05364385162791299, -1.0, -0.29968491160191, -0.26121019262912554, 0.5129667266156068, 0.10880818515197807, 0.4379832746273899, 0.4578336686592514, 0.6410152438624929, 0.04505745487156291, 0.009320840498865012, -0.020523154351014667, 0.15769522391129656, -0.06258071892051822, -0.32998771255367343,-0.026800743748199102, 0.06143781619783842, 0.2638120455179165, -0.22637137386687595, -0.11251678451158974, 0.16285541165079287, 0.3448661028020803, 0.11032425579585989, -0.202603672236863, 0.20275884364997449, 0.021536080257503056, 0.09769910839919703, 0.04065443581198678, 0.016374516757634818, -0.04941672542098832, 0.004623650619335018, -0.03698279648337918, 0.04542971088281707, 0.029894380679541355, 0.006709929607947772, 0.1259699234656209, -0.13024251384970847, 0.011410810975428844, -0.1649990258899369, -0.05988959070358722, -0.0077795105884423486, 0.2589272069617304, 0.05920341373526867, 0.02012209142441832, 0.04875346310567188, -0.003247853249949209, -0.019200657794131867, -0.0036123823687246345, 0.0929339537345263)

df2<-data.frame(R)
df2$Combinaison <-c(replicate(length(C40C43) , "C40C43"),replicate(length(C64C63) , "C64C63"),replicate(length(C26C32) , "C26C32"),replicate(length(C21C24) , "C21C24"),replicate(length(C51C52) , "C51C52"))
#colnames(df2 ) <- c("ref","R","Experiment")

#########clustering using the elbow method to find the optimal nbr of clusters

contributingpairs<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/5SFREQUENCY_COUNTex_contributingpairsforrankedhelices.csv')
target<-c('C64-C63','C44','C22','C54','C26-C32','C36','C56','C55','C58','C27','C33','C31','C28','C38','C59','C23','C52-C51','C53','C49','C46','C45','C40-C43','C35','C68','C37','C47','C24-C21')
Rankedhelices<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/5SFREQUENCY_COUNTex_helices.csv')


contributingpairs<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/addriboswitch-apoB_FREQUENCY_COUNTex_contributingpairsforrankedhelices.csv')
target<-c('C120','C195','C59','C164-C163','C50','C64','C162','C86-C92','C60','C48-C45','C52-C43','C53-C49','C87','C72-C86','C70-C98','C186','C44','C56','C173-C182','C178-C180','C179-C178','C179-C180','C54','C174','C156','C166','C131','C58','C137-C144','C178-C177','C74','C179-C177','C88','C159','C85-C91','C169-C172','C177-C176','C145','C110','C57','C143','C71','C180-C181','C160','C76','C175','C170-C173','C138','C178-C176','C126','C112','C93','C80','C151-C154','C81-C96','C106-C101','C128','C188','C144-C147','C123','C114','C109','C99','C113','C142','C165','C67','C97-C82','C119-C116','C117-C111','C83','C75','C40','C39','C122-C121','C118-C117','C155','C132-C134','C73','C190','C103','C68-C84','C115','C133','C152','C100','C141','C69','C107-C102','C148','C90','C51','C47','C66','C62','C33','C129','C41','C149','C46','C63','C61')
Rankedhelices<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/addriboswitch-apoB_FREQUENCY_COUNTex_helices.csv')

contributingpairs<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/addriboswitch-holo_FREQUENCY_COUNTex_contributingpairsforrankedhelices.csv')
target<-c('C60','C50','C120','C64','C195','C52-C43','C59','C53-C49','C76','C87','C57','C70-C98','C86-C92','C112','C48-C45','C40','C88','C72-C86','C93','C115','C44','C164-C163','C119-C116','C143','C39','C186','C174','C170-C173','C145','C179-C180','C159','C131','C178-C180','C179-C178','C178-C177','C180-C181','C81-C96','C74','C179-C177','C142','C173-C182','C85-C91','C166','C122-C121','C106-C101','C190','C126','C160','C177-C176','C103','C151-C154','C162','C83','C137-C144','C178-C176','C169-C172','C175','C100','C97-C82','C128','C138','C188','C56','C156','C33','C110','C41','C152','C113','C107-C102','C69','C75','C114','C141','C80','C51','C123','C148','C71','C144-C147','C73','C165','C117-C111','C109','C90','C133','C68-C84','C47','C62','C99','C58','C46','C67','C132-C134','C149','C54','C118-C117','C66','C155','C129','C63','C61')
Rankedhelices<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/addriboswitch-holo_FREQUENCY_COUNTex_helices.csv')


contributingpairs<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/RMRPFREQUENCY_COUNTex_contributingpairsforrankedhelices.csv')
target<-c('C417-C418','C73-C72','C258','C71-C73','C301','C230-C232','C74','C311','C223','C298','C300','C234','C295','C71-C72','C226','C323','C265-C266','C351','C207','C262','C296','C104-C100','C312','C261-C266','C320','C222','C208-C210','C302','C225','C386-C385','C403','C208-C206','C378-C377','C81','C364-C363','C384-C385','C390-C384','C212','C218','C386-C384','C214','C362-C374','C390-C385','C336','C276','C216-C217','C188-C183','C221','C410-C411','C237','C305','C387-C385','C330','C377-C374','C284','C291-C289','C60','C281','C177','C398','C210-C206','C80','C291-C288','C299-C308','C170','C325','C290-C287','C346','C180','C187-C188','C421','C365','C133-C134','C205-C209','C386-C387','C376','C386-C388','C65','C286','C187-C183','C310-C308','C171','C95','C402-C405','C75','C247','C215','C363-C358','C405-C406','C290-C292','C405-C404','C395','C241-C240','C402-C404','C307-C309','C93','C354','C372','C407-C405','C343','C141-C145','C219','C292-C287','C289-C288','C133-C141','C55','C88','C105','C253','C352','C404-C406','C132-C134','C168-C162','C154-C162','C151-C154','C249-C248','C407-C406','C146-C148','C154-C148','C334','C84','C193-C192','C173','C128','C120','C396-C392','C151-C162','C154-C168','C133-C132','C148-C142','C304','C291-C293','C146-C142','C299-C306','C127','C151-C148','C174','C67','C64','C356-C358','C151-C146','C176','C244-C243','C132-C141','C309-C303','C399','C203-C202','C107','C114','C66','C152-C155','C102','C155-C169','C108-C109','C197','C196-C195','C155-C149','C192-C191','C397-C396','C132-C145','C317-C319','C152-C149','C339','C193-C191','C189-C201','C113','C220','C192-C202','C141-C144','C147-C143','C273','C90','C152-C163','C86-C87','C147-C149','C152-C147','C146-C139','C307-C303','C85','C191-C202','C70-C69','C169-C163','C143-C149','C96','C155-C163','C348','C256-C268','C118','C106-C110','C181-C183','C367','C196-C194','C195-C194','C130','C147-C140','C424-C423','C308-C306','C82','C310-C299','C274','C175','C275-C272','C157','C238','C187-C186','C115','C124','C242-C243','C199-C198','C387-C388','C103','C139-C142','C156-C164','C48','C101','C68','C83','C159-C158','C263','C125','C131-C123','C116','C278','C373','C122','C129-C137','C190','C381','C143-C140','C137-C142','C138-C143','C137-C139','C213','C187-C181','C62','C236','C138-C140','C135-C139','C136-C140','C112-C111','C92','C135-C137','C235','C98-C97','C135-C129','C138-C123','C285','C126','C138-C136','C131-C136','C179','C76','C89','C150-C144','C61','C338-C331','C136-C123','C370','C333-C332','C153-C144','C282','C260','C99-C97','C277-C275','C357','C144-C145','C150-C145','C245','C332-C331','C272-C271','C313','C77','C371','C277-C279','C178','C153-C150','C98-C99','C255-C254','C166-C160','C326','C165-C160','C270-C257','C150-C161','C360','C270-C269','C293-C288','C119','C167-C165','C186-C181','C166-C167','C255-C267','C166-C161','C166-C165','C359','C59','C153-C161','C246-C264','C51','C167-C161','C250','C322','C153-C167','C63','C52','C341-C340','C57','C267-C271','C335','C280','C56','C58','C347','C324','C337','C315-C314','C318')
Rankedhelices<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/RMRPFREQUENCY_COUNTex_helices.csv')


contributingpairs<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/RMRPincelFREQUENCY_COUNTex_contributingpairsforrankedhelices.csv')
target<-c('C417-C418','C258','C323','C73-C72','C71-C73','C311','C205-C209','C214','C48','C71-C72','C290-C292','C226','C403','C290-C287','C292-C287','C399','C301','C300','C330','C253','C55','C317-C319','C178','C179','C296','C223','C61','C60','C326','C212','C381','C230-C232','C89','C52','C291-C288','C177','C291-C289','C265-C266','C357','C398','C75','C351','C56','C51','C244-C243','C92','C68','C424-C423','C346','C320','C421','C370','C295','C410-C411','C122','C397-C396','C57','C359','C309-C303','C215','C289-C288','C307-C303','C85','C261-C266','C305','C237','C312','C291-C293','C174','C364-C363','C360','C334','C276','C208-C206','C74','C241-C240','C407-C406','C304','C395','C352','C242-C243','C65','C199-C198','C402-C405','C307-C309','C218','C333-C332','C210-C206','C325','C90','C371','C67','C208-C210','C256-C268','C407-C405','C402-C404','C59','C378-C377','C310-C308','C127','C260','C81','C180','C373','C76','C64','C124','C249-C248','C298','C293-C288','C262','C119','C377-C374','C135-C129','C235','C339','C367','C197','C222','C83','C133-C134','C137-C139','C106-C110','C362-C374','C70-C69','C128','C129-C137','C80','C332-C331','C277-C279','C299-C308','C104-C100','C356-C358','C236','C135-C139','C135-C137','C396-C392','C151-C148','C277-C275','C137-C142','C154-C162','C363-C358','C82','C154-C168','C139-C142','C146-C139','C168-C162','C118','C116','C151-C162','C88','C247','C348','C308-C306','C313','C151-C146','C390-C384','C84','C151-C154','C405-C406','C284','C115','C132-C134','C146-C148','C96','C63','C322','C338-C331','C154-C148','C275-C272','C365','C103','C282','C405-C404','C133-C132','C107','C77','C66','C62','C278','C404-C406','C299-C306','C133-C141','C246-C264','C285','C93','C148-C142','C203-C202','C195-C194','C384-C385','C101','C302','C216-C217','C171','C219','C159-C158','C113','C193-C191','C191-C202','C190','C189-C201','C220','C267-C271','C196-C195','C386-C384','C354','C270-C257','C192-C191','C196-C194','C132-C141','C270-C269','C192-C202','C125','C141-C145','C147-C149','C126','C193-C192','C286','C152-C147','C86-C87','C187-C183','C152-C149','C146-C142','C310-C299','C155-C169','C143-C140','C132-C145','C156-C164','C188-C183','C157','C147-C140','C114','C221','C143-C149','C170','C213','C274','C372','C108-C109','C141-C144','C386-C387','C187-C188','C181-C183','C138-C123','C387-C385','C58','C376','C98-C99','C390-C385','C155-C149','C147-C143','C187-C181','C131-C123','C152-C155','C99-C97','C273','C130','C95','C225','C136-C123','C169-C163','C234','C341-C340','C245','C186-C181','C152-C163','C98-C97','C386-C385','C155-C163','C187-C186','C131-C136','C105','C387-C388','C272-C271','C138-C143','C138-C140','C136-C140','C238','C386-C388','C207','C173','C138-C136','C102','C153-C144','C255-C267','C144-C145','C150-C144','C153-C161','C336','C150-C145','C343','C165-C160','C153-C150','C112-C111','C153-C167','C150-C161','C166-C161','C167-C165','C166-C160','C255-C254','C166-C167','C175','C166-C165','C250','C167-C161','C335','C318','C347','C263','C280','C120','C315-C314','C281','C176','C324','C337')
Rankedhelices<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/RMRPincelFREQUENCY_COUNTex_helices.csv')


contributingpairs<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/Nuy77FREQUENCY_COUNTex_contributingpairsforrankedhelices.csv')
target<-c('C10-C7','C33','C14','C17','C21','C28','C13','C18-C20','C8-C3','C25','C38','C37','C22','C29','C12','C30','C36','C11','C31','C1-C6','C35','C32','C23','C19-C15','C24','C16-C15','C4','C9','C26-C27','C34','C5')
Rankedhelices<- read.csv(file = '/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2021-publication-withoutfiltering_Feb/Nuy77FREQUENCY_COUNTex_helices.csv')


M<-mean(Rankedhelices$averageRLtdiff)
sigma<-sd(Rankedhelices$averageRLtdiff)
          
library(magrittr)
library(dplyr)
library(tidyverse)
Rankedhelices <-Rankedhelices  %>%
  mutate(Statt= case_when(Rankedhelices$averageRLtdiff >=M+2*sigma  ~ '>= M+2sigma',
                          Rankedhelices$averageRLtdiff >=M+1*sigma & Rankedhelices$averageRLtdiff <M+2*sigma ~ '>=M+sigma & <M+2sigma',
                          Rankedhelices$averageRLtdiff >=M & Rankedhelices$averageRLtdiff <M+1*sigma ~ '>=M & <M+sigma',
                          Rankedhelices$averageRLtdiff >=M-sigma & Rankedhelices$averageRLtdiff <M~ ' >=M-sigma& <M',
                          Rankedhelices$averageRLtdiff <M-sigma ~ ' <M-sigma& <M' ))
#install.packages("tidyverse")       


set.seed(123)
library(purrr)
# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(Rankedhelices$averageRLtdiff, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)
N<-0.1*max(wss_values)
thresh<-1+sum(  wss_values >=N )
thresh
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

k <-kmeans(Rankedhelices$averageRLtdiff, centers=thresh) #Create 5 clusters, Remove columns 1 and 2
#k$centers #Display&nbsp;cluster centers
#table(k$cluster) #Give a count of data points in each cluster
#k$withinss
k$cluster
# xlim(5, 20)+.
Newdata<-data.frame(target,k$cluster,Rankedhelices$Statt)
Distributionlengths<- read.csv(file ='/home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/Feb2021_profiling/Nuy77Profiling.csv')
#write.table(Newdata,"/home/user/Documents/Alain/finalfigures/RMRPincellclustering")
draw<-merge(Distributionlengths, Newdata, by.x="Helix", by.y="target")
draw[order(draw$k.cluster),]
contributingpairs<-merge(contributingpairs, Newdata, by.x="Helix", by.y="target")
contributingpairs$Helix <- factor(contributingpairs$Helix, levels=unique(target))

setEPS()
postscript("/home/user/Documents/Alain/finalfigures/clusterNuy77combinaisonhelices.eps", width =30.0, height = 6.0)



library(ggplot2)
p<-ggplot(contributingpairs,aes(x=1,y=Rforpairs,color=as.factor(k.cluster)))+#as.factor(Helix))) +
  geom_violin(colour = "grey50")+
  geom_boxplot(width=0.1,colour = "grey50")+
  stat_summary(fun.y=mean, geom="point", size=2, color="black")+
  geom_jitter( size = 1.5)+ylim(-1,3)+
  
  facet_wrap(~as.factor(contributingpairs$Helix),nrow = 1)

p<-p+ theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x=element_blank())
#p<-p+ geom_hline(yintercept=M, linetype="dashed", color = "red")+ geom_hline(yintercept=M+sigma, linetype="dashed", color = "grey")+ geom_hline(yintercept=M-sigma, linetype="dashed", color = "grey")+ geom_hline(yintercept=M+2*sigma, linetype="dashed", color = "grey")
p
dev.off()
#############classification baased on sigma and mean of all Rc values
setEPS()
postscript("/home/user/Documents/Alain/finalfigures/Mean_sigmaclusterNuy77combinaisonhelices.eps", width =30.0, height = 6.0)


g<-ggplot(contributingpairs,aes(x=1,y=Rforpairs,color=as.factor(Rankedhelices.Statt)))+#as.factor(Helix))) +
  geom_violin(colour = "grey50")+
  geom_boxplot(width=0.1,colour = "grey50")+
  stat_summary(fun.y=mean, geom="point", size=2, color="black")+
  geom_jitter( size = 1.5)+ylim(-1,3)+
  
  facet_wrap(~as.factor(contributingpairs$Helix),nrow = 1)

g<-g+ theme(legend.title= element_blank(),legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x=element_blank())
#p<-p+ geom_hline(yintercept=M, linetype="dashed", color = "red")+ geom_hline(yintercept=M+sigma, linetype="dashed", color = "grey")+ geom_hline(yintercept=M-sigma, linetype="dashed", color = "grey")+ geom_hline(yintercept=M+2*sigma, linetype="dashed", color = "grey")
g

dev.off()

#install.packages("Rcpp") 
