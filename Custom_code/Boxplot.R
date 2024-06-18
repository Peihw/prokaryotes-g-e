#Boxplot
library(openxlsx)
dup80=dup[!is.na(dup$`Redundancy.ration(Set-80)`),]
phy=as.data.frame(table(dup80$Phylum))
phya=as.data.frame(table(dup$Phylum))
phy$all=phya$Freq[match(phy$Var1,phya$Var1)]
phy$ratio=round(phy$Freq/phy$all*100,1)
phy$name80=paste(phy$Var1,"(n=",phy$Freq,", ",phy$ratio,"%)",sep = "")
phy$Var1=factor(phy$Var1,levels = phylist)
phy=na.omit(phy)
dup80$name80=phy$name80[match(dup80$Phylum,phy$Var1)]
phy=phy[order(phy$Var1),]
dup80$name80=factor(dup80$name80,levels = phy$name80)
library(ggplot2)
ggplot(dup80[!is.na(dup80$name80),], aes(x=name80,y=`Redundancy.ration(Set-80)`,fill=name80)) + 
       geom_boxplot(outlier.shape = 1,outlier.size = 1) + 
       scale_y_continuous(trans = log10_trans(),expand=c(0,0))+
       expand_limits(y=c(0,1))+theme_classic()+
       theme(axis.text.x  = element_text(angle=60, hjust=1),axis.text = element_text(size=9))+ 
       labs(x = "",y="Redundancy ratio")+scale_color_gradientn(colours = rainbow(32)) + 
       guides(fill=FALSE)+stat_summary(fun=mean, geom="point", aes(group=name80),color="brown", size=1,shape=6)


