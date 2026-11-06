#directory

#Cage Survival Analysis
#clearing workspace
rm(list=ls())
#unloading any packages
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=FALSE))

#libraries
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggExtra)
library(lvmisc)

cage<-read.csv("Cage23.csv")
cage$Treatment<-cage$Treatment%>%ordered(levels=c("Control","Amitraz"))


#proportion
cage<-cage%>%mutate(n.Caged=as.numeric(n.Caged),
                    n.Survived=as.numeric(n.Survived),
                    p.Surv=n.Survived/n.Caged,
                    asin.Surv=asin(sqrt(p.Surv)),
                    Panel=ifelse(Treatment=="Control","Control","Treatment")
                    )
#descriptives
cage$n.Caged%>%range(na.rm=TRUE)
cage$n.Caged%>%median(na.rm=TRUE)
cage%>%group_by(Colony,Treatment)%>%summarise(N=n())%>%print(n=36)

cage%>%
  group_by(Treatment)%>%
  summarise(Mean=mean(p.Surv,na.rm=TRUE),
            N=n()-sum(is.na(p.Surv)),
            SE=sd(p.Surv,na.rm=TRUE)/sqrt(N),
            Med=median(p.Surv,na.rm=TRUE),
            p.25=quantile(p.Surv,na.rm=TRUE,probs=.25),
            p.75=quantile(p.Surv,na.rm=TRUE,probs=.75),
            IQR=p.75-p.25)

hist(cage$p.Surv)
hist(cage$asin.Surv)
 
# test for a treatment effect on the proportion of survived.
fit1<-cage%>%lm(asin.Surv~Treatment*n.Caged,data=.)
anova(fit1)
confint(fit1)
plot_model(fit1)
#initial model revealed no significances.
library(lme4)
library(lmerTest)
fit2<-cage%>%lmer(asin.Surv~Treatment*n.Caged+(1|Colony),data=.)
anova(fit2)
plot_model(fit2)
#using the random effects model doesn't change anything so we'll leave it go.

#vis
p1a<-
ggplot(cage,aes(y=p.Surv,x=Treatment))+
  geom_boxplot(aes(fill=Treatment),alpha=.6,outlier.alpha=0)+
  geom_point(aes(color=Treatment),position=position_jitterdodge())+
  theme_minimal()+
  theme(axis.title.y=element_text(size=20,face="bold"),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14,face="bold",angle=-30,hjust=0),
        axis.title.x=element_blank(),
        legend.position="none",
        plot.margin=margin(.25,1,.25,.25,unit="cm"))+
  scale_fill_manual(values=c("forestgreen","firebrick"))+
  scale_color_manual(values=c("forestgreen","firebrick"))+
  labs(y="Proportion drones survived")

p1a
ggsave(plot=p1a,filename="F1_Cage23.png", height= 7, width=7, units="in")



#Immature reproductive & physiological quality
rm(list=ls())
#unloading any packages
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=FALSE))
#libraries

library(tidyverse)
library(readxl)
library(broom)
library(RColorBrewer)
library(ggpubr)

Immatures<-read_csv("Immatures23.csv")
Immatures<-Immatures%>%mutate(Colony=as.character(Colony),
                              Treatment=ordered(Treatment,levels=c("Control","Amitraz")))

#linear models
Fixed.Immatures<-
  Immatures%>%
  pivot_longer(cols=c(4:7),names_to="Parameter",values_to="Value")%>%
  group_by(Parameter,Year)%>%
  nest()%>%
  mutate(models=lapply(data,function(df)anova(lm(Value~Treatment+Colony,data=df))),results=map(models,tidy))%>%
  unnest(results)%>%
  subset(select=-c(3,4,7,8))%>%
  pivot_wider(names_from=term,values_from=c(4:6))%>%
  subset(select=-c(8,11))
Fixed.Immatures
#what we see is that colony Is a highly significant factor across all parameters and years. There are differences between the two treatments as well, but I would think that a random effect model is more appropriate for what we’re doing here.

#descriptives
Immatures%>%
  group_by(Treatment,Colony)%>%
  summarise(N=n())%>%
  group_by(.)%>%
  na.omit()%>%
  summarise(Min=min(N),Max=max(N),Mean=mean(N),Med=median(N))

Des.Immatures<- Immatures%>%
  pivot_longer(c(4:7),names_to="Parameter",values_to="Value")%>%
  group_by(Parameter,Year,Treatment)%>%
  summarise(N=n()-sum(is.na(Value==TRUE)),
            Mean=mean(Value,na.rm=TRUE),
            se=sd(Value,na.rm=TRUE)/sqrt(N),
            Med=median(Value,na.rm=TRUE),
            Quant55=quantile(Value,probs=.55,na.rm=TRUE),
            Quant95=quantile(Value,probs=.95,na.rm=TRUE))%>%
  mutate(Panel=ifelse(Treatment=="Control","Control","Treatment"))
Des.Immatures
#colony is significant in several of the parameters, but is not inherently interesting for our analyses (of course drones are likely to differ by colony) and while this might become more useful in studying the multivariate data, for these analyses, we just want to control for it. We’ll therefore assess these models using mixed effects with Colony as a random factor.

library(lme4)
library(car)


Rand.Immatures<-Immatures%>%
  pivot_longer(c(4:7),names_to="Parameter",values_to="Value")%>%
  group_by(Parameter,Year)%>%
  nest()%>%
  mutate(models=lapply(data,function(df)Anova(lmer(Value~Treatment+(1|Colony),data=df))),results=map(models,tidy))%>%
  unnest(results)%>%
  subset(select=-c(3,4))

Rand.Immatures<-Des.Immatures%>%
  group_by(Parameter,Year)%>%
  summarise(Quant95=mean(Quant95))%>%
  full_join(Rand.Immatures,.)%>%
  group_by()%>%
  mutate(p.rank=rank(p.value, ties.method="last"),
         p.crit=(p.rank/n()*0.1),
         sig=ifelse(p.value<=p.crit,"*","NS"))
Rand.Immatures

#including colony as a random effect eliminates the significant impact of treatment from the main-effects fixed models above. I interpret this to mean that colony variation is high enough that the trends seen in the treatment are swamped out.
detach("package:car",unload=FALSE)

#figure
l.Immatures<-Immatures%>%
  pivot_longer(cols=c(4:7),names_to="Parameter",values_to="Value")%>%
  mutate(Panel=ifelse(Treatment=="Control","Control","Treatment"))
immatures.labs<-c("Body (mg)","Testes (mg)","Head (mm)","Thorax (mm)")
names(immatures.labs)<-c("m.Body","m.Testes","w.Head","w.Thorax")

D1<- l.Immatures%>%ggplot(aes(x=Treatment,y=Value))+
  geom_violin(aes(fill=Treatment),alpha=.4)+
  geom_boxplot(aes(fill=Treatment),alpha=.4)+
  geom_text(data=Des.Immatures,aes(x=Treatment,y=Quant55,label=paste(sprintf("%0.2f",round(Mean,2)),sprintf("%0.2f",round(se,2)),sep="±")),vjust=0,size=2)+
  geom_text(data=Rand.Immatures,aes(x=1.5,y=Quant95,label=paste("P",sprintf("%0.2f",round(p.value,2)),sep="=")),vjust=0,size=3)+
  facet_wrap(facets=vars(Parameter),scales="free",labeller=labeller(Parameter=immatures.labs))+
  theme_minimal()+
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=14,face="bold"),
        strip.background=element_rect(fill="slategray2",color=NA),
        strip.text=element_text(size=12,face="bold"),
        legend.position="none",
        axis.title=element_blank())+
  scale_fill_manual(values=c("forestgreen","firebrick"))+
  scale_color_manual(values=c("forestgreen","firebrick"))
D1
ggsave("Immature23.png",height=7,width=7,unit="in")
ggsave("Immature23.pdf",height=7,width=7,unit="in", device=pdf)

#Summary Table
full_join(
Des.Immatures%>%subset(select=c(1:6)),
Rand.Immatures%>%subset(select=c(1,2,4,6)))%>%
  write_csv("Immatures_Stats.csv")

#Sample selection
Sel.Immatures<-Immatures%>%
  group_by(Treatment,Colony)%>%
  summarise(Testes=mean(m.Testes,na.rm=TRUE))%>%
  group_by(Treatment)%>%mutate(Rank=rank(Testes))%>%
  filter(Rank>=max(Rank)-1|Rank<=2)%>%
  mutate(pRank=ifelse(Rank>2,"Best","Worst"))%>%
  group_by(pRank,Treatment,Colony)%>%
  summarise(Testes=mean(Testes))%>%
  group_by(pRank,Treatment)%>%
  arrange(desc(Testes),.by_group = TRUE)

#figure
target<-Sel.Immatures$Colony%>%paste(collapse = "|")

Fil.Immatures<-Immatures%>%
  filter(str_detect(Colony,target))%>%
  left_join(.,Sel.Immatures)%>%
  na.omit()%>%
  mutate(Colony=factor(Colony,levels=Sel.Immatures$Colony))

Fil.Immatures%>%
  ggplot(aes(x=Colony,y=m.Testes))+
  geom_boxplot(aes(fill=Treatment),alpha=.4)+
  geom_point(aes(color=Treatment),position=position_dodge2(width=.3))+
  facet_wrap(facets=vars(pRank),scale="free_x",strip.position="bottom")+
  theme_minimal()+
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=14,face="bold"),
        strip.background=element_rect(fill="white",color="black"),
        strip.text=element_text(size=12,face="bold"),
        legend.position="inside",
        legend.position.inside = c(.75,.8),
        legend.text = element_text(size=12,face="bold"),
        legend.title=element_blank(),
        axis.title=element_text(size=14,face="bold"))+
  scale_fill_manual(values=c("forestgreen","firebrick"))+
  scale_color_manual(values=c("forestgreen","firebrick"))+
  labs(y="Testes mass (mg)")
ggsave("Sel_Immatures.png")
ggsave("Sel_Immatures.pdf")

#Mature drone analysis
#clearing prior data
rm(list=ls())
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=FALSE))

library(tidyverse)
library(readxl)
library(broom)
library(RColorBrewer)
library(ggpubr)

Matures<-read_csv("Mature23.csv")
Matures<-Matures%>%
  mutate(Treatment=ordered(recode(Treatment,"Treatment"="Amitraz"),levels=c("Control","Amitraz")),
         av.Sperm=asin(sqrt(v.Sperm)),
         lt.Sperm=log(t.Sperm+1),
         ll.Sperm=log(l.Sperm+1),
         ld.Sperm=log(d.Sperm+1),
         Colony=as.character(Colony)
  )%>%
  filter(is.na(Treatment)==FALSE)

#Preliminary Figure
#Sperm traits
Matures%>%
  subset(select=c(Treatment,Colony, av.Sperm,lt.Sperm,ll.Sperm,ld.Sperm))%>%
  pivot_longer(c(3:6),names_to="Trait",values_to='value')%>%
  ggplot(aes(x=Treatment,y=value))+
  geom_violin(aes(fill=Treatment),alpha=.4)+
  geom_boxplot(aes(fill=Treatment),alpha=.4)+
  facet_wrap(facets=vars(Trait),scales="free_y")

#Morpho-repro traits
Matures%>%
  subset(select=c(Treatment,Colony, m.Body,w.Head,w.Thorax,l.SVs,l.MGs))%>%
  pivot_longer(c(3:7),names_to="Trait",values_to='value')%>%
  ggplot(aes(x=Treatment,y=value))+
  geom_violin(aes(fill=Treatment),alpha=.4)+
  geom_boxplot(aes(fill=Treatment),alpha=.4)+
  facet_wrap(facets=vars(Trait),scales="free_y")


#not very impressive; right now, it appears there’s little difference in the matures
Fixed.Mature<- 
  Matures%>%
  subset(select=c(Colony,Treatment,m.Body,w.Thorax,w.Head,l.SVs,l.MGs,av.Sperm,ll.Sperm,ld.Sperm,lt.Sperm))%>%
  pivot_longer(cols=c(3:ncol(.)),names_to="Parameter",values_to="Value")%>%
  group_by(Parameter)%>%
  nest()%>%
  mutate(models=lapply(data,function(df)anova(lm(Value~Treatment+Colony,data=df))),results=map(models,tidy))%>%
  unnest(results)%>%
  subset(select=-c(2,3,6,7))%>%
  pivot_wider(names_from=term,values_from=c(3:5))%>%
  subset(select=-c(7,10))%>%
  group_by()
Fixed.Mature
#Fixed effects show extreme effect of colony as well as a few treatment effects.

Des.Mature<- 
  Matures%>%
  subset(select=c(Colony,Treatment,m.Body,w.Thorax,w.Head,l.SVs,l.MGs,av.Sperm,ll.Sperm,ld.Sperm,lt.Sperm))%>%
  pivot_longer(cols=c(3:ncol(.)),names_to="Parameter",values_to="Value")%>%
  group_by(Parameter,Treatment)%>%
  summarise(N=n()-sum(is.na(Value==TRUE)),
            Mean=mean(Value,na.rm=TRUE),
            se=sd(Value,na.rm=TRUE)/sqrt(N),
            Med=median(Value,na.rm=TRUE),
            Quant55=quantile(Value,probs=.55,na.rm=TRUE),
            Quant95=quantile(Value,probs=.95,na.rm=TRUE))%>%
  mutate(Panel=ifelse(Treatment=="Control","Control","Treatment"))
Des.Mature

#Complete analysis (of the non-factor traits)
library(lme4)
library(car)
Rand.Mature<- 
  Matures%>%
  subset(select=c(Colony,Treatment,m.Body,w.Thorax,w.Head,l.SVs,l.MGs,av.Sperm,ll.Sperm,ld.Sperm,lt.Sperm))%>%
  pivot_longer(cols=c(3:ncol(.)),names_to="Parameter",values_to="Value")%>%
  group_by(Parameter)%>%
  nest()%>%
  mutate(models=lapply(data,function(df)Anova(lmer(Value~Treatment+(1|Colony),data=df))),results=map(models,tidy))%>%
  unnest(results)%>%
  subset(select=-c(2,3,4))%>%
  group_by()%>%
  mutate(p.rank=rank(p.value, ties.method="last"),
         p.crit=(p.rank/n()*0.1),
         Sig=ifelse(p.value<=p.crit,"SIG","NS"))

Rand.Mature

#and what we find is that including colony as a random factor in the models eliminates any significances that we see in the data using only the fixed factors. My interpretation of this is that the main-effects, fixed models show that there’s a non-independent impact of colony and treatment on drone parameters, but that these differences are driven by between-colony variation rather than between-treatment variation overall.
detach("package:car",unload=FALSE)

#figure
#Labels
Labels<-
  full_join(Des.Mature,Rand.Mature)%>%
  full_join(.,
            tibble(Parameter=unique(.$Parameter),
                   Label=c("Sperm viability (asin-%)",
                           "Mucus glands (mm)",
                           "Seminal vesicles (mm)",
                           "Dead sperm (ln)",
                           "Live sperm (ln)",
                           "Total sperm (ln)",
                           "Body (mg)",
                           "Head (mm)",
                           "Thorax (mm)"),
                   Section=c("Sperm","Repro","Repro","Sperm","Sperm","Sperm","Morpho","Morpho","Morpho")))
Labels

#Pt 1: morphological traits: box & violin plot
m.Matures<-Matures%>%
  subset(select=c(Colony,Treatment,m.Body,w.Thorax,w.Head))%>%
  pivot_longer(c(3:5),names_to="Parameter")%>%
  na.omit()
m.Labels<-Labels%>%
  filter(Section=="Morpho")
m.Titles<-unique(m.Labels$Label)
names(m.Titles)<-unique(m.Labels$Parameter)
m.Sig<-m.Labels%>%
  group_by(Parameter,p.value,Sig)%>%
  summarise(Quant95=mean(Quant95))%>%
  group_by()

P1<-ggplot(data=m.Matures,aes(x=Treatment,y=value))+
  geom_violin(aes(fill=Treatment),alpha=.4)+
  geom_boxplot(aes(fill=Treatment),alpha=.4)+
  geom_text(data=m.Labels,
            aes(x=Treatment,y=Quant55,
                label=paste(sprintf("%0.2f",round(Mean,2)),
                            sprintf("%0.2f",round(se,2)),sep="±")),
            vjust=0,size=3)+
  geom_text(data=m.Sig,
            aes(x=1.5,y=Quant95,
                label=paste("P",sprintf("%0.2f",round(p.value,2)),sep="=")),
            vjust=0,size=4)+
  facet_wrap(facets=vars(Parameter),scales="free_y",labeller=labeller(Parameter=m.Titles))+
  theme_minimal()+
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=14,face="bold"),
        strip.background=element_rect(fill="slategray2",color=NA),
        strip.text=element_text(size=12,face="bold"),
        plot.title=element_text(size=22,face="bold",hjust=.5),
        axis.title=element_blank(),
        legend.position = "none")+
  scale_fill_manual(values=c("forestgreen","firebrick"))+
  scale_color_manual(values=c("forestgreen","firebrick"))

#Pt 2: Sperm Traits
s.Matures<-Matures%>%
  subset(select=c(Colony,Treatment,av.Sperm,ll.Sperm,ld.Sperm,lt.Sperm))%>%
  pivot_longer(c(3:6),names_to="Parameter")%>%
  na.omit()
s.Labels<-Labels%>%
  filter(Section=="Sperm")
s.Titles<-unique(s.Labels$Label)
names(s.Titles)<-unique(s.Labels$Parameter)
s.Sig<-s.Labels%>%
  group_by(Parameter,p.value,Sig)%>%
  summarise(Quant95=mean(Quant95))%>%
  group_by()

P2<-ggplot(data=s.Matures,aes(x=Treatment,y=value))+
  geom_violin(aes(fill=Treatment),alpha=.4)+
  geom_boxplot(aes(fill=Treatment),alpha=.4)+
  geom_text(data=s.Labels,
            aes(x=Treatment,y=Quant55,
                label=paste(sprintf("%0.2f",round(Mean,2)),
                            sprintf("%0.2f",round(se,2)),sep="±")),
            vjust=0,size=3)+
  geom_text(data=s.Sig,
            aes(x=1.5,y=Quant95,
                label=paste("P",sprintf("%0.2f",round(p.value,2)),sep="=")),
            vjust=0,size=4)+
  facet_wrap(facets=vars(Parameter),scales="free_y",labeller=labeller(Parameter=s.Titles))+
  theme_minimal()+
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=14,face="bold"),
        strip.background=element_rect(fill="slategray2",color=NA),
        strip.text=element_text(size=12,face="bold"),
        plot.title=element_text(size=22,face="bold",hjust=.5),
        axis.title=element_blank(),
        legend.position = "none")+
  scale_fill_manual(values=c("forestgreen","firebrick"))+
  scale_color_manual(values=c("forestgreen","firebrick"))

#Pt 3: Reproductive gland traits
r.Matures<-Matures%>%
  subset(select=c(Colony,Treatment,l.SVs,l.MGs))%>%
  pivot_longer(c(3:4),names_to="Parameter")%>%
  na.omit()
r.Labels<-Labels%>%
  filter(Section=="Repro")
r.Titles<-unique(r.Labels$Label)
names(r.Titles)<-unique(r.Labels$Parameter)
r.Sig<-r.Labels%>%
  group_by(Parameter,p.value,Sig)%>%
  summarise(Quant95=mean(Quant95))%>%
  group_by()

P3<-ggplot(data=r.Matures,aes(x=Treatment,y=value))+
  geom_violin(aes(fill=Treatment),alpha=.4)+
  geom_boxplot(aes(fill=Treatment),alpha=.4)+
  geom_text(data=r.Labels,
            aes(x=Treatment,y=Quant55,
                label=paste(sprintf("%0.2f",round(Mean,2)),
                            sprintf("%0.2f",round(se,2)),sep="±")),
            vjust=0,size=3)+
  geom_text(data=r.Sig,
            aes(x=1.5,y=Quant95,
                label=paste("P",sprintf("%0.2f",round(p.value,2)),sep="=")),
            vjust=0,size=4)+
  facet_wrap(facets=vars(Parameter),scales="free_y",
             labeller=labeller(Parameter=r.Titles))+
  theme_minimal()+
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=14,face="bold"),
        strip.background=element_rect(fill="slategray2",color=NA),
        strip.text=element_text(size=12,face="bold"),
        plot.title=element_text(size=22,face="bold",hjust=.5),
        axis.title=element_blank(),
        legend.position = "none")+
  scale_fill_manual(values=c("forestgreen","firebrick"))+
  scale_color_manual(values=c("forestgreen","firebrick")) 


ggarrange(ggarrange(P1,P3,labels=c("A","C"),ncol = 1),P2,labels=c("","B"))
ggsave("Matures23.png",height=14,width=14,units="in")
ggsave("Matures23.pdf",height=14,width=14,units="in",device="pdf")

#summary table
Matures%>%
  subset(select=c(Colony,Treatment,m.Body,w.Head,w.Thorax,l.SVs,l.MGs,t.Sperm,v.Sperm))%>%
  group_by(Treatment,Colony)%>%
  na.omit()%>%
  summarise(N=n())%>%
  group_by()%>%
  summarise(Max=max(N),Min=min(N),Med=median(N))

#Summary Table
Matures%>%subset(select=c(Treatment,m.Body,w.Head,w.Thorax,l.SVs,l.MGs,t.Sperm,l.Sperm,d.Sperm,v.Sperm))%>%
  pivot_longer(cols=c(2:ncol(.)),names_to="Parameter",values_to="Value")%>%
  group_by(Treatment,Parameter)%>%
  summarise(N=n()-sum(is.na(Value)==TRUE),Mean=mean(Value,na.rm=TRUE),SE=sd(Value,na.rm=TRUE)/sqrt(N))%>%
  write_csv("Mature_Descriptives.csv")

Rand.Mature%>%
  write_csv("Matures_Stats.csv")

#Sample selection START HERE
#Sample selection
Sel.Matures<-Matures%>%
  group_by(Treatment,Colony)%>%
  summarise(Sperm=mean(l.Sperm,na.rm=TRUE))%>%
  group_by(Treatment)%>%mutate(Rank=rank(Sperm))%>%
  filter(Rank>=max(Rank)-1|Rank<=2)%>%
  mutate(pRank=ifelse(Rank>2,"Best","Worst"))%>%
  group_by(pRank,Treatment,Colony)%>%
  summarise(Sperm=mean(Sperm))%>%
  group_by(pRank,Treatment)%>%
  arrange(desc(Sperm),.by_group = TRUE)%>%
  filter(Colony!="8")

#figure
target<-Sel.Matures$Colony%>%paste(collapse = "|")

Fil.Matures<-Matures%>%
  filter(str_detect(Colony,target))%>%
  left_join(.,Sel.Matures)%>%
  filter(is.na(Ejaculated)==TRUE)%>%
  subset(select=c(Colony,Treatment,pRank,v.Sperm,t.Sperm,l.Sperm,d.Sperm))%>%
  mutate(Colony=factor(Colony,levels=Sel.Matures$Colony))

Fil.Matures%>%
  ggplot(aes(x=Colony,y=l.Sperm))+
  geom_boxplot(aes(fill=Treatment),alpha=.4,outliers = FALSE)+
  geom_point(aes(color=Treatment),position=position_dodge2(width=.3))+
  facet_wrap(facets=vars(pRank),scale="free_x",strip.position="bottom")+
  theme_minimal()+
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=14,face="bold"),
        strip.background=element_rect(fill="white",color="black"),
        strip.text=element_text(size=12,face="bold"),
        legend.position="inside",
        legend.position.inside = c(.75,.8),
        legend.text = element_text(size=12,face="bold"),
        legend.title=element_blank(),
        axis.title=element_text(size=14,face="bold"))+
  scale_fill_manual(values=c("forestgreen","firebrick"))+
  scale_color_manual(values=c("forestgreen","firebrick"))+
  labs(y="Live Sperm")
ggsave("Sel_Matures.png")
ggsave("Sel_Matures.pdf")


#Colony ranking for selection of microbial samples
#Mature
m.Rank<-
  Matures%>%group_by(Treatment,Colony)%>%summarise(Valid.Samples=n()-sum(is.na(l.Sperm)==TRUE),Sperm=mean(l.Sperm,na.rm=TRUE)) %>%group_by(Treatment)%>%mutate(mRank=rank(Sperm))

m.Select<-left_join(m.Rank,Matures)

Hoops.Codes<-read_excel("2023-BaC-Consolidated.xlsx",sheet="2022-Testeepro Mature")%>%subset(select=c(1,4)) %>%mutate(Code="Hoops")
names(Hoops.Codes)<- c("QD","Code","GN")

Trick.Codes<-read_excel("2023-BaC-Consolidated.xlsx",sheet="2023-Amitrestes Mature")%>%subset(select=c(2,5))%>%mutate(Code="Trick")
names(Trick.Codes)<-c("QD","Code","GN")
Codes<-rbind(Hoops.Codes,Trick.Codes)

mLexi<-left_join(mLexi,Codes)

#Immature
i.Rank<-Immatures%>%group_by(Year,Treatment,Colony)%>%summarise(Valid.Samples=n()-sum(is.na(m.Testes)==TRUE),Testes=mean(m.Testes,na.rm=TRUE))%>% group_by(Year,Treatment)%>%mutate(iRank=rank(Testes))


iLexi<-left_join(i.Rank,Immatures)%>%filter(is.na(m.Testes)==FALSE)

write_csv(mLexi,"mLexi.csv")
write_csv(iLexi,"iLexi.csv")


##Pathogens & Parasites
#Mite count analyses

#clear any BS (if starting from an open workspace)
rm(list=ls())
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=FALSE))

library(tidyverse)
library(readxl)
library(broom)
library(RColorBrewer)
library(ggpubr)
library(rstatix)

Mites<-read.csv("Mites23.csv")
Mites<-Mites%>%filter(Period!="Post-Formic")


#Ordering
Mites$Period<-ordered(Mites$Period,levels=c("Start of Study","Mid Study","End of Study","Post-Formic"))

Mites<-Mites%>%mutate(asin.Mites=asin(sqrt(p.Mites)))
Mites<-Mites%>%mutate(Treatment=ifelse(Mites$Year=="2022",ifelse(Mites$Colony<21,"Terrapro","Control"),ifelse(Mites$Colony<10|(Mites$Colony>20&Mites$Colony<31),"Control","Amitraz")))
Mites$Treatment<-ordered(Mites$Treatment,levels=c("Control","Amitraz"))

Mites<-Mites%>%filter(is.na(p.Mites)!=TRUE)

#Mite figure (no stats)
p2<-ggplot(Mites,aes(x=Period,y=p.Mites))+
  geom_boxplot(aes(fill=Treatment),alpha=.4,outlier.alpha=0)+
  geom_point(aes(color=Treatment),position=position_jitterdodge())+
  theme_minimal()+
  labs(x=NULL,y="Mite Proportion")+
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14,face="bold",angle=-30,hjust=0),
        strip.background=element_rect(fill="slategray2",color=NA),
        strip.text=element_text(size=14,face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=12,face="bold"),
        legend.position="inside",
        legend.position.inside=c(.25,.7),
        plot.margin=margin(.25,1,.25,.25,unit="cm"))+
  scale_fill_manual(values=c("forestgreen","firebrick"))+
  scale_color_manual(values=c("forestgreen","firebrick"))



#summary table
Mites%>%
  summarise(N=n()-sum(is.na(p.Mites)),
            Mean=mean(p.Mites,na.rm=TRUE),
            SE=sd(p.Mites,na.rm=TRUE)/sqrt(N))

Mites%>%group_by(Period)%>%
  summarise(N=n()-sum(is.na(p.Mites)),
            Mean=mean(p.Mites,na.rm=TRUE),
            SE=sd(p.Mites,na.rm=TRUE)/sqrt(N))

Mites%>%group_by(Treatment,Period)%>%
  summarise(N=n()-sum(is.na(p.Mites)),
            Mean=mean(p.Mites,na.rm=TRUE),
            SE=sd(p.Mites,na.rm=TRUE)/sqrt(N))

#a rise in mites over time, but highly variable
Mites%>%group_by(Period,Treatment)%>%identify_outliers(asin.Mites)
#no extreme outliers

ggqqplot(Mites,"asin.Mites",facet.by=c("Treatment","Period"))
#normality pretty good, potential problems at the edges; possibly acrsin-transformation hasn’t done its job.

#Repeated measures test
anova_test(data=Mites,dv=asin.Mites,wid=Colony,within=Period,between=Treatment)%>%get_anova_table()

p2
#ggsave("Mites23.png",height=7,width=7,unit="in")
#ggsave("Mites23.pdf",height=7,width=7,unit="in",device=pdf)

#Pathogens
Virus<-read_csv("Virus23.csv")%>%
  mutate(Treatment=ifelse(Colony<10|(Colony>20&Colony<31),"Control","Amitraz"),
         Treatment=ordered(Treatment,levels=c("Control","Amitraz")))%>%
  subset(select=c(1,11:13,2:10))

Virus[c(5:13)]<-log(Virus[c(5:13)]+1)

##Did detection rates vary over treatment period/across treatments?
#detection rates
#By period
Per.det.Virus<-Virus%>%pivot_longer(cols=c(5:13),names_to="Pathogen",values_to="D")%>%
  mutate(D=ifelse(D>0,1,0))%>%
  group_by(Year,Period,Treatment,Pathogen)%>%
  summarise(N=n()-sum(is.na(D)),D=sum(D,na.rm=TRUE),ND=N-D)%>%
  subset(select=-N)%>%
  pivot_longer(cols=c(5,6),names_to = "Detect",values_to="N")%>%
  mutate(Period=ordered(Period,levels=c("Pretreatment","Posttreatment")))%>%
  group_by(Treatment,Pathogen)%>%
  nest()%>%
  mutate(tabs=lapply(data,function(df)xtabs(N~Detect+Period,data=df)))%>% 
  mutate(ORS=lapply(tabs,function(x)x+.5),
         ORS=lapply(ORS,function(x)as.data.frame(ORS)),
         ORS=lapply(ORS,function(df)pivot_wider(df,values_from=3,names_from = c(1,2))))%>%
  unnest(ORS)%>%
  mutate(ODDS=(D_Pretreatment/ND_Pretreatment)/(D_Posttreatment/ND_Posttreatment),
         test=lapply(tabs,fisher.test),results=map(test,tidy))%>%
  unnest(results)%>%
  subset(select=c(1:2,9,11,12))%>%
  group_by()%>%
  mutate(p.rank=rank(p.value, ties.method="last"),
         p.crit=(p.rank/n()*0.1),
         sig=ifelse(p.value<=p.crit,"SIG","NS"))

Per.det.Virus

Virus%>%
  mutate(D=ifelse(NUNIV>0,1,0))%>%
  group_by(Treatment,D)%>%
  summarise(N=sum(D))%>%
  pivot_longer(cols=c(2,3),names_to = "Detect",values_to="N")%>%
  xtabs(N~Detect+Treatment,data=.)%>%
  +.5%>%
  as.data.frame()%>%
  pivot_wider(values_from=3,names_from = c(1,2))%>%
  mutate(OR=(D_Control/N_Control)/(D_Amitraz/N_Amitraz))


#By Treatment
Trt.det.Virus<-Virus%>%pivot_longer(cols=c(5:13),names_to="Pathogen",values_to="D")%>%
  mutate(D=ifelse(D>0,1,0))%>%
  group_by(Year,Period,Treatment,Pathogen)%>%
  summarise(N=n()-sum(is.na(D)),D=sum(D,na.rm=TRUE),ND=N-D)%>%
  subset(select=-N)%>%
  pivot_longer(cols=c(5,6),names_to = "Detect",values_to="N")%>%
  mutate(Period=ordered(Period,levels=c("Pretreatment","Posttreatment")))%>%
  group_by(Period,Pathogen)%>%
  nest()%>%
  mutate(tabs=lapply(data,function(df)xtabs(N~Detect+Treatment,data=df)))%>% 
  mutate(ORS=lapply(tabs,function(x)x+.5),
         ORS=lapply(ORS,function(x)as.data.frame(ORS)),
         ORS=lapply(ORS,function(df)pivot_wider(df,values_from=3,names_from = c(1,2))))%>%
  unnest(ORS)%>%
  mutate(ODDS=(D_Control/ND_Control)/(D_Amitraz/ND_Amitraz),
         test=lapply(tabs,fisher.test),results=map(test,tidy))%>%
  unnest(results)%>%
  subset(select=c(1:2,9,11,12))%>%
  group_by()%>%
  mutate(p.rank=rank(p.value, ties.method="last"),
         p.crit=(p.rank/n()*0.1),
         sig=ifelse(p.value<=p.crit,"*","NS"))

Trt.det.Virus

#Did the repeated measure of colony viral load differ by treatment?
#anova_test doesn't have a tidy method, so we built one
#this test includes zeros as meaningfully zero, rather than the standard, "non-detect"
stat.Virus<-
  Virus%>%
  pivot_longer(cols=c(5:13),names_to="Pathogen",values_to = "ln.Copies")%>%
  group_by(Pathogen)%>%
  nest()%>%
  mutate(models=lapply(data,function(df)anova_test(data=df,dv=ln.Copies,wid=Colony,within=Period,between=Treatment)),
         results=lapply(models,function(models)as_tibble(models)),
         results=lapply(results,function(results)pivot_wider(results,values_from = c(2:7),names_from = Effect))
  )%>%
  unnest(results)%>%
  subset(select=c(1,4,7,10:15))%>%
  group_by()%>%
  mutate(r_Treatment=rank(p_Treatment, ties.method="last"),
         c_Treatment=(r_Treatment/n()*0.1),
         s_Treatment=ifelse(p_Treatment<=c_Treatment,"*","NS"),
         r_Period=rank(p_Period, ties.method="last"),
         c_Period=(r_Period/n()*0.1),
         s_Period=ifelse(p_Period<=c_Period,"*","NS"),
         `r_Treatment:Period`=rank(`p_Treatment:Period`, ties.method="last"),
         `c_Treatment:Period`=(`r_Treatment:Period`/n()*0.1),
         `s_Treatment:Period`=ifelse(`p_Treatment:Period`<=`c_Treatment:Period`,"*","NS"))

stat.Virus

#detection table
des.Virus<-Virus%>%pivot_longer(cols=c(5:13),names_to="Pathogen",values_to="ln.Copies")%>%
  group_by(Pathogen,Treatment,Period)%>%
  mutate(det.Copies=ifelse(ln.Copies==0,NA,ln.Copies))%>%
  summarise(Total=n(),
            Detected=sum(ifelse(is.na(det.Copies)==FALSE,1,0)),
            Global.mean=mean(ln.Copies,na.rm=TRUE),
            Global.SE=sd(ln.Copies,na.rm=TRUE)/sqrt(n()),
            Detected.Mean=mean(det.Copies,na.rm=TRUE),
            Detected.SE=sd(det.Copies,na.rm=TRUE)/sqrt(Detected),
  )

Per.det.Virus%>%filter(p.value<0.05)
Trt.det.Virus%>%filter(p.value<0.05)

#Figures
path.labs<-c("ABPV","BQCV","CBPV","DWV-A","DWV-B","IAPV","LSV","Vairimorpha spp.","Trypanosomes","Pre-treatment","Post-treatment")
names(path.labs)<-c("ABPV","BQCV","CBPV","DWVA","DWVB","IAPV","LSVU","NUNIV","TUNIV","Pretreatment","Posttreatment")

#Detection Rate
#Detection figure
P4a<-Virus%>%pivot_longer(cols=c(5:13),names_to="Pathogen",values_to="D")%>%
  mutate(D=ifelse(D>0,1,0),
         Detect=case_when(
           D==0~"Not-detected",
           D==1~"Detected"),
         Detect=ordered(Detect,levels=c("Not-detected","Detected")),
         Period=ordered(Period,levels=c("Pretreatment","Posttreatment"))
         )%>%
  ggplot()+
  geom_bar(aes(x=Treatment,fill=Detect),color="grey50",position = "fill")+
  geom_text(data=des.Virus,aes(x=Treatment,y=Detected/Total,label=paste(Detected,Total,sep="/")),vjust=0)+
  scale_fill_manual(values=c("grey90","indianred"))+
  theme_minimal()+
  facet_grid(cols=vars(Pathogen),rows=vars(Period),
             labeller=as_labeller(path.labs))+
  theme(strip.background=element_rect(fill="slategray2"),
        strip.text=element_text(face="bold",size=11),
        axis.title=element_text(face="bold",size=16),
        axis.text.x=element_text(face="bold",size=12,angle=-30, hjust=0,vjust=0),
        axis.text.y=element_text(size=10),
        legend.position="top",
        legend.text=element_text(face="bold",size=16),
        )+
  labs(y="Proportion Detected",fill=NULL,x=NULL)

#Loads Figure
fig.Virus<-
  full_join(
    Virus%>%
      pivot_longer(cols=c(5:13),names_to="Pathogen",values_to="Copies"),
    rbind(
  stat.Virus%>%
  subset(select=c(Pathogen,s_Treatment))%>%
    mutate(Period="Pretreatment"),
  stat.Virus%>%
    subset(select=c(Pathogen,s_Treatment))%>%
    mutate(Period="Posttreatment")
   )%>%
  mutate(Alpha=ifelse(s_Treatment=="*",1,.35))
  )

des.Virus<-full_join(des.Virus,
          rbind(
            stat.Virus%>%
              subset(select=c(Pathogen,s_Treatment))%>%
              mutate(Period="Pretreatment"),
            stat.Virus%>%
              subset(select=c(Pathogen,s_Treatment))%>%
              mutate(Period="Posttreatment")
          )%>%
            mutate(Alpha=ifelse(s_Treatment=="*",.85,.25)))
P4b<-
ggplot()+
  geom_col(data=des.Virus,aes(x=Treatment,y=Global.mean,fill=Treatment,alpha=Alpha),position=position_dodge())+
  geom_errorbar(data=des.Virus,aes(x=Treatment,ymin=Global.mean-Global.SE,ymax=Global.mean+Global.SE,color=Treatment,alpha=Alpha),width=.25)+
  geom_point(data=fig.Virus,aes(color=Treatment,y=Copies,x=Treatment,alpha=Alpha),position=position_jitterdodge())+
  facet_grid(rows=vars(Period),cols=vars(Pathogen),
             labeller=as_labeller(path.labs))+
  theme_minimal()+
  theme(strip.background=element_rect(fill="slategray2"),
        strip.text=element_text(face="bold",size=11),
        axis.title=element_text(face="bold",size=16),
        axis.text.x=element_text(face="bold",size=12,angle=-30, hjust=0,vjust=0),
        axis.text.y=element_text(size=10),
        legend.position = "none",
  )+
  scale_fill_manual(values=c("forestgreen","firebrick"))+
  scale_color_manual(values=c("forestgreen","firebrick"))+
  labs(y="Copy number (ln)",x=NULL)

ggarrange(P4a,P4b,labels=c("A.","B."),ncol = 1,heights = c(1,.9))
ggsave("Virus.png",width=14,height=10,units="in")
ggsave("Virus.pdf",width=14,height=10,units="in")
