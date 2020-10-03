

####HOW TO MAKE LAMARC G GRAPH
library(car)
library(multcomp)
library(MASS)

read.csv("~/Desktop/diet_by_g6.csv")
diet_by_g<-read.csv("~/Desktop/diet_by_g6.csv")
diet_by_g<-diet_by_g[1:37,]
diet_by_g
hist(diet_by_g$G118)
hist(log(diet_by_g$G118))




#####THIS IS TO TEST DIFFERENT DISTRIBUTIONS
fit1<-fitdistr(diet_by_g$G118, "normal")
fit2<-fitdistr(diet_by_g$G118, "lognormal")
fit3<-fitdistr(round(diet_by_g$G118), "negative binomial")
fit4<-fitdistr(round(diet_by_g$G118), "geometric")
fit5<-fitdistr(round(diet_by_g$G118), "exponential")
fit6<-fitdistr(round(diet_by_g$G118), "weibull")



AIC(fit1,fit2,fit3,fit4,fit5,fit6)
shapiro.test(diet_by_g$G)

glm0<-lm(log(G+118)~DIET, diet_by_g)
glm01<-lm(log(G+118)~DIET+POP, diet_by_g)
AIC(glm0,glm01)
Anova(glm0)
car::anova(glm01)
summary(glm0)
plot(glm0)
shapiro.test(log(diet_by_g$G+118))
hist(log(diet_by_g$G+118))
 
glm1<-glm(G118~DIET, family=negative.binomial(1), diet_by_g)
glm2<-glm(G118~DIET+POP, family=negative.binomial(1), diet_by_g)
glm3<-glm(G118~DIET*POP, family=negative.binomial(1), diet_by_g)
glm4<-glm(G118~POP, family=negative.binomial(1), diet_by_g)
AIC(glm1,glm2, glm3, glm4)

Anova(glm2)
summary(glm2)

anova(glm1,glm2, test="Chisq")#test for significance of POP
anova(glm4,glm2, test="Chisq")#test for significance of DIET


plot(G~CLASS, diet_by_g)

library(ggplot2)
ggplot(diet_by_g, aes(CLASS, G118, fill=DIET)) + geom_boxplot() +scale_fill_manual(values=c("white","gray")) +theme_bw()  +theme(legend.background = element_rect())+theme(panel.grid.major = element_line())

+ theme(panel.background=element)

means<-tapply(diet_by_g$G,list(diet_by_g$POP,diet_by_g$DIET), mean)
SEs<-tapply(diet_by_g$G,list(diet_by_g$POP,diet_by_g$DIET), sd)/sqrt(tapply(diet_by_g$G,list(diet_by_g$POP,diet_by_g$DIET), length))


interaction.plot(diet_by_g$POP, diet_by_g$DIET, diet_by_g$G118)




fit1<-fitdistr(diet_by_g$G, "normal")
hist(diet_by_g$G)

fit4<-fitdistr(round(diet_by_g$G), "geometric")
hist(diet_by_g$G)
AIC(fit1,fit4)


R2<-read.csv("~/Desktop/Workbook11.csv")
R2
hist(R2$SGRTP)
fit1<-fitdistr(R2$P, "normal")
fit2<-fitdistr(R2$P, "lognormal")
fit3<-fitdistr(round(R2$P), "negative binomial")
fit4<-fitdistr(round(R2$P), "geometric")
fit5<-fitdistr(round(R2$P), "exponential")
AIC(fit1,fit2,fit3,fit4,fit5)

fit100<-fitdistr(R2$LOGP, "normal")

R2GLM<-glm(P~DIET, family=gaussian, R2)
R2GLM1<-glm(P~POPULATION, family=gaussian, R2)
R2GLM2<-glm(P~ DIET + POPULATION, family=gaussian, R2)
AIC(R2GLM, R2GLM1, R2GLM2)

anova(R2GLM1, R2GLM2, test="Chisq")
interaction.plot(R2$POPULATION, R2$DIET, R2$P)

