library(MASS)
library(car)
bel_pan2020 <- read.csv("~/bp_R/final_data.csv")

bel_pan2020 <- transform(bel_pan2020,Diet=factor(Diet))
bel_pan2020 <- transform(bel_pan2020,Habitat=factor(Habitat))
bel_pan2020 <- transform(bel_pan2020,Stratum=factor(Stratum))
bel_pan2020 <- transform(bel_pan2020,Family=factor(Family))

###DETERMINE BEST DISTRIBUTION FOR DA
fit1<-fitdistr(bel_pan3$DA, "normal")
fit2<-fitdistr(bel_pan3$DA+.000005, "lognormal")
fit3<-fitdistr(bel_pan3$DA+.000005, "gamma")
fit4<-fitdistr(bel_pan3$DA+.000005, "weibull")
fit5<-fitdistr(bel_pan3$DA+.000005, "exponential")
AIC(fit1,fit2,fit3, fit4,fit5)
#Gamma is best


###TEST FOR SIGNIFICANCE IN SINGLE PARAMETERS
model2020_stratum<-glm(DA+.000005~Stratum, family="Gamma",data=bel_pan2020)
Anova(model2020_stratum)
summary(model2020_stratum)

model2020_habitat<-glm(DA+.000005~Habitat, family="Gamma",data=bel_pan2020)
Anova(model2020_habitat)
summary(model2020_habitat)

model2020_diet<-glm(DA+0.00005~Diet, family="Gamma",data=bel_pan2020)
Anova(model2020_diet)
summary(model2020_diet)

model2020_mass<-glm(DA+.000005~Mass, family="Gamma",data=bel_pan2020)
Anova(model2020_mass)
summary(model2020_mass)


###TUKEY TEST FOR MULTIPLE COMPARISONS
library(multcomp)

model2020_diet<-glm(DA+0.00005~Diet, family="Gamma",data=bel_pan2020)
posthoc_diet<-glht(model2020_diet, linfct=mcp(Diet="Tukey"))
summary(posthoc_diet)


####Let's include family as a random effect
library(glmmADMB)

model100 <-glmmadmb(DA+.000005~Habitat+(1|Family), family="Gamma",data=bel_pan2020)
model101<-glmmadmb(DA+.000005~Stratum+(1|Family), family="Gamma",data=bel_pan2020)
model102<-glmmadmb(DA+.000005~Diet+(1|Family), family="Gamma",data=bel_pan2020)
model103<-glmmadmb(DA+.000005~Mass+(1|Family), family="Gamma",data=bel_pan2020)
model104<-glmmadmb(DA+.000005~Habitat+Stratum+(1|Family), family="Gamma",data=bel_pan2020)
model105<-glmmadmb(DA+.000005~Habitat+Diet+(1|Family), family="Gamma",data=bel_pan2020)
model106<-glmmadmb(DA+.000005~Habitat+Mass+(1|Family), family="Gamma",data=bel_pan2020)
model107<-glmmadmb(DA+.000005~Stratum+Diet+(1|Family), family="Gamma",data=bel_pan2020)
model108<-glmmadmb(DA+.000005~Stratum+Mass+(1|Family), family="Gamma",data=bel_pan2020)
model109<-glmmadmb(DA+.000005~Diet+Mass+(1|Family), family="Gamma",data=bel_pan2020)
model110<-glmmadmb(DA+.000005~Habitat+Stratum+Diet+(1|Family), family="Gamma",data=bel_pan2020)
model111<-glmmadmb(DA+.000005~Habitat+Stratum+Mass+(1|Family), family="Gamma",data=bel_pan2020)
model112<-glmmadmb(DA+.000005~Habitat+Diet+Mass+(1|Family), family="Gamma",data=bel_pan2020)
model113<-glmmadmb(DA+.000005~Stratum+Diet+Mass+(1|Family), family="Gamma",data=bel_pan2020)
model114<-glmmadmb(DA+.000005~Habitat+Stratum+Diet+Mass+(1|Family), family="Gamma",data=bel_pan2020)

AIC(model100,model101,model102,model103,model104,model105,model106,model107,model108,model110,model111,model112,model113,model114)

library(AICcmodavg)
bel_pan2020 <- transform(bel_pan2020,Diet=factor(Diet))


posthoc102<-glht(model102, linfct=mcp(Diet="Tukey"))
summary(posthoc102)

AIC(model102)
AICc(model102, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = NULL)


###REDO ANOVA ON DIETXStratum
model2020_diet2<-glm(DA+0.00005~Diet*Stratum, family="Gamma",data=bel_pan2020)
Anova(model2020_diet2)
summary(model2020_diet2)
AICc(model2020_diet2, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = NULL)




####Four panel plot
dev.new(width=4, height=15)

pdf(file='plot_amazing7.pdf', width=3.4, height=8.7)
par(mfrow=c(4,1), mai=c(0.25,0.5,0.25,0.05))

temp<-tapply(bel_pan2020$DA, bel_pan2020$Diet, mean)
temp<-temp[c("Frug/Nect","Arthropod","Mixed")]
temp.SE<-tapply(bel_pan2020$DA, bel_pan2020$Diet, sd)/sqrt(tapply(bel_pan2020$DA, bel_pan2020$Diet, length))
barplot(temp, ylab="DA (net nucleotide divergence)", xlab="Diet classification", xaxt="n", ylim=c(0,0.04))
abline(h=0)
axis(1, at=c(0.7, 1.9, 3.1), labels=c("Frug/Nect","Arthropod","Mixed-Diet"))


arrows(x0=c(0.7, 1.9, 3.1), x1=c(0.7, 1.9, 3.1), y0=temp-temp.SE, y1=temp+temp.SE, angle=90, length=.05, code=3)

segments(x0=0.7, x1=1.9, y0=0.018, y1=0.018)

text(x=0.7, y=0.010, label="a", cex=1)
text(x=1.9, y=0.019, label="b", cex=1)
text(x=3.1, y=0.025, label="b", cex=1)


temp<-tapply(bel_pan2020$DA, bel_pan2020$Habitat, mean)
temp.SE<-tapply(bel_pan2020$DA, bel_pan2020$Habitat, sd)/sqrt(tapply(bel_pan2020$DA, bel_pan2020$Habitat, length))
barplot(temp, ylab="DA (net nucleotide divergence)", xlab="Habitat preference", xaxt="n", ylim=c(0,0.045))
arrows(x0=c(0.7, 1.9), x1=c(0.7, 1.9), y0=temp-temp.SE, y1=temp+temp.SE, angle=90, length=.05, code=3)
abline(h=0)
axis(1, at=c(0.7, 1.9), labels=c("Forest","Open/Edge"))

temp<-tapply(bel_pan2020$DA, bel_pan2020$Stratum, mean)
temp.SE<-tapply(bel_pan2020$DA, bel_pan2020$Stratum, sd)/sqrt(tapply(bel_pan2020$DA, bel_pan2020$Stratum, length))
barplot(temp, ylab="DA (net nucleotide divergence)", xlab="Canopy stratum", xaxt="n", ylim=c(0,0.045))
arrows(x0=c(0.7, 1.9), x1=c(0.7, 1.9), y0=temp-temp.SE, y1=temp+temp.SE, angle=90, length=.05, code=3)
abline(h=0)
axis(1, at=c(0.7, 1.9), labels=c("Canopy","Understory"))

plot(DA~Mass, data=bel_pan2020, ylab="DA (net nucleotide divergence)", xlab="Body mass (g)", ylim=c(0,0.045), bty="l")
glm1<-glm(DA+0.00005~Mass, data=bel_pan2020, family="Gamma")
predict.data<-predict(glm1, newdata=data.frame(Mass=1:200), type="response", se.fit=T)
lines(x=1:200, y=predict.data$fit)
lines(x=1:200, y=predict.data$fit+predict.data$se.fit, lty=2)
lines(x=1:200, y=predict.data$fit-predict.data$se.fit, lty=2)

dev.off()