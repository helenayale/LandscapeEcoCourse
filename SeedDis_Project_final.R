setwd("D:/IntensiveCourse/project")
dispersal<-read.csv(file="dispersal.csv",head=T)
dispersal$proportion<-(dispersal$end/dispersal$start)
head(dispersal)
str(dispersal)

type<-as.character(dispersal[,2])
type[type=="forest"]<-"trees"
type
dispersal$type<-as.factor(type)
vector<-as.character(dispersal[,3])
vector[vector=="walk"]<-"human"
vector
dispersal$vector<-as.factor(vector)

boxplot(proportion~vector, data=dispersal, ylab="Proportion of seeds reaching 2m",xlab="vector")
boxplot(proportion~type:vector, data=dispersal, ylab="Proportion of seeds reaching 2m",xlab="type.vector")
boxplot(proportion~person,data=dispersal)

dispersal$dis_rate<-cbind(sucess=dispersal$end, fail=dispersal$start-dispersal$end)
str(dispersal)
summary(dispersal)


ANOVA_vec<-glm(dis_rate~vector, family="binomial",data=dispersal)
summary(ANOVA_vec)
plot(ANOVA_vec)

ANOVA_null<-glm(dis_rate~1, family="binomial",data=dispersal)
anova(ANOVA_vec,ANOVA_null,test="Chisq")

ANOVA_type<-glm(dis_rate~vector+type+vector:type, family="binomial",data=dispersal)
summary(ANOVA_type)
plot(ANOVA_type)

anova(ANOVA_type,ANOVA_vec, test="Chisq")

library(lme4)
ANOVA_full<-glmer(dis_rate~vector+type+vector:type+(1|person/plot), family="binomial",data=dispersal)
summary(ANOVA_full)
plot(fitted(ANOVA_full),resid(ANOVA_full))
qqnorm(resid(ANOVA_full))
qqline(resid(ANOVA_full))

# predictions
plogis(0.6013)
plogis(0.6013-1.5915)
plogis(0.6013-0.1734)
plogis(0.6013-1.5915-0.1734+0.1762)

# test of overdispersion
library(blmeco)
dispersion_glmer(ANOVA_full)

library(MuMIn)
r.squaredGLMM(ANOVA_full)

# model selection

options(na.action=na.fail)
dis_dredge<-dredge(ANOVA_full)

# best model
best.dis<-get.models(dis_dredge,1)[[1]]
summary(best.dis)
plot(fitted(best.dis),resid(best.dis))
qqnorm(resid(best.dis))
qqline(resid(best.dis))

# predictions
plogis(0.5144)
plogis(0.5144-1.5041)


# R2m:how much of the variance of prediction rate is explained by the fixed effects
# R2c:how much of the variance of prediction rate is explained by the fixed and random effects
r.squaredGLMM(best.dis)

ANOVA_null1<-glmer(dis_rate~1+(1|person/plot),family="binomial",data=dispersal)
anova(best.dis,ANOVA_null1)
