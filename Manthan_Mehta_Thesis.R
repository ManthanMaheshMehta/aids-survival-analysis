library(survminer)
library(survival)
library(ggplot2)
library(randomForestSRC)
library(ggRandomForests)
library(partykit)
library(gbm)
library(pec)
library(survAUC)
aids_data <- read.csv("AIDS_Clinical_Trials_Study_175.csv")


#Kaplan Meier Curves
fit_str2 <- survfit(Surv(time,cid)~str2,data = aids_data)
ggsurvplot(fit_str2,data = aids_data,pval = TRUE,conf.int = TRUE,
           risk.table = FALSE,title="Kaplan-Meier Curve by Antiretroviral History",
           xlab="Time",ylab="Survival Probability",legend.title="History")

fit_symptom <- survfit(Surv(time, cid) ~ symptom, data = aids_data)
ggsurvplot(fit_symptom, data = aids_data, pval = TRUE, conf.int = TRUE,
           risk.table = F, title = "Kaplan-Meier Curve by Symptom Status",
           xlab = "Time", ylab = "Survival Probability", legend.title = "Symptom")

fit_oprior<-survfit(Surv(time,cid)~oprior,data=aids_data)
ggsurvplot(fit_oprior,data=aids_data,pval=T,conf.int=T,
           risk.table=F,title="Kaplan-Meier Curve by Prior Non-ZDV Therapy",
           xlab="Time",ylab="survival probability",legend.title="Prior Therapy")


#Cox proportional hazard model
cox_model <- coxph(Surv(time,cid)~trt+age+wtkg+karnof+race+gender+
                     cd40+cd420+cd80+cd820+str2+symptom+
                     offtrt+oprior+z30,data=aids_data)

supremum_test<-cox.zph(cox_model)

print(supremum_test)


cor(aids_data$z30,aids_data$str2)

cor(aids_data$str2, aids_data$oprior)

cox_model_z30 <- coxph(Surv(time, cid) ~ age + wtkg + karnof + race + gender +
                         cd80 + cd820 + symptom + oprior + z30 +
                         strata(trt, offtrt),
                       data = aids_data)
summary(cox_model_z30)

cox_model_str2 <- coxph(Surv(time, cid) ~ age + wtkg + karnof + race + gender +
                          cd80 + cd820 + symptom + oprior + str2 +
                          strata(trt, offtrt),
                        data = aids_data)
summary(cox_model_str2)


#Random Survival Forest
rsf_model_full<-rfsrc(Surv(time,cid)~age+karnof+cd420+wtkg+cd40+race+
                        cd80+cd820+symptom+oprior+trt+offtrt+age+gender+z30,
                        data=aids_data,
                        ntree=1000,
                      importance = T,
                      na.action = "na.impute")
print(rsf_model_full)

plot(gg_vimp(rsf_model_full))+
  ggtitle("Variable Importance-Random Survival Forest")

rsf_model <- rfsrc(Surv(time, cid) ~ age + karnof + cd420 + wtkg +
                     cd40 + cd80 + cd820 + symptom + oprior  + trt + offtrt,
                   data = aids_data,
                   ntree = 1000,
                   importance = TRUE,
                   na.action = "na.impute")
print(rsf_model)

plot(gg_vimp(rsf_model)) +
  ggtitle("Variable Importance - Random Survival Forest")


aids_data$SurvObj <- with(aids_data,Surv(time,cid))

ctree_model<-ctree(SurvObj~age + karnof + cd420 + wtkg +
                     cd40 + cd80 + cd820 + symptom + oprior  + trt + offtrt,
                   data = aids_data)

plot(ctree_model,main="Interpretable Survival Tree")


aids_data$SurvObj <- with(aids_data, Surv(time, cid))

gbm_model <- gbm(
  formula = SurvObj ~ age + wtkg + karnof + race + gender + cd80 + cd820 + symptom + oprior + z30 + cd40 + cd420 + trt + offtrt,
  data = aids_data,
  distribution = "coxph",
  n.trees = 3000,
  interaction.depth = 3,
  shrinkage = 0.01,
  n.minobsinnode = 15,
  cv.folds = 5,
  verbose = TRUE
)

best_iter <- gbm.perf(gbm_model, method = "cv")

summary(gbm_model, n.trees = best_iter)

summary(gbm_model, n.trees = best_iter, plotit = TRUE, cBars = length(gbm_model$var.names))



rsf_pred <- predict(rsf_model, newdata = aids_data)$survival
rsf_timepoint <- which.min(abs(rsf_model$time.interest - 1000))  


rsf_surv_probs <- rsf_pred[, rsf_timepoint]

# Calculating  C-index
rsf_cindex <- UnoC(Surv(aids_data$time, aids_data$cid), Surv(aids_data$time, aids_data$cid), 1 - rsf_surv_probs)
print(rsf_cindex)


gbm_pred <- predict(gbm_model, newdata = aids_data, n.trees = best_iter, type = "response")

#  C-index for GBM
gbm_cindex <- UnoC(Surv(aids_data$time, aids_data$cid), Surv(aids_data$time, aids_data$cid), 1 - gbm_pred)
print(gbm_cindex)

