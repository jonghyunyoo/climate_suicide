source("scripts/functions.R")


### Bootstrap CI for main specification ###

#Load US data
data_rok <- read.csv("inputs/SuicideData_ROK.csv")

###########################################################
#The below code has been commented out because it
#takes ~30 minutes to run. The stored data
#from thesese bootstrap runs will be read-in below.
#Uncomment out this section to run bootstraps and replace
#file `BootstrapMainModel_rok.csv' with updated data.      
###########################################################      
# run bootstrap of base model for ROK (county-month + state-year FE). this takes a while. 
# cc <- unique(data_rok$fips)
# set.seed(1)
#     out <- c()
#     for (i in 1:1000) {
#       tryCatch( { samp <- data.frame(fips=sample(cc,length(cc),replace=T))  #use tryCatch() because occasionally weird samples throw an error
#           subdata <- inner_join(data_rok,samp)
#           reg <- summary(felm(rate_adj ~ tmean + prec | fipsmo + stateyear, weights = subdata$popw, data = subdata))$coefficients[c("tmean"),"Estimate"]
#           out <- c(out,reg)
#       }, error=function(e){})
#       print(i)
#       }
#     write_csv(data.frame(x = 1:length(out), est = out),path = "inputs/bootstrap_runs/BootstrapMainModel_rok.csv")


### Run alternative specifications ###

#ROK
data_rok$yearmonth <- data_rok$yr*100+data_rok$month
xx=-20:30

mod1 <- felm(rate_adj ~ tmean + prec | fipsmo + stateyear | 0 | fips , data=data_rok, weights=data_rok$popw)  # stateyear FE
yy = data.frame(xx,stateyear=coef(mod1)[1]*xx)
mod2 <- felm(rate_adj ~ tmean + prec | fipsmo + stateyear | 0 | fips , data=data_rok)  # state-year no weights
yy = data.frame(yy,noweight=coef(mod2)[1]*xx)
mod3 <- felm(rate_adj ~ tmean + prec | fipsmo + year | 0 | fips , data=data_rok, weights=data_rok$popw)  # year FE
yy = data.frame(yy,year=coef(mod3)[1]*xx)
mod4 <- felm(rate_adj ~ tmean + prec + as.factor(state)*yr | fipsmo + year | 0 | fips , data=data_rok, weights=data_rok$popw)  # year FE + time trend
yy = data.frame(yy,yearTT=coef(mod4)[1]*xx)
mod5 <- felm(rate_adj ~ tmean + prec | fipsmo + yearmonth | 0 | fips , data=data_rok, weights=data_rok$popw)  # year-month FE
yy = data.frame(yy,yearmonth=coef(mod5)[1]*xx)
