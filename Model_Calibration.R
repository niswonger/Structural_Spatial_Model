# This script is designed to calibrate a structural spacial model with educational choice
# It uses census data 
##############################################################################
# Set up Libraries 
##############################################################################
library(readr)
library(data.table)
library(ggplot2)
library(xtable)
library(foreign)
library(usmap)
library(multiwayvcov)
library(AER)
setwd('~/Documents/Harvard/Research/College_Response')
##############################################################################
# Bring in Census Data
##############################################################################
# This data is taken from 1980,1990,2000 census and 5 year ACS 2010 from IPUMS
# Pull in Census data
# dt.census <- fread("Data/Census/20190802_wide_data_from_1980_2000_census.csv",
#                    select = c("YEAR","SAMPLE","SERIAL", "PERNUM","PERWT","WKSWORK1","EMPSTAT","UHRSWORK",
#                               "STATEFIP","COUNTYFIP","METAREA","MIGRATE5","MIGMET5","MIGPLAC5","AGE","EDUC","EDUCD","OCC2010","IND","INCWAGE",
#                               "INCTOT","RENT", "VALUEH"))
# saveRDS(dt.census,'Data/census_cal.RDS')
dt.census <- readRDS('Data/census_cal.RDS')[AGE => 25 & AGE <= 55 & WKSWORK1 >= 47 & UHRSWORK > 34]
##############################################################################
# Calculate Endogenous Parameter Values
##############################################################################
# Change MIGMET5 to be non-detailed
dt.census[MIGMET5>0,MIGMET5 := substr(MIGMET5,1,nchar(MIGMET5)-1)]
# In order to get the values of parameters that existed at the time of moving we need to "unmove" people 
dt.census[,state_5 := ifelse(MIGPLAC5 == 990, STATEFIP, 
                             ifelse(MIGRATE5 > 1, MIGPLAC5,
                                    ifelse(MIGRATE5 == 1,STATEFIP,NA)))] # Get trailing state if stayed in same place, new place if moved, or NA if not available
dt.census[,met_5 := ifelse(MIGRATE5 == 1, METAREA, 
                           ifelse(MIGRATE5 > 1, MIGMET5,NA)) ]
# 50% of 1980 data do not have previous metro area but they seem fairly balanced 
# summary(dt.census[MIGPLAC5 == 0  & MIGRATE5 != 1])
# summary(dt.census[!(MIGPLAC5 == 0  & MIGRATE5 != 1) & YEAR == 1980])
dt.census[ ,list(sum(PERWT*(MIGRATE5 %in% c(1)))/sum(PERWT*(MIGRATE5 > 0))), by = list(YEAR)]
# Double PERWT for 1980 to have consistent population based weighting across periods
dt.census[YEAR == 1980, PERWT:= 2*PERWT]
# Get dummy for college
dt.census[,col:= as.numeric(EDUC >= 10)]
# Get data for the MSA level, pnc, phc, w_l, w_h, L_c, del_c, eta_c
dt.met <- dt.census[state_5 <= 56, list(pnc = sum(PERWT*INCTOT*(1-col))/sum(PERWT*(1-col)),
                                        phc = 12*sum(PERWT*RENT*(RENT>1))/sum(PERWT*(RENT>1)),
                                        frac_col = sum(PERWT*col)/sum(PERWT),
                                        w_h = sum(PERWT*INCTOT*col)/sum(PERWT*col),
                                        w_l =  sum(PERWT*INCTOT*(1-col))/sum(PERWT*(1-col)),
                                        L = sum(PERWT),
                                        L_l = sum(PERWT*(1-col)),
                                        L_h = sum(PERWT*col)), by = list(YEAR, state_5, met_5)]
# Find consistent MSAs
dt.constant_MSAs <- dt.met[,list(.N), by = list(state_5,met_5)][N==3]
# Focus on MSAs that are consistent throughout the period 
dt.met_const <- merge(dt.met, dt.constant_MSAs, c('state_5','met_5'))
##############################################################################
# Low Skilled Indifference
# In this section we attempt to calculate amenity values and the parameter of lambda which make people indifferent.
# First we assume that amenities are constant in MSAs overtime to jointly estimate amenities and spending on local good. 
##############################################################################
# get joint state, msa value
dt.met_const[,MSA := paste0(state_5,'_',met_5)]
l.model <- lm(phc ~ pnc + as.factor(MSA)-1, dt.met_const, weights = dt.met_const$L_l)
vcov_firm <- cluster.vcov(l.model, dt.met_const$MSA)
dt.results <- data.table(tidy(coeftest(l.model,vcov_firm )))
# lambda should be equal to 1 minus the coefficient on the nontradeable good
lambda <- 1-dt.results[term == 'pnc']$estimate
# Amenities should just be equal to the coefficient on the MSA dummy 
a <- dt.results[grepl('[0-9]',term), list(MSA = substr(term,15,nchar(term)), estimate)][order(estimate)]
# Using a single value of lambda we should be able to get some sense of how amenities may be changing over time. 
dt.met_const[,amenity := phc - (1-lambda)*pnc]
a_t <- dcast(dt.met_const, MSA ~ YEAR, value.var = 'amenity')
a_msa <- dt.met_const[,list(a_sd = sd(amenity), a_m = mean(amenity)), by = list(MSA)]
##############################################################################
# High Skilled Optimization
# Here we invert the equation for migration to get the elasticity of the penalty function
##############################################################################
# Get information from skilled movers: 
dt.moves <- dt.census[!is.na(state_5) & col == 1 & AGE > 25 & AGE<= 35 ,
                      list(movers = sum(PERWT)), by = list(from = paste0(state_5,'_',met_5), 
                                                               to = paste0(STATEFIP,'_',METAREA),
                                                               YEAR)]
dt.moves_f <- merge(dt.moves[,list(MSA=from,to,movers,YEAR)], dt.met_const, by = c('MSA','YEAR'))
setnames(dt.moves_f,c('MSA','to'), c( 'from','MSA'))
dt.moves_t <- merge(dt.moves_f, dt.met_const, by = c('MSA','YEAR'), suffixes = c('_f','_t'))
# Get the number of movers who stayed put
dt.moves_stay <- merge(dt.moves_t[from == MSA, list(movers), by = list(YEAR, from)], dt.moves_t,by = c('from','YEAR'), suffixes = c('_f','_t'), all.y = T)
# Invert the probability equation to get the choices of where people are going
dt.moves_stay[,p_res := log(movers_f)-log(movers_t)-(w_h_f -w_h_t) + lambda*(pnc_f -pnc_t)+(phc_f-phc_t)-(amenity_f-amenity_t)]
# label by year and from 
dt.moves_stay[,from_t := paste0(from,'_',YEAR)]
##########################################
# test
dt.test <- dt.moves_stay[,sum(movers_t*(MSA!=from)/sum(movers_t)), by = list(frac_col_f)]
ggplot(dt.test, aes(frac_col_f,V1))+geom_point() + geom_smooth(method = 'lm')
##################################################
# Get mean of this residual by year and from county
dt.res <- dt.moves_stay[p_res>0,list(mean_p_res = sum(p_res*movers_t)/sum(movers_t)), by = list(from_t,del_f = frac_col_f)]
# Merge in the total number of movers 
dt.movers_f <- dt.moves_stay[from != MSA, list(movers = sum(movers_t)), by = list(from_t = paste0(from,'_',YEAR))]
dt.res_col_m <- merge(dt.movers_f, dt.res, by  = 'from_t')
# Get variables of interest for regression
dt.res_col_m[,log_del_f := log(del_f)]
dt.res_col_m[,log_res := log(mean_p_res)]
dt.res_col_m[,from := substr(from_t, 1, nchar(from_t)-5)] # Get from again so that we can cluster by MSA to account for the correlation in time
# Run the regression
l.model <- lm(log_res ~ log_del_f, dt.res_col_m, weights = dt.res_col_m$movers)
vcov_firm <- cluster.vcov(l.model, dt.res_col_m$from)
dt.results <- data.table(tidy(coeftest(l.model,vcov_firm )))
# Back out parameters
gamma_v <- dt.results$estimate[2]
xi <- exp(dt.results$estimate[1])






ggplot(dt.res_col_m, aes(del_f, mean_p_res)) + geom_point() + geom_smooth(method = 'lm')





