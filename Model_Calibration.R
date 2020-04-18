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
library(tidyr)
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
dt.census <- readRDS('Data/census_cal.RDS')
##############################################################################
# Get Origin Location (Location 5 YEARS AGO)
##############################################################################
# Change previous MSA (MIGMET5) to be non-detailed
dt.census[MIGMET5>0,MIGMET5 := substr(MIGMET5,1,nchar(MIGMET5)-1)]
# In order to get the values of parameters that existed at the time of moving we need to "unmove" people 
dt.census[,state_5 := ifelse(MIGPLAC5 == 990, STATEFIP, 
                             ifelse(MIGRATE5 > 1, MIGPLAC5,
                                    ifelse(MIGRATE5 == 1,STATEFIP,NA)))] # Get trailing state if stayed in same place, new place if moved, or NA if not available
dt.census[,met_5 := ifelse(MIGRATE5 == 1, METAREA, 
                           ifelse(MIGRATE5 > 1, MIGMET5,NA)) ]
##############################################################################
# Correct for missing data in 1980
##############################################################################
# 50% of 1980 data do not have previous metro area but they seem fairly balanced 
# summary(dt.census[MIGPLAC5 == 0  & MIGRATE5 != 1])
# summary(dt.census[!(MIGPLAC5 == 0  & MIGRATE5 != 1) & YEAR == 1980])
dt.census[ ,list(sum(PERWT*(MIGRATE5 %in% c(1)))/sum(PERWT*(MIGRATE5 > 0))), by = list(YEAR)]
# Double PERWT for 1980 to have consistent population based weighting across periods
dt.census[YEAR == 1980, PERWT:= 2*PERWT]
##############################################################################
# Deflate using CPI
#############################################################################
dt.cpi <- fread('Local_Market_Effect/Data/CPI/CPI_change.csv')[,list(YEAR=Year, rel_CPI = `Change from1982-1984`)]
# Get desired Dollar Year
basis <- dt.cpi[YEAR == 2018]$rel_CPI
# Deflate wage for changes in CPI
dt.census <- merge(dt.census,dt.cpi,'YEAR')
dt.census[,INCTOT:= INCTOT/rel_CPI*basis]
dt.census[,RENT:= RENT/rel_CPI*basis]
#############################################################################
# Investigate Differences in Subsetting
##############################################################################
# Get dummy for college
dt.census[,col:= as.numeric(EDUC >= 10)]
# Make subgroups
dt.census_prime <- dt.census[AGE >= 25 & AGE <= 55 & YEAR < 2010 & state_5<=56]
dt.census_prime_ft <- dt.census_prime[WKSWORK1 >= 47 & UHRSWORK > 34]
# Check if fraction of full time employment of low skill is increasing in the fraction of skilled
dt.comp <- dt.census_prime[,list(pop_ft = sum(PERWT*(WKSWORK1 > 47 & UHRSWORK > 34)),
                                 pop_total = sum(PERWT),
                                 frac_col = sum(PERWT*(col))/sum(PERWT),
                                 pop_lt_ft = sum(PERWT*(WKSWORK1 > 47 & UHRSWORK > 34)*(1-col)),
                                 pop_lt = sum(PERWT*(1-col)),
                                 frac_l_ft = sum(PERWT*(WKSWORK1 > 47 & UHRSWORK > 34)*(1-col))/
                                                   sum(PERWT*(1-col))), by = list(YEAR, state_5, met_5)]
ggplot(dt.comp[!is.na(state_5)], aes(frac_col, frac_l_ft,weight=pop_total, color = as.factor(YEAR))) + geom_point(aes(size = pop_total)) + geom_smooth(method='lm')
summary(lm(frac_col ~ frac_l_ft, dt.comp[!is.na(met_5)], weights = dt.comp[!is.na(met_5)]$pop_total))

# dt.comp[,list(lambda = sum(pop_lt_ft)/sum(pop_total)), by = list(YEAR)] # get lambda based on full time workers
# dt.comp[,list(lambda = sum(pop_lt)/sum(pop_total)), by = list(YEAR)]# get lambda based on all workers
dt.comp[,list(lambda = sum(pop_lt)/sum(pop_total))]
summary(lm(pop_lt ~ pop_total-1, dt.comp[!is.na(met_5)], weights = dt.comp[!is.na(met_5)]$pop_total))
summary(lm(pop_lt_ft ~ pop_total-1, dt.comp[!is.na(met_5)], weights = dt.comp[!is.na(met_5)]$pop_total))
##############################################################################
# Calculate Endogenous Parameter Values
##############################################################################
# Get data for the MSA level, pnc, phc, w_l, w_h, L_c, del_c, eta_c
dt.met <- dt.census_prime_ft[state_5 <= 56, list(pnc = sum(PERWT*INCTOT*(1-col))/sum(PERWT*(1-col)),
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
ggplot(dt.met_const, aes(phc,pnc))+geom_point()+geom_smooth(method = 'lm')
dt.met_const[pnc>65000]
# get joint state, msa value
dt.met_const[,MSA := paste0(state_5,'_',met_5)]
l.model <- lm(phc ~ pnc + as.factor(MSA), dt.met_const, weights = dt.met_const$L_l)
vcov_firm <- cluster.vcov(l.model, dt.met_const$MSA)
dt.results <- data.table(tidy(coeftest(l.model,vcov_firm )))
# lambda should be equal to 1 minus the coefficient on the nontradeable good
lambda <- 1-dt.results[term == 'pnc']$estimate
# Alternative calculation
lambda2 <- 1-dt.met_const[,list(sum(L_h)/sum(L))]$V1
lambda3 <- 1- data.table(tidy(lm(L_h~L,dt.met_const,weights = dt.met_const$L)))[term == 'L']$estimate
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
dt.moves <- dt.census_prime_ft[!is.na(state_5) & col == 1 & AGE > 25 & AGE<= 35 ,
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
# Get year to control for changes over time
dt.res_col_m[,YEAR := substr(from_t,nchar(from_t)-3,nchar(from_t))]
# Run the regression
l.model <- lm(log_res ~ log_del_f + as.factor(YEAR)-1, dt.res_col_m, weights = dt.res_col_m$movers)
vcov_firm <- cluster.vcov(l.model, dt.res_col_m$from)
dt.results <- data.table(tidy(coeftest(l.model,vcov_firm )))
# Again separating years
# l.model <- lm(log_res ~ log_del_f, dt.res_col_m[YEAR ==1980], weights = dt.res_col_m[YEAR ==1980]$movers)
# vcov_firm <- cluster.vcov(l.model, dt.res_col_m[YEAR ==1980]$from)
# dt.results <- data.table(tidy(coeftest(l.model,vcov_firm )))
# 
# l.model <- lm(log_res ~ log_del_f, dt.res_col_m[YEAR ==1990], weights = dt.res_col_m[YEAR ==1990]$movers)
# vcov_firm <- cluster.vcov(l.model, dt.res_col_m[YEAR ==1990]$from)
# dt.results <- data.table(tidy(coeftest(l.model,vcov_firm )))
# 
# l.model <- lm(log_res ~ log_del_f, dt.res_col_m[YEAR ==2000], weights = dt.res_col_m[YEAR ==2000]$movers)
# vcov_firm <- cluster.vcov(l.model, dt.res_col_m[YEAR ==2000]$from)
# dt.results <- data.table(tidy(coeftest(l.model,vcov_firm )))
# Back out parameters
gamma_v <- dt.results[term=='log_del_f']$estimate
xi <- dt.results[term!='log_del_f'][,list(YEAR = substr(term,nchar(term)-3,nchar(term)),xi = exp(estimate))]
xtable(dt.results)
ggplot(dt.res_col_m, aes(del_f, mean_p_res, color = YEAR)) + geom_point() + geom_smooth(method = 'lm')
##############################################################################
# Convert impact into dollars
##############################################################################
# Get standard deviation of college fraction
dt.sd <- dt.met[,list(sd = sd(frac_col),
                      mean = mean(frac_col),
                      p20 = quantile(frac_col,.20),
                      p80 = quantile(frac_col,.80)), by = list(YEAR = as.character(YEAR))]
# Merge in xi
dt.sd <- merge(dt.sd,xi,by = 'YEAR')
# Convert a 1sd change into dollars 
dt.sd[,impact := xi*((mean+sd)^gamma_v)-xi*(mean^gamma_v) ]
dt.sd[,impact2 := xi*(p80^gamma_v)-xi*(p20^gamma_v) ]
# Compare to typical variation
dt.met_f <- merge(dt.met_const,xi[,list(xi,YEAR = as.numeric(YEAR))],by = 'YEAR') 
dt.met_f[,numerator := w_h-frac_col^gamma_v-lambda*pnc-phc+amenity]
dt.met_f[,list(sd(numerator)), by = list(YEAR= as.factor(YEAR))]
##############################################################################
# Now Use Wage to pull out elasticity with respect to 
##############################################################################
dt.wage_premium <- dt.met_const[,list(log_w = log(w_h), L, log_frac_col = log(frac_col)), by = list(YEAR, state_5, met_5)]
dt.wage_premium[,MSA := paste0(state_5,'_',met_5)]
l.model <- lm(log_w ~ log_frac_col+as.factor(YEAR)+as.factor(MSA)-1 , dt.wage_premium, weights = dt.wage_premium$L)

vcov_firm <- cluster.vcov(l.model, dt.wage_premium$MSA)
dt.results <- data.table(tidy(coeftest(l.model,vcov_firm )))
gamma_w <- dt.results[term=='log_frac_col']$estimate
theta_t <- dt.results[grepl( 'YEAR',term)][,list(YEAR = substr(term,nchar(term)-3,nchar(term)),theta_t = exp(estimate))]
theta_c <- dt.results[grepl( 'MSA',term)][,list(YEAR = substr(term,15,nchar(term)),theta_c = exp(estimate))]

xtable(dt.results)
##############################################################################
# Convert impact into dollars
##############################################################################
# Get standard deviation of college fraction
dt.sd2 <- merge(dt.sd,theta_t,by = 'YEAR')
# Convert a 1sd change into dollars 
dt.sd2[,impact_w := theta_t*(p80^gamma_w)-theta_t*(p20^gamma_w) ]
# Compare to typical variation
dt.met_f <- merge(dt.met_const,xi[,list(xi,YEAR = as.numeric(YEAR))],by = 'YEAR') 
dt.met_f[,numerator := w_h-frac_col^gamma_v-lambda*pnc-phc+amenity]
dt.met_f[,list(sd(numerator)), by = list(YEAR= as.factor(YEAR))]



