library(dagitty)
library(dplyr)
library(survival)
################################################
# Define data-generazing model (DAG) == g
# including exposure and network-variables (a-p)
g <- dagitty('dag{
            Exposure -> b[beta=.15] Exposure -> l[beta=.15] Exposure -> o[beta=.15];
             a -> b[beta=.3] b -> c[beta=.2] c -> d[beta=.5] d -> e[beta=.3] e -> f[beta=.6]
             b -> k [beta=.3] k -> l [beta=.3] l -> m [beta=.5] ;
             m -> o [beta=.3] o -> p [beta=.5];
             a -> v [beta=.2] ;
             c -> v [beta=.5] ;
             d -> v [beta=.3] ;
             e -> m [beta=.5] ;
             k -> o [beta=.3] ;
             o -> p [beta=.65] }')
################################################
# Simulate data based on the DAG g
x <- simulateSEM(g, N=20000)

################################################
# Define simulation function for survival time
# based on the Gompertz distribution

survsim.cw <- function(object,IV1=X1,IV2=X2,IV3=X3,beta1,beta2,beta3)
{
    n<-20000

    # Variation of these parameters will modify the distribution of survival times
    # (and thus the prevalence of the event at a given censoring date)
    U = round(runif(1, 0, 1),2)  	#Parameter of Gompertz distribution -> between 0 and 1
    alpha = 0.2138				      	#Parameter of Gompertz distribution
    lambdaEvent =  0.7*10 **(-11)	#Parameter of Gompertz distribution

    X1<-object[,IV1]
    X2<-object[,IV2]
    X3<-object[,IV3]
    RV1<-rnorm(n,mean=10,sd=1)
    RV2<-rnorm(n,mean=10,sd=1)
    #print(X1[1:10])
    #print(X2[1:10])
    #print(X3[1:10])

    expBetaX = exp (beta1*X1+beta2*X2+beta3*X3+ RV1 + RV2)

    fup_time = (1 / alpha) * log(1 - ((alpha * log(U)) / (lambdaEvent * expBetaX)))


    #find the first x% event times; here one third
    sort_time = sort(fup_time )

    # simulate different censoring scenarios by varying
    # five<-n/20 and p5 = (sort_time[five]+sort_time[five+1])/2
    # or   third<-n/3 and   p33 = (sort_time[third]+sort_time[third+1])/2
    ten<-n/10
    p10 = (sort_time[ten]+sort_time[ten+1])/2


    # calculate censoring variable
    cens = matrix(0,n,)
    for(i in 1:n){
        if (fup_time[i] < p10){cens[i]<-1}
        else {fup_time[i] <-p10}
    }
    cens<-as.numeric(cens)
    #print(mean(fup_time))
    print(max(fup_time))
    #print(min(fup_time))
    print(sum(cens))
    simSurv<-bind_cols(data.frame(fup_time,cens))
}
#Simulate survival times that depend on specific variables in the DAG-based data in x
SurvTime<-c()
DAG_SURV_EXP<-data.frame()
SURV<-c()
## generate a survival object where survtime depends on influential variables (IV) with specific betas
SurvTime<-data.frame(survsim.cw(object=x,IV1="c" ,IV2= "k",IV3= "p",beta1=0.5 ,beta2=0.5 ,beta3= -0.5))
## merge with exposure & network data
x_SurvTime<-bind_cols(x,SurvTime)
### check: are specified network variables associated with survival-time
### fit<-coxph(Surv(fup_time,cens)~c+k+p,data =x_SurvTime )
### summary(fit)


#Option 1: add non-informative variables
#x_SurvTime1<-x_SurvTime%>%mutate(NoEffect1=rnorm(20000,mean=0,sd=1),ID=c(1:20000))

##Option 2: creating a case-cohort
## draw random subsample (subcohort)
#x_SurvTime_SC<-x_SurvTime%>%sample_n(5000, replace=FALSE)
## ...and oversample (all) cases
#x_SurvTime_Case<-x_SurvTime%>%dplyr::filter(fup_time<max(x_SurvTime$fup_time))
## prepare for Prentice weighted analysis by setting defining a difference between start and stop
## close to 0 for external cases and setting an indicator for subcohort membership sc = 0
#x_SurvTime_extCase<-setdiff(x_SurvTime_Case,x_SurvTime_SC)%>%mutate(sc=0,start=(fup_time-0.0005),stop_t=fup_time)
## ... and the start and end of follow-up as start and stop for all members of the subcohort and  sc = 1
#x_SurvTime_SC<-x_SurvTime_SC%>%mutate(sc=1,start=0,stop_t=fup_time)
## merge subcohort and external cases
#x_SurvTime_CC<-bind_rows(x_SurvTime_SC,x_SurvTime_extCase)
#rm(x_SurvTime_Case)
#rm(x_SurvTime_extCase)

