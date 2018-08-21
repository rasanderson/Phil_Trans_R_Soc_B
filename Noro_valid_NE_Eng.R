# Norovirus simulation model
  
set.seed(123)

library(lubridate)
library(deSolve)
library(mc2d)
library(FME)
library(tidyr)
library(ggplot2)
library(lattice)

rm(list=ls())

####### BASIC MODEL STARTING CONDITIONS AND PARAMETERS ##########################
#                    cohort demographics to analyse
#                            NORTHUMBERLAND

humanoids<-read.table("data/humanoids_3yr_end.csv", header=TRUE, sep=",")
humanoids<-data.frame(humanoids)

num_prim_schools <- 4
num_sec_schools  <- 2
num_hospitals    <- 2
num_carehome     <- 2                          

prim_school_dis_stat <- matrix(0,nrow=num_prim_schools,ncol=2)
sec_school_dis_stat  <- matrix(0,nrow=num_sec_schools,ncol=2)
hosp_dis_stat        <- matrix(0,nrow=num_hospitals,ncol=2)
carehome_dis_stat    <- matrix(0,nrow=num_carehome,ncol=2)
cohort     <- nrow(humanoids)
rep_period <- 36
ext_month  <- 0
years      <- 3   

length_run=((365*years)+ext_month)-2   # extra for symptom initiation

max_age  <- 80  # upper age limit
sens_run <- 10

#FIXED MODEL PARAMETERS 

week_counter <- 1
week_temp    <- 0
week_rain    <- 0

###  OTHER PARAMETERS
dayno   <- 1 # monday
weekend <- 0


### TRANSMISSION PROCESSES 
trans_process <- 7       ###       
max_immune_time <- 365

expose_rate<-matrix(0,nrow=trans_process, ncol=1)
lev1             <- 0.010    
lev2             <- 0.002    
oyster_cont_rate <- 0.079   

env_cont_trans_risk <- 0.001
env_part_dose       <- 50

# Setup matrices to store results etc.
NORO_cases     <- matrix(0,nrow=length_run+rep_period+ext_month+10,ncol=8)
NORO_SES_class <- matrix(0,nrow=length_run+rep_period+ext_month+10,ncol=7)
NORO_OVERALL   <- matrix(0,nrow=length_run+rep_period+ext_month+10,ncol=(sens_run*8)+1)
for(q in 1:length_run){
    NORO_OVERALL[q,1]=q
}
institutes <- num_prim_schools+num_sec_schools+num_hospitals

LHS_sample3 <- c(0.001400064, 0.004659069, 0.003631911, 0.001109124, 0.00329738,
                 11.00014, 0.003594912, 0.001907125,   0.002060451,
      1320.305,        253.0543)
print(LHS_sample3)
 


################################# social setting hospital admission #################
num_families <- max(humanoids[,14])

# population over 65 in care homes
prop_carehome   <- 2535/71144
hosp_admin_risk <- 0.0008584973      # per person per day


########################### FUNCTIONS ###############################################
alloc_carehome<-function(num_carehome){
            arbit=1/num_carehome
            home=runif(1)
            home=1+trunc(home/arbit)

            return(home)
            }

alloc_hospital<-function(num_hospitals, hosp_admin_risk){
            # am I in hospital ?
            hospital=0
            admit=runif(1)
            if(admit<=hosp_admin_risk){
                   arbit=1/num_hospitals
                   hospital=runif(1)
                   hospital=1+trunc(hospital/arbit)
                   }
            return(hospital)
            }

immunity_challenge<-function(max_immune_time,immune_stat)
                      {
                      new_immune_stat=0
                      
                      if(immune_stat>0)
                           {
                           grad=1-exp(max_immune_time*log(0.99))

                           new_immune_stat<-grad*immune_stat
                           }
  
                     if(new_immune_stat<0.0001)
                       {
                       new_immune_stat=0
                       }
                      return(new_immune_stat)
                      }

hand_contam<-function()
            {
            prob<-rbeta(1, 1.78, 41.1, ncp = 0) 
            prob1<-rbeta(1, 0.6, 2.3, ncp = 0)  
            prob3<-rpert(1, 0.0, 0.4,1)   
            return(prob*prob1*prob3)
            }

weib_surv<-function(N0,t)
            {
            prop=0
            if(N0>0)
              {
               t=t*24
               b=58.19
               n=0.53
               lNt<-log(N0)- ((1/2.303)*((t/b)^n)) 
               Nt=exp(lNt)
               prop=Nt/N0
               
               }      
            return(prop)
            }

age_increment<-function(age, immune)
               {
               age=age+1
                 
               if(age>=81)
                 {
                 age=1
                 immune=0
                 }               
               return(c(age, immune))
               }


age_classer<-function(age)
             {
             age_class=0
             if(age<=5.99)
               {
               age_class=1
               }

            if(age>6 & age <=10.0)
               {
               age_class=2
               }
           if(age>10.1 & age <=15.0)
               {
               age_class=3
               }
           if(age>15.1 & age <=20.0)
               {
               age_class=4
               }
           if(age>20.1 & age <=40.0)
               {
               age_class=5
               }
           if(age>40.1 & age <=60.0)
               {
               age_class=6
               }
           if(age>60)
               {
               age_class=7
               }
             return(age_class)
             }


illness<-function(dose)
              {
              onoverill=0
              if(dose>0)
                {
                ill=0
                Asym =99.2046
                lrc=-0.6959 
                onoverill<-Asym*(1-exp(-exp(lrc)*(log(dose))))
                onoverill=onoverill/100
                }
              
              return(onoverill)
              }


oysters<-function(time, oyster_cont_rate, lev1,lev2)      
         {
         dosed=0                      
         if(time>365)
           {
           deno=trunc(time/365)
           time=time-(deno*365)
           }
         if(time >273 | time < 90)  
           {
           eat=eat_oysters(lev1,lev2) 
           if(eat>0)
             {
              is_it_cont<-runif(1)  
              if(is_it_cont<=oyster_cont_rate)   
                 {
                 dosed=1                 
                 }
               }
            }
         return(dosed)                

         }


Oyster_contamination<-function(day, G_type)
            {
            cal_day<-day %% 365
           
            cal_month<-1+cal_day %/% 30.4
            cosmo<-cos(2*pi*cal_month/12)
            sinmo<-sin(2*pi*cal_month/12)
            if(G_type==1)
              {
               predo<-cosmo*-0.515+sinmo*-1.246+0.379
               predo<-exp(predo)/(1+exp(predo))   
              }
           if(G_type==2)
              {
               predo<-cosmo*-0.342+sinmo*-1.028+1.02
               predo<-exp(predo)/(1+exp(predo))
              }

           return(predo)
            }


Oyster_part_weight<-function(day, G_type)
            {
            cal_day<-day %% 365
           
            cal_month<-1+cal_day %/% 30.4
            cosmo<-cos(2*pi*cal_month/12)
            sinmo<-sin(2*pi*cal_month/12)
            if(G_type==1)
              {
               predo<-exp(cosmo*-0.42+sinmo*-1.10+3.4)
              
              }
           if(G_type==2)
              {
               predo<-exp(cosmo*-0.55+sinmo*-1.93+3.24)
              
              }
           predo=predo*25
           return(predo)
            }

eat_oysters<-function(lev1,lev2)
             {
            tot1=0
            eat1<-runif(1)
            if(eat1<=lev1)
               {
               eat2<-runif(1)
               if(eat2<=0.1428)
                 {
                 tot1=1
                 }
               }
            lev3=lev1+lev2
            if(eat1>lev1 & eat1<=lev3)
              {
              yesno<-runif(1)  
              if(yesno<=0.033)
                {
                 tot1=1
                 }
              }    
              return(tot1)
           }



part_dose<-function(day, exposed, env_part_dose)
           {
           particles=0
           # 1 pschool; 2 is sschool; 3 hospital 4 family; 5 chome; 6 food; 7 env
           if(exposed >0 & exposed <=5)
             {
              fraction_of_room=0.01
              droplet=0.001  # fraction of ml inhaled
              particles=fraction_of_room*droplet*170000000
              }
           if(exposed==6)
              {
              G_type=2   
              particles=Oyster_part_weight(day, G_type)           
              }
          if(exposed==7)
              {
              particles=runif(1)*env_part_dose
              }
           return(particles)
           }

contagion<-function(x, data, expose, proc)
           {
           df <- data[[x]]
           if(x==1) return(df)
              inf_status<-sum(df[,6])
              susc_num<-(nrow(df))-inf_status
              if(inf_status>=1 & susc_num >=1)
                {
                df[,15]=1   
                for(k in 1:nrow(df))
                   {
                   if(df[k,6]==0)
                     {
                     potential_exposed<-(expose/susc_num)*inf_status   
                     piddle<-runif(1)
                     if(piddle<=potential_exposed)                             
                       {
                       df[k,13]=proc
                       } 
                    }
                  }
                }
           return(df)
          }


# Function setup complete

sens_run=10                                                                                       
for(ko in 1:sens_run)
  {
   humanoids<-read.table("data/humanoids_3yr_end.csv", header=TRUE, sep=",")
   humanoids<-data.frame(humanoids)
   NORO_cases<-matrix(0,nrow=length_run+rep_period+ext_month+10,ncol=8)
   NORO_ages<-matrix(0,nrow=length_run+rep_period+ext_month+10,ncol=7)
   NORO_PRIM<-matrix(0,nrow=length_run+rep_period+ext_month+10,ncol=num_prim_schools+1)
   NORO_SEC<-matrix(0,nrow=length_run+rep_period+ext_month+10,ncol=num_sec_schools+1)
   NORO_HOSPITAL<-matrix(0,nrow=length_run+rep_period+ext_month+10,ncol=num_hospitals+1)
   NORO_CAREHOME<-matrix(0,nrow=length_run+rep_period+ext_month+10,ncol=num_carehome+1)
   expose_rate<-matrix(0,nrow=trans_process, ncol=1)
   LHS_sample3 <- c(0.001400064, 0.004659069, 0.003631911, 0.001109124, 0.00329738,
                    11.00014, 0.003594912, 0.001907125,   0.002060451,
                    1320.305,        253.0543)
   
  for(lh in 1:5)
     {
     expose_rate[lh,1]=LHS_sample3[lh]
     }
    illness_duration=round(LHS_sample3[6])
    lev1=LHS_sample3[7]
    lev1=LHS_sample3[8]
    env_cont_trans_risk=LHS_sample3[9]
    env_part_dose=round(LHS_sample3[10])
    max_immune_time=round(LHS_sample3[11])


  year_count=0
  year_trip=1
  weather_year_count=1
  for(time in 1:length_run)
    {
    cat("run is ",ko,"  time is ",time," \n")
    weekend=0
    if(dayno==6 | dayno == 7)
      {
      weekend=1
      }
     
    NORO_cases[time,8]=time

     rand_proc<-sample(trans_process, replace=FALSE)
      for(counter in 1: trans_process)
         {  
         proc=rand_proc[counter]
         if(proc==1) 
              {              
              if(weekend==0)
                {
                humanoids.primaryschool<-split(humanoids, humanoids[,8])
                primschool.out<-lapply(1:length(humanoids.primaryschool), contagion, data=humanoids.primaryschool, expose=expose_rate[1,1],proc=proc)
                humanoids<-Reduce("rbind",primschool.out)
                }
              }
          if(proc==2) 
              {              
              if(weekend==0)
                {
                humanoids.secondaryschool<-split(humanoids, humanoids[,9])
                secschool.out<-lapply(1:length(humanoids.secondaryschool), contagion, data=humanoids.secondaryschool, expose=expose_rate[2,1],proc=proc)
                humanoids<-Reduce("rbind",secschool.out )
                }
              }
          if(proc==3) 
              {
               humanoids.hospital<-split(humanoids, humanoids[,10])
               hospital.out<-lapply(1:length(humanoids.hospital), contagion, data=humanoids.hospital, expose=expose_rate[3,1],proc=proc)
               humanoids<-Reduce("rbind",hospital.out)              
              }
          if(proc==4) 
              {                     
               humanoids.family<-split(humanoids, humanoids[,14])
               family.out<-lapply(1:length(humanoids.family), contagion, data=humanoids.family, expose=expose_rate[4,1],proc=proc)
               humanoids<-Reduce("rbind", family.out)
              }
          if(proc==5) 
              {                     
               humanoids.carehome<-split(humanoids, humanoids[,11])
               carehome.out<-lapply(1:length(humanoids.carehome), contagion, data=humanoids.carehome, expose=expose_rate[5,1],proc=proc)
               humanoids<-Reduce("rbind", carehome.out)
              }
          } # end of proc

   pick_this<-sample(nrow(humanoids), replace=FALSE)
   for(ind in 1:nrow(humanoids))
      {
      me<-pick_this[ind]
     tripper=0 
     if(humanoids[me,10]>0)
        {
        humanoids[me,12]=humanoids[me,12]-1
        if( humanoids[me,12]<=0)
          {
          humanoids[me,12]=0
          humanoids[me,10]=0
          }
        }
      if(humanoids[me,10]==0)
         {
         humanoids[me,10]<-alloc_hospital(num_hospitals, hosp_admin_risk)
        
         if(humanoids[me,10]>0)
            {
            humanoids[me,12]=7     
            }
         }


      oyster_cont_rate<-Oyster_contamination(time,2)
      eat_oyster<-oysters(time, oyster_cont_rate,lev1,lev2)
      if(eat_oyster>0) 
        {
        humanoids[me,13]<-6                               
        }

      if(humanoids[me,13]==0  & humanoids[me,15]>0 & humanoids[me,6]==0)
                 {             
                 do_I_touch_stuff<-runif(1)
                 act_cont_risk=env_cont_trans_risk
                 if(do_I_touch_stuff<act_cont_risk)
                    {                   
                    humanoids[me,13]=7  
                    humanoids[me,15]=env_part_dose                                          
                    }
                 }             

      if(humanoids[me,13]>0 & humanoids[me,6]==0)
                 {   
                 dose=0
                 dose<-part_dose(time, humanoids[me,13], env_part_dose)          
                 get_ill<-illness(dose)
                 if(humanoids[me,5]<= get_ill)
                     { 
                     proc_rec<-humanoids[me,13]                    
                     NORO_cases[time,proc_rec]=NORO_cases[time,proc_rec]+1
                     humanoids[me,5]=1    # set to ill
                     humanoids[me,6]=1    # I am now ill
                     }
                 }
       if(humanoids[me,13]==0)
                 {                
                 immunity_stat=humanoids[me,5]
                 if(humanoids[me,5]>0)
                    {
                    tempe<-immunity_challenge(max_immune_time,immunity_stat)
                    humanoids[me,5]=tempe
                    }
                 
                 }
        

      if(humanoids[me,15]>0)
        {
        cont_leve<-humanoids[me,15]
        humanoids[me,15]=humanoids[me,15]*weib_surv(cont_leve,1)
        if(cont_leve<=18)
           {
           humanoids[me,15]=0  # if <18 there is zero contamination
           }
        }


       if(year_trip==365)
         {      
            ague<- humanoids[me,2]  
            imm<-humanoids[me,5]             
            increment<-age_increment(ague,imm)
            humanoids[me,2]<-increment[1]
        
            if(increment[2]==1)
              {
              humanoids[me,5]=increment[2]
               humanoids[me,8]=0
              } 
            if(humanoids[me,2]==5)
              {
              humanoids[me,8]<-round(runif(1)*num_prim_schools)    
              humanoids[me,9]<-0
              }
            if(humanoids[me,2]==11)
              {
              humanoids[me,9]<-round(runif(1)*num_sec_schools)    
              humanoids[me,8]<-0
              }
            if(humanoids[me,2]==18)
              {
              humanoids[me,8]<-0    
              humanoids[me,9]<-0
              }
            if(humanoids[me,2]>65)
              {
              humanoids[me,8]<-0    
              humanoids[me,9]<-0
              care=runif(1)
              if(care<=prop_carehome)
                 {
                 home_num=alloc_carehome(num_carehome)
                 humanoids[me,11]=home_num
                 }
 
              }
            if(humanoids[me,2]>81)
              {
              humanoids[me,]=0
              not_found=1
              while(not_found==1)
                   {
                   fam_pick<-round(runif(1)*max(humanoids[,14]))
                   if(humanoids[fam_pick,2]<=40 & humanoids[fam_pick,2]>16 & humanoids[fam_pick,3]==2)
                     {
                     humanoids[me,14]=fam_pick
                     gen_pick<-runif(1)
                     gens=2
                     if(gen_pick<=0.5)
                       {
                       gens=1
                       }
                     humanoids[me,3]=gens
                     humanoids[me,2]=0
                     not_found=0
                     }
                   } 
              }
          }
        # check on disease status
        if(humanoids[me,6]>0 & humanoids[me,7]>=1)
              {
              humanoids[me,7]=humanoids[me,7]+1
              if(humanoids[me,7]>illness_duration) 
                 {
                 humanoids[me,7]=0
                 humanoids[me,6]=0
                 humanoids[me,13]=0
                 }
              }
        if(humanoids[me,6]>0 & humanoids[me,7]==0)
              {
              humanoids[me,7]=1
              }
        # end of year count and incrementation for this individual  
       }
      # end of all individuals

    prim_time=num_prim_schools+1
    sec_time=num_sec_schools+1
    hosp_time=num_hospitals+1
    care_time=num_carehome+1
    for(kok in 1:num_prim_schools)
        {
        outbreak<-subset(humanoids, humanoids[,8]==kok)
        NORO_PRIM[time,kok]=sum(outbreak[,6])
        NORO_PRIM[time,prim_time]=time
        }  
   for(kok in 1:num_sec_schools)
        {
        outbreak<-subset(humanoids, humanoids[,9]==kok)
        NORO_SEC[time,kok]=sum(outbreak[,6])
        NORO_SEC[time,sec_time]=time
        }  

   for(kok in 1:num_hospitals)
        {
        outbreak<-subset(humanoids, humanoids[,10]==kok)
        NORO_HOSPITAL[time,kok]=sum(outbreak[,6])
        NORO_HOSPITAL[time,hosp_time]=time
        }  

   for(kok in 1:num_carehome)
        {
        outbreak<-subset(humanoids, humanoids[,11]==kok)
        NORO_CAREHOME[time,kok]=sum(outbreak[,6])
        NORO_CAREHOME[time,care_time]=time
        }  

     year_trip=year_trip+1
     if(year_trip==366)
       {
        year_count=year_count+1
        year_trip=1
       }
    dayno=dayno+1
    if(dayno==8)
       {
       dayno=1
       }
      cat("NORO_cases up to today ", sum(humanoids[,6]), "\n")

      # for age classes on each day
    for(pp in 1:nrow(humanoids))
       {
       if(humanoids[pp,6]==1 & humanoids[pp,7]==1)
          { 
          age_years<-humanoids[pp,2]
          age_class=age_classer(age_years)
          NORO_ages[time,age_class]=NORO_ages[time,age_class]+1
          }
       }
   for(hum in 1:nrow(humanoids))
      {
      if(humanoids[hum,6]>0 & humanoids[hum,7]==1)
        {
        transtype=humanoids[hum,13]
        sens_TRANS=((ko-1)*7)+1+transtype
        NORO_OVERALL[time,sens_TRANS]=NORO_OVERALL[time,sens_TRANS]+1
        }
      }


   #end of time
   }
} #end of rep loop


#####################  SAVING WORK ################################
print("Final saving session...")
save.image("outputs/valid_rep_NE_Eng.RData")
write.csv(NORO_OVERALL, "outputs/NORO_FINAL_REPS_NE_Eng.csv")
write.csv(NORO_ages,"outputs/NORO_age_1yr_FINAL_REPS_NE_Eng.csv")
write.csv(NORO_PRIM,"outputs/NORO_PRIM_FINAL_REPS_NE_Eng.csv")
write.csv(LHS_sample3,"outputs/NORO_LHS_samples_FINAL_REPS_NE_Eng.csv")
write.csv(NORO_SEC,"outputs/NORO_SEC_FINAL_REPS_NE_Eng.csv")
write.csv(NORO_HOSPITAL,"outputs/NORO_HOSPITAL_FINAL_REPS_NE_Eng.csv")
