# F1 and F2 Plasticity Comparision
# Code by Baker and Sommer

# Packages 
library(dplyr)
library(lmerTest)
library(lubridate)

# Behavior 

Behavior_2023 <- read.csv("Data/Behavior_Raw.csv")
# Raw X, Y data in given in 4cm x 4cm grid coordinates. Converting to cm
Behavior_2023$Y <- Behavior_2023$Y*4
Behavior_2023$X <- Behavior_2023$Y*4
Behavior_2023$Distance <- Behavior_2023$Distance*4

# Start time should be adjusted because it was often hard to find grasshoppers in the first couple time steps before 9am leading to patchy coverage
Behavior_2023$Time <- hms(Behavior_2023$Time)
Behavior_2023 <- subset(Behavior_2023,Time>(hms("8:40:00")))


# Calculating Mean Canopy Height for that Day
Behavior_2023_Terraria_Mean <- aggregate(Y ~ Population +
                                            Predation.Treatment +
                                            Generation +
                                            Day +
                                            Terraria, 
                                          FUN=mean, data=Behavior_2023)
Behavior_2023_Terraria_Mean

# Because there was only one grasshopper for each terrarium we can also calculate the "cumulative distance" traveled by the grasshopper in the X,Y plane. I'm not sure how robust this metric is, but cool to look at!

Behavior_2023_Terraria_Distance<-aggregate(Distance ~ Population + 
                                              Predation.Treatment + 
                                              Generation + 
                                              Day + 
                                              Terraria, 
                                            FUN=sum, data=Behavior_2023)
Behavior_2023_Terraria_Distance


# Respiration 

## Data Processing

Respiration_Full<-read.csv( "Data/Respriation_2023_F1_F2_Raw.csv")
Respiration_Full<-Respiration_Full[,-1]

Tempcatlist <- data.frame(matrix(ncol = 1, nrow = 2))
Tempcatlist[,1]<-c("25C","30C")
# For internal consistency with measurements 2022 I measured these grasshoppers at a 2min dwell time unfortunate this lead to frequent saturation of the censor making the 35C data patchy. Here I include only the data for 25C and 30C

Final_output_all_temps <- data.frame(matrix(ncol = 9, nrow = 126))

for (Z in 1:2){
  Final_output_Per_Temp <- data.frame(matrix(ncol = 8, nrow = 63))
  for(R in 1:9){
    Raw<-subset(Respiration_Full, Run==R & Tempcat==Tempcatlist[Z,1])
    
    names(Raw)[1] <- "time"
    names(Raw)[65] <- "respiration"
    output<-subset(Raw,time>18& time<33&respiration>0)
    Logger_Pro_Calcd_Resp<-output[,65]
    Logger_Pro_Calcd_Resp
    
    df <- data.frame(matrix(ncol = 6, nrow = 7))
    
    colnames(df) <- c('delay_low','delay_high','end','Stop_manual','Temp_stop_14mins','Logger_Pro_Calcd')
    df$delay_low<-c(18,20,22,24,26,28,30)
    df$delay_low<-df$delay_low+0.2
    df$delay_high<-c(18,20,22,24,26,28,30)
    df$delay_high<-df$delay_high+1.5
    df$end<-c(18,20,22,24,26,28,30)
    df$end<-df$end+2
    
    names(Raw)[1] <- "time"
    names(Raw)[2] <- "CO2"
    names(Raw)[4] <- "temp"
    names(Raw)[15] <- "flow"
    
    
    for(x in 1:7) {
      
      
      Pulse_and_Background<-subset(Raw,time>df[x,1]&time<=df[x,2])
      Background<-subset(Raw,time>df[x,2]&time<=df[x,3])
      
      X<-(sum(Pulse_and_Background$CO2*(1/60))-(sum(Background$CO2*(1/60))/.5)*1.3)
      df[x,4]<-(((X*1000)/1000000)*mean(Pulse_and_Background$flow))/14
      
      Incubation<-subset(Raw,time>(df[x,3]-16)&time<=(df[x,3]))
      df[x,5]<-mean(Incubation$temp)
    }
    
    df[,6]<-(Logger_Pro_Calcd_Resp)/14
    df<-df[,4:6]
    df$chamber<-c(2,3,4,5,6,7,8)
    
    df$Run<-R
    Final_output_Per_Temp[(((R)*7)-6):((R)*7),]<-df
    
    Final_output_Per_Temp$TempCat<-Tempcatlist[Z,1]
  }
  Final_output_all_temps[((Z*63)-62):((Z*63)),]<-Final_output_Per_Temp[1:63,]
}
Final_output_all_temps<-Final_output_all_temps[,4:9]
colnames(Final_output_all_temps) <- c('Chamber','Run','Manual_Resp_Calc','Incubation_Temp','Logger_Pro_Calcd','Temp_Catagorical')

# Logger Pro Automatically Calculates the integral for respiration (Logger_Pro_Calcd) I also calculated the values manually (Manual_Resp_Calc)

Respriation_Mass_2023<-read.csv("Data/Respriation_Mass_2023_F1_F2.csv")
Respriation_Mass_2023<-arrange(Respriation_Mass_2023, Run,Chamber)
Respriation_Mass_2023_times2<-rbind(Respriation_Mass_2023,Respriation_Mass_2023)
Respiration_Baker_2023_Final<-cbind(Final_output_all_temps,Respriation_Mass_2023_times2)
Respiration_Baker_2023_Final

write.csv(Respiration_Baker_2023_Final,"Respiration_F1_F2_2023_Summary_Stats.csv")



## Sample Analysis 

Respiration_Baker_2023_Final<-read.csv("Respiration_F1_F2_2023_Summary_Stats.csv")

# Calculating Mass Specific Resp
Respiration_Baker_2023_Final$Mass_Specific_Resp<-Respiration_Baker_2023_Final$Logger_Pro_Calcd/Respiration_Baker_2023_Final$Mass

# Assigning Regions
A<-subset(Respiration_Baker_2023_Final, Population=="FN"| Population=="SC"|Population=="MC"|Population=="UP")
A$Region<-"West"
B<-subset(Respiration_Baker_2023_Final, Population=="YM"| Population=="HF"|Population=="SP"|Population=="DC")
B$Region<-"East"

Respiration_Baker_2023_Final<-rbind(A,B)

# Removing Saturated and Dead Grasshoppers
Respiration_Baker_2023_Final<-subset(Respiration_Baker_2023_Final, Exclude!=1)
Respiration_Baker_2023_Final

# Prep For Analysis
Respiration_Baker_2023_Final$Individual<-as.character(Respiration_Baker_2023_Final$Individual)
Respiration_Baker_2023_Final$Generation<-as.character(Respiration_Baker_2023_Final$Generation)


# Here I just treated population as a fixed effect instead of nesting it within Region due to issues with convergence. Also due to limited replication there is no special treatment of grasshoppers depending on their cage or orgin
M1<-lmer(Mass_Specific_Resp~Generation*Population*Temp_Catagorical+(1|Individual),data=Respiration_Baker_2023_Final)
anova(M1)
qqnorm(resid(M1))
qqline(resid(M1))
shapiro.test(resid(M1))
leveneTest(resid(M1) ~ Generation*Population*Temp_Catagorical, data = Respiration_Baker_2023_Final ) 
