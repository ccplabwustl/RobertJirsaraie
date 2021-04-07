
SCALE <- function(x){(x-min(x))/(max(x)-min(x))}
library(data.table)
library(tidyverse)
library(ggplot2) 
library(ggpubr)
library(dplyr)
library(plyr)
library(scales)
set.seed(777)

######
### Create the Simulated Dataset To Test Allocation Methods
######

method1 <- matrix(nrow = 100, ncol = 8)
colnames(method1) <- c("EQU-Unemploy","EQU-Location","EFF-Unemploy","EFF-Social",
                       "EQU-Unemploy-Scaled","EQU-Location-Scaled",
                       "EFF-Unemploy-Scaled","EFF-Social-Scaled")
method2 <- matrix(nrow = 100, ncol = 8)
colnames(method2) <- c("EQU-Unemploy","EQU-Location","EFF-Unemploy","EFF-Social",
                       "EQU-Unemploy-Scaled","EQU-Location-Scaled",
                       "EFF-Unemploy-Scaled","EFF-Social-Scaled")
method3 <- matrix(nrow = 100, ncol = 8)
colnames(method3) <- c("EQU-Unemploy","EQU-Location","EFF-Unemploy","EFF-Social",
                       "EQU-Unemploy-Scaled","EQU-Location-Scaled",
                       "EFF-Unemploy-Scaled","EFF-Social-Scaled")
method4 <- matrix(nrow = 100, ncol = 8)
colnames(method4) <- c("EQU-Unemploy","EQU-Location","EFF-Unemploy","EFF-Social",
                       "EQU-Unemploy-Scaled","EQU-Location-Scaled",
                       "EFF-Unemploy-Scaled","EFF-Social-Scaled")
method5 <- matrix(nrow = 100, ncol = 8)
colnames(method5) <- c("EQU-Unemploy","EQU-Location","EFF-Unemploy","EFF-Social",
                       "EQU-Unemploy-Scaled","EQU-Location-Scaled",
                       "EFF-Unemploy-Scaled","EFF-Social-Scaled")
method6 <- matrix(nrow = 100, ncol = 8)
colnames(method6) <- c("EQU-Unemploy","EQU-Location","EFF-Unemploy","EFF-Social",
                       "EQU-Unemploy-Scaled","EQU-Location-Scaled",
                       "EFF-Unemploy-Scaled","EFF-Social-Scaled")

###Location

for (i in 1:100) {
  df <- as.data.frame(c(rep("Production",20), rep("FoodService",20), rep("Retail",20), rep("Transportation",20), rep("Construction",20)))
  names(df)[1] <- "Business" ; df$Business <- as.factor(df$Business)
LOCATION<-NULL 
for (PROB in c(1, 12, 16, 10, 9)){
	LOCATION<-c(LOCATION,rbinom(20, 2, PROB/20))
}

###Hours

HOURS<-NULL 
for (PROB in c(1, 12, 16, 10, 9)){
	HOURS<-c(HOURS,rbinom(20, 1, PROB/20))
	}

###Social

SOCIAL<-NULL 
for (MEANxSD in c("1500x1000", "750x400", "1000x1000", "3500x1000", "900x500")){
	MEANxSD<-unlist(strsplit(MEANxSD,"x"))
	SOCIAL<-c(SOCIAL,rnorm(20, as.numeric(MEANxSD[1]), as.numeric(MEANxSD[2])))
}

###Employed

EMPLOYED<-NULL 
for (MEANxSD in c("400x250", "50x50", "100x150", "200x500", "400x200")){
	MEANxSD<-unlist(strsplit(MEANxSD,"x"))
	EMPLOYED<-c(EMPLOYED,rnorm(20, as.numeric(MEANxSD[1]), as.numeric(MEANxSD[2])))
}

df[,c("Location","Hours","Social","Employed")] <- c(LOCATION, HOURS, SOCIAL, EMPLOYED)
df$Location<-recode(df$Location, "0" = "Rural", "1" = "Suburb", "2" = "Urban", .default = NA_character_)
df$Hours<-recode(df$Hours, "0" = "Day", "1" = "Night", .default = NA_character_)
df[,c("Social","Employed")] <- abs(round(df[,c("Social","Employed")]))
m1<-lm(Employed~Social,data=df) ; df$Residuals<-m1$residuals

######
### Run Every Allocation Method
######

###Sim-1: Employees
df[,"Incl_Sim1"]<-0
df[order(df$Employed,decreasing=TRUE)[1:50],"Incl_Sim1"]<-1

###Sim-2: Employees + Social
df[,"Incl_Sim2"]<-0
df[order(df$Residuals)[1:50],"Incl_Sim2"]<-1

###Sim-3: Employees + Social with Hall Passes
df[,"Incl_Sim3"]<-0
rdf<-df[-c(order(df$Employed,decreasing=TRUE)[0:25]),]
df[which(df$Residuals %in% rdf[order(rdf$Residuals,decreasing=TRUE)[1:50],"Residuals"]),"Incl_Sim3"]<-1

###Sim-4: Will 50/50
df[,"Incl_Sim4"]<-0
df$sEmployed <- SCALE(df$Employed)
df$SOCIAL=1/df$Social ; df$sSocial<-SCALE(df$SOCIAL)
df$Sim4 <- (df$sEmployed*0.5) + (df$sSocial*0.5)
df<-df[order(df$Sim4, decreasing=TRUE),] ; df[1:50,"Incl_Sim4"]<-1

###Sim-5: Will 75/25
df[,"Incl_Sim5"]<-0
df$sEmployed <- SCALE(df$Employed) 
df$SOCIAL=1/df$Social ; df$sSocial<-SCALE(df$SOCIAL)
df$Sim5 <- (df$sEmployed*0.75) + (df$sSocial*0.25)
df<-df[order(df$Sim5, decreasing=TRUE),] ; df[1:50,"Incl_Sim5"]<-1
df<-subset(df,select = -c(sEmployed,SOCIAL,sSocial,Sim4,Sim5)) ; df<-df[order(as.numeric(row.names(df))),]

###Sim-6: Will 25/75
df[,"Incl_Sim6"]<-0
df$sEmployed <- SCALE(df$Employed) 
df$SOCIAL=1/df$Social ; df$sSocial<-SCALE(df$SOCIAL)
df$Sim6 <- (df$sEmployed*0.25) + (df$sSocial*0.75)
df<-df[order(df$Sim6, decreasing=TRUE),] ; df[1:50,"Incl_Sim6"]<-1
df<-subset(df,select = -c(sEmployed,SOCIAL,sSocial,Sim6)) ; df<-df[order(as.numeric(row.names(df))),]

######
### Evaluate the Dispearsion of Burdens for Each Industry
######

###Number of Open Business By Industry
COUNT<-ddply(df, "Business", summarise, grp.mean=mean(Incl_Sim1)*20) ; names(COUNT)[2]<-"Incl_Sim1"
COUNT<-merge(COUNT,ddply(df, "Business", summarise, grp.mean=mean(Incl_Sim2)*20),by="Business") ; names(COUNT)[3]<-"Incl_Sim2"
COUNT<-merge(COUNT,ddply(df, "Business", summarise, grp.mean=mean(Incl_Sim3)*20),by="Business") ; names(COUNT)[4]<-"Incl_Sim3"
COUNT<-merge(COUNT,ddply(df, "Business", summarise, grp.mean=mean(Incl_Sim4)*20),by="Business") ; names(COUNT)[5]<-"Incl_Sim4"
COUNT<-merge(COUNT,ddply(df, "Business", summarise, grp.mean=mean(Incl_Sim5)*20),by="Business") ; names(COUNT)[6]<-"Incl_Sim5"
COUNT<-merge(COUNT,ddply(df, "Business", summarise, grp.mean=mean(Incl_Sim6)*20),by="Business") ; names(COUNT)[7]<-"Incl_Sim6"

###Percentage of Unemployed By Industry
PERCENT<-COUNT ; PERCENT[,2:ncol(PERCENT)]<-NA 
PERCENT$Business<-as.character(PERCENT$Business) ;PERCENT[nrow(PERCENT)+1,"Business"]<-"Total"
for (SIM in names(df)[grep("Incl_Sim",names(df))]){
	INDEX=0
	for (IND in unique(df$Business)){
		DF<-df[which(df$Business == IND),]
		PCT=(sum(DF[which(DF[,SIM] == 0),"Employed"])/sum(DF$Employed))*100
		INDEX=INDEX+1 ; PERCENT[INDEX,SIM]<-PCT
	}
	PCT=(sum(df[which(df[,SIM] == 0),"Employed"])/sum(df$Employed))*100
	INDEX=INDEX+1 ; PERCENT[INDEX,SIM]<-PCT
}

### Equal Distributions of Industries
EMPLOY<-COUNT ; EMPLOY[,2:ncol(EMPLOY)]<-NA ; UNEMPLOY<-EMPLOY
for (SIM in names(df)[grep("Incl_Sim",names(df))]){
	INDEX=0
	for (IND in unique(df$Business)){
		DF<-df[which(df$Business == IND),]
		SUM=sum(DF[which(DF[,SIM] == 0),"Employed"])
		INDEX=INDEX+1 ; UNEMPLOY[INDEX,SIM]<-SUM
	}
	INDEX=0
	for (IND in unique(df$Business)){
		DF<-df[which(df$Business == IND),]
		SUM=sum(DF[which(DF[,SIM] == 1),"Employed"])
		INDEX=INDEX+1 ; EMPLOY[INDEX,SIM]<-SUM
	}
}
colnames(EMPLOY) <- paste("EMPLOY", colnames(EMPLOY), sep = "_") ; names(EMPLOY)[1]<-"Business"
colnames(UNEMPLOY) <- paste("UNEMPLOY", colnames(UNEMPLOY), sep = "_") ; names(UNEMPLOY)[1]<-"Business"
FINAL1<-NULL ; FINAL1<-merge(EMPLOY,UNEMPLOY,by=c("Business"))

### Equal Distributions of Location
LOCAT<-as.data.frame(matrix(NA, nrow=3, ncol=7)) ; names(LOCAT)[1]<-"Location"
LOCAT[1,1]<-unique(df$Location)[1] ; LOCAT[2,1]<-unique(df$Location)[2] ; LOCAT[3,1]<-unique(df$Location)[3] 
names(LOCAT)[2:ncol(COUNT)]<-names(COUNT)[2:ncol(COUNT)] ; UNLOCAT<-LOCAT
for (SIM in names(df)[grep("Incl_Sim",names(df))]){
	INDEX=0
	for (JOINT in unique(df$Location)){
		DF<-df[which(df$Location == JOINT),]
		SUM=sum(DF[which(DF[,SIM] == 0),"Employed"])
		INDEX=INDEX+1 ; UNLOCAT[INDEX,SIM]<-SUM
	}
	INDEX=0
	for (JOINT in unique(df$Location)){
		DF<-df[which(df$Location == JOINT),]
		SUM=sum(DF[which(DF[,SIM] == 1),"Employed"])
		INDEX=INDEX+1 ; LOCAT[INDEX,SIM]<-SUM
	}
}
colnames(LOCAT) <- paste("EMPLOY", colnames(LOCAT), sep = "_") ; names(LOCAT)[1]<-"Location"
colnames(UNLOCAT) <- paste("UNEMPLOY", colnames(UNLOCAT), sep = "_") ; names(UNLOCAT)[1]<-"Location"
FINAL2<-NULL ; FINAL2<-merge(LOCAT,UNLOCAT,by=c("Location"))

### Aggreate Equality For Business and Location
EQUAL<-COUNT ; EQUAL<-EQUAL[c(1:2),] ; EQUAL[,1]<-c("Unemployment","Location")
names(EQUAL)[1]<-"Metric" ; EQUAL[,2:ncol(EQUAL)]<-NA 
for (IND in names(df)[grep("Incl_Sim",names(df))]){
	COL1=which(names(FINAL1)==paste0("EMPLOY_",IND)) 
	COL2=which(names(FINAL1)==paste0("UNEMPLOY_",IND)) 
	MOD1 <- chisq.test(cbind(FINAL1[,COL1],FINAL1[,COL2]))
	STAT1  <- MOD1$statistic
	MOD2 <- chisq.test(cbind(FINAL2[,COL1],FINAL2[,COL2]))
	STAT2  <- MOD2$statistic
	EQUAL[1,IND]<-STAT1 ; EQUAL[2,IND]<-STAT2
}

######
### Run Projected Social Interactions Over Time
######

TRAJ<-as.data.frame(matrix(0,nrow=1,ncol=3))
for (SIM in names(df)[grep("Incl_Sim",names(df))]){
	OPEN<-df[which(df[SIM]==1),] ; PROB<-1
	for (GEN in 1:52){
		WEEK<-c() 
		for (ROW in 1:50){
			if (OPEN[ROW,"Business"]=="Production"){
				SOCIAL<-rnorm(1, 1500, 1000)
			}
			if (OPEN[ROW,"Business"]=="FoodService"){
				SOCIAL<-rnorm(1, 750, 400)
			}
			if (OPEN[ROW,"Business"]=="Retail"){
				SOCIAL<-rnorm(1, 1000, 1000)
			}
			if (OPEN[ROW,"Business"]=="Transportation"){
				SOCIAL<-rnorm(1, 3500, 1000)
			}
			if (OPEN[ROW,"Business"]=="Construction"){
				SOCIAL<-rnorm(1, 900, 500)
			}
			WEEK<-c(WEEK,SOCIAL)
		}
		if (GEN == 1){
			ADD<-data.frame(V1=c(SIM),V2=GEN,V3=(sum(WEEK)*PROB))
		} else {
			ADD<-data.frame(V1=c(SIM),V2=GEN,V3=(sum(WEEK)*PROB)+TRAJ[nrow(TRAJ),3])
		}
		TRAJ<-rbind(TRAJ,ADD)
	}
}
TRAJ<-TRAJ[-1,]

######
### Aggregate Both Efficency and Equality Outcome Metrics
######

### Raw Values For Each Outcome Metric
names(PERCENT)[1]<-"Metric" ; OUTCOMES<-rbind(EQUAL,PERCENT[nrow(PERCENT),]) ; OUTCOMES[4,]<-NA
for (SIM in 1:6){
	OUTCOMES[4,SIM+1]<-TRAJ[52*SIM, 3]
}
OUTCOMES<-as.data.frame(t(OUTCOMES)[-1,]) ; OUTCOMES[,1:4]<-lapply(OUTCOMES[,1:4],as.numeric)
names(OUTCOMES)<-c("EQU-Unemploy","EQU-Location","EFF-Unemploy","EFF-Social") 

### Scaled Outcomes For Standardized Comparison
RANKED<-OUTCOMES
for (COL in 1:4){
	RANKED[,COL]<- SCALE(OUTCOMES[,COL])
}
RANKED$Average<-rowMeans(RANKED)
for (ROW in 1:6){
	RANKED[ROW,"SD"]<-sd(RANKED[ROW,1:4])
}

method1[i,] = c(as.matrix(OUTCOMES)[1,], as.matrix(RANKED)[1,1:4])
method2[i,] = c(as.matrix(OUTCOMES)[2,], as.matrix(RANKED)[2,1:4])
method3[i,] = c(as.matrix(OUTCOMES)[3,], as.matrix(RANKED)[3,1:4])
method4[i,] = c(as.matrix(OUTCOMES)[4,], as.matrix(RANKED)[4,1:4])
method5[i,] = c(as.matrix(OUTCOMES)[5,], as.matrix(RANKED)[5,1:4])
method6[i,] = c(as.matrix(OUTCOMES)[6,], as.matrix(RANKED)[6,1:4])

print(i)
}


method1 <- as.data.frame(method1)
method2 <- as.data.frame(method2)
method3 <- as.data.frame(method3)
method4 <- as.data.frame(method4)
method5 <- as.data.frame(method5)
method6 <- as.data.frame(method6)
method1$Method <- "1"
method2$Method <- "2"
method3$Method <- "3"
method4$Method <- "4"
method5$Method <- "5"
method6$Method <- "6"

### Group Methods
Allocation <- rbind(method1, method2, method3, method4, method5, method6)
colnames(Allocation) <- c("EQU_Unemploy","EQU_Location","EFF_Unemploy", "EFF_Social" ,        
                          "EQU_Unemploy_Scaled", "EQU_Location_Scaled",
                          "EFF_Unemploy_Scaled", "EFF_Social_Scaled","Method")
Allocation <- mutate(Allocation,
                     avgEQUEFF = (EQU_Unemploy_Scaled+EQU_Location_Scaled+EFF_Unemploy_Scaled+EFF_Social_Scaled)/4)

AllocationEQU <- mutate(Allocation,
                     avgEQUEFF = ((EQU_Unemploy_Scaled*0.4)+(EQU_Location_Scaled*0.4)+(EFF_Unemploy_Scaled*0.1)+(EFF_Social_Scaled*0.1)))

AllocationEFF <- mutate(Allocation,
                     avgEQUEFF = ((EQU_Unemploy_Scaled*0.1)+(EQU_Location_Scaled*0.1)+(EFF_Unemploy_Scaled*0.4)+(EFF_Social_Scaled*0.4)))
######
### Create Figures of Evaluation Metrics
######

### Unemployment Efficiency Box plot for type of business
ggplot(data = Allocation, aes(x = Method, y = EQU_Unemploy/1000)) +
  geom_boxplot()+
  labs(title = "Boxplot of Unemployment Fairness for Business Type",
       x = "Allocation Method", 
       y = "Unemployment Fairness (Thousands)")
### Unemployment Efficiency Box plot for location
ggplot(data = Allocation, aes(x = Method, y = EQU_Location/1000)) +
  geom_boxplot() +
  labs(title = "Boxplot of Unemployment Fairness for Business Location",
       x = "Allocation Method", 
       y = "Unemployment Fairness (Thousands)")
### Unemployment Box plot
ggplot(data = Allocation, aes(x = Method, y = EFF_Unemploy)) +
  geom_boxplot()+
  labs(title = "Boxplot of Unemployment by Allocation Method",
       x = "Allocation Method", 
       y = "Unemployment")
### Social Interactions Box plot
ggplot(data = Allocation, aes(x = Method, y = EFF_Social/1000000)) +
  geom_boxplot()+
  labs(title = "Boxplot of Social Interactions by Allocation Method",
       x = "Allocation Method", 
       y = "Number of Social Interactions After One Year (Millions)")
### Social Interactions Box plot
ggplot(data = Allocation, aes(x = Method, y = EFF_Social/1000000)) +
  geom_boxplot()+
  labs(title = "Boxplot of Social Interactions by Allocation Method",
       x = "Allocation Method", 
       y = "Number of Social Interactions After One Year (Millions)")
### Social Interactions Box plot
ggplot(data = Allocation, aes(x = Method, y = EFF_Social/1000000)) +
  geom_boxplot()+
  labs(title = "Boxplot of Social Interactions by Allocation Method",
       x = "Allocation Method", 
       y = "Number of Social Interactions After One Year (Millions)")
### Average EFF EQU box plot
ggplot(data = AllocationEQU, aes(x = Method, y = avgEQUEFF)) +
  geom_boxplot()+
  labs(title = "Boxplot of Average Efficiency and Equity Score by Allocation Method",
       x = "Allocation Method", 
       y = "Average Efficiency and Equity Score")
### Average EFF EQU box plot
ggplot(data = AllocationEQU, aes(x = Method, y = avgEQUEFF)) +
  geom_boxplot()+
  labs(title = "Boxplot of Weighted Equity over Efficiency Score by Allocation Method",
       x = "Allocation Method", 
       y = "Weighted Equity (80%) over Efficiency (20%) Score")
### Average EFF EQU box plot
ggplot(data = AllocationEFF, aes(x = Method, y = avgEQUEFF)) +
  geom_boxplot()+
  labs(title = "Boxplot of Weighted Efficiency over Equity Score by Allocation Method",
       x = "Allocation Method", 
       y = "Weighted Efficiency (80%) over Equity (20%) Score")

### PRELIM-Hist of Employees by Business
mu <- ddply(df, "Business", summarise, grp.mean=mean(Employed))
ggplot(df, aes(x=Employed, color=Business, fill=Business)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.40) +
	geom_vline(data=mu, aes(xintercept=grp.mean, color=Business),linetype="longdash",size=2) +
	scale_color_manual(values=c("#d11141","#00b159", "#00aedb","#f37735","#ffc425")) +
	scale_fill_manual(values=c("#d11141","#00b159", "#00aedb","#f37735","#ffc425")) +
	theme_classic()

### PRELIM-Hist of Social Interaction by Business
mu <- ddply(df, "Business", summarise, grp.mean=mean(Social))
ggplot(df, aes(x=Social, color=Business, fill=Business)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.40) +
	geom_vline(data=mu, aes(xintercept=grp.mean, color=Business),linetype="longdash",size=2) +
	scale_color_manual(values=c("#d11141","#00b159", "#00aedb","#f37735","#ffc425")) +
	scale_fill_manual(values=c("#d11141","#00b159", "#00aedb","#f37735","#ffc425")) +
	theme_classic()

### PRELIM-Scatterplot of Emlpoyees by Social Interactions
ggplot() + 
	geom_point(size=2,alpha=0.8,aes(x=Employed,y=Social,color=Business),df) +
	geom_smooth(method = "lm",se=F,size=2,alpha=0.9, color="#000000", aes(x=Employed,y=Social),df) +
	scale_color_manual(values=c("#d11141","#00b159","#00aedb","#f37735","#ffc425","#000000")) +
	theme_classic()

### Unemployment Figure
PERCENT$Business<-factor(PERCENT$Business,levels=c("Production", "FoodService", "Retail","Transportation","Construction","Total"))
LONG <- PERCENT %>% pivot_longer(cols = 2:ncol(PERCENT), names_to = "Allocaton", values_to = "Unemployment")
LONG$Allocaton<-gsub("Incl_Sim","",LONG$Allocaton)
ggplot() +  
	geom_point(LONG,mapping=aes(x=Allocaton,y=Unemployment,color=Business),size=3,alpha=0.9,stat = "identity") +
	#geom_line(LONG,mapping=aes(x=Allocaton,y=Unemployment,color=Business,group=Business),size=0.5,alpha=0.35) +
	scale_color_manual(values=c("#d11141","#00b159","#00aedb","#f37735","#ffc425","#000000")) +
	theme_classic()

### Rate of Change Based on Open Businesses
TRAJ$V1<-gsub("Incl_Sim","",TRAJ$V1) ; TRAJ$V3<-TRAJ$V3/1000000
ggplot() +  
	geom_point(TRAJ,mapping=aes(x=V2,y=V3,color=V1),size=1,alpha=0.9,stat = "identity") +
	#geom_line(TRAJ,mapping=aes(x=V2,y=V3,color=V1,group=V1),size=1,alpha=0.5) +
	scale_color_manual("Allocaton",values=c("#d11141","#00b159","#00aedb","#f37735","#ffc425","#632a06","#000000","#858585")) +
	xlab("Time (Weeks)") + ylab("Number of Social Interactions (Millions)") + 
	theme_classic()

### Outcome Figures for Efficency and Equality
OUTCOMES$Allocation<-row.names(OUTCOMES)
LONG <- OUTCOMES %>% pivot_longer(cols = 1:ncol(OUTCOMES)-1, names_to = "Outcome", values_to = "Value")
LONG$Allocation<-gsub("Incl_Sim","",LONG$Allocation)
ggplot() +  
	geom_point(LONG,mapping=aes(x=Allocation,y=Value,color=Allocation),size=3,alpha=0.9,stat = "identity") +
	#geom_line(LONG,mapping=aes(x=Allocation,y=Value,color=Allocation,group=Outcome),size=0.5,alpha=0.35) +
	scale_color_manual("Allocaton",values=c("#d11141","#00b159","#00aedb","#f37735","#ffc425","#632a06","#000000","#858585")) +
	facet_wrap(~ Outcome, scales = "free") + 
	theme_classic()

### Outcome Figures for Efficency and Equality
RANKED$Allocation<-row.names(RANKED)
LONG <- RANKED %>% pivot_longer(cols = 1:ncol(RANKED)-1, names_to = "Outcome", values_to = "Value")
LONG$Allocation<-gsub("Incl_Sim","",LONG$Allocation)
ggplot() +  
	geom_point(LONG,mapping=aes(x=Allocation,y=Value,color=Outcome),size=3,alpha=0.9,stat = "identity") +
	#geom_line(LONG,mapping=aes(x=Allocation,y=Value,color=Outcome,group=Outcome),size=0.5,alpha=0.2) +
	scale_color_manual("Allocaton",values=c("#000000","#d11141","#00b159","#00aedb","#f37735","#ffc425","#632a06","#858585")) +
	theme_classic()


########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######