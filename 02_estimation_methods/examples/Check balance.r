##############################################
# Read in packages and set working directory #
##############################################

# Read in packages #
options(java.parameters = "- Xmx1024m")
library(foreign)
library(xlsx)
library(gdata)
library(gtools)
library(MatchIt)
library(lme4)

# Set working directory #
setwd("C:/Users/mcitkowicz/American Institutes for Research in the Behavioral Sciences/TFA AC Evaluation - !Methods Paper/simulation/code")

####################
# Method 1 (units) #
####################

# Read in data #
dat_full <- read.csv("test/student data.csv")
dat_M1 <- read.csv("test/matched data (M1).csv")

# Check balance #
source("Check balance function.r")

balance.M1 <- fun.balance(dat_full = dat_full,
						   dat_match = dat_M1,
						   tx.var = "D",
						   vars = c("Z_k","r_k","W_jk","u_jk","X_jk","X_ijk","U_ijk"))

balance.M1

write.csv(balance.M1,file="test/matching results (M1).csv",row.names=FALSE)

#######################
# Method 2 (clusters) #
#######################

# Read in data #
dat_full <- read.csv("test/student data.csv")
dat_M2 <- read.csv("test/matched data (M2).csv")

# Merge into student data #
dat_M2 <- merge(dat_full,dat_M2[c("teacher_id","distance","subclass","weights")],by=c("teacher_id"))

# Check balance #
source("Check balance function.r")

balance.M2 <- fun.balance(dat_full = dat_full,
						   dat_match = dat_M2,
						   tx.var = "D",
						   vars = c("Z_k","r_k","W_jk","u_jk","X_jk","X_ijk","U_ijk"))

balance.M2

write.csv(balance.M2,file="test/matching results (M2).csv",row.names=FALSE)

####################
# Method 4 (units) #
####################

# Read in data #
dat_full <- read.csv("test/student data.csv")
dat_M4 <- read.csv("test/matched data (M4).csv")

# Check balance #
source("Check balance function.r")

balance.M4 <- fun.balance(dat_full = dat_full,
						   dat_match = dat_M4,
						   tx.var = "D",
						   vars = c("Z_k","r_k","W_jk","u_jk","X_jk","X_ijk","U_ijk"))

balance.M4

write.csv(balance.M4,file="test/matching results (M4).csv",row.names=FALSE)

#######################
# Method 5 (clusters) #
#######################

# Read in data #
dat_full <- read.csv("test/student data.csv")
dat_M5 <- read.csv("test/matched data (M5).csv")

# Merge into student data #
dat_M5 <- merge(dat_full,dat_M5[c("teacher_id","distance","subclass","weights")],by=c("teacher_id"))

# Check balance #
source("Check balance function.r")

balance.M5 <- fun.balance(dat_full = dat_full,
						   dat_match = dat_M5,
						   tx.var = "D",
						   vars = c("Z_k","r_k","W_jk","u_jk","X_jk","X_ijk","U_ijk"))

balance.M5

write.csv(balance.M5,file="test/matching results (M5).csv",row.names=FALSE)

####################
# Method 7 (units) #
####################

# Read in data #
dat_full <- read.csv("test/student data.csv")
dat_M7 <- read.csv("test/matched data (M7).csv")

# Check balance #
source("Check balance function.r")

balance.M7 <- fun.balance(dat_full = dat_full,
						   dat_match = dat_M7,
						   tx.var = "D",
						   vars = c("Z_k","r_k","W_jk","u_jk","X_jk","X_ijk","U_ijk"))

balance.M7

write.csv(balance.M7,file="test/matching results (M7).csv",row.names=FALSE)

#######################
# Method 8 (clusters) #
#######################

# Read in data #
dat_full <- read.csv("test/student data.csv")
dat_M8 <- read.csv("test/matched data (M8).csv")

# Merge into student data #
dat_M8 <- merge(dat_full,dat_M8[c("teacher_id","distance","subclass","weights")],by=c("teacher_id"))

# Check balance #
source("Check balance function.r")

balance.M8 <- fun.balance(dat_full = dat_full,
						   dat_match = dat_M8,
						   tx.var = "D",
						   vars = c("Z_k","r_k","W_jk","u_jk","X_jk","X_ijk","U_ijk"))

balance.M8

write.csv(balance.M8,file="test/matching results (M8).csv",row.names=FALSE)

#####################
# Method 10 (units) #
#####################

# Read in data #
dat_full <- read.csv("test/student data.csv")
dat_M10 <- read.csv("test/matched data (M10).csv")

dat_M10 <- m_10

# Create weights for comparison students selected more than once #
for(i in unique(dat_M10$student_id)){
  dat_M10[which(dat_M10$student_id==i),"N"] <- nrow(dat_M10[which(dat_M10$student_id==i),])
}
dat_M10$weights <- 1/dat_M10$N

# Drop the duplicate matched cases #
dat_M10 <- dat_M10[!duplicated(dat_M10[c("student_id")]),]

# Merge into student data #
dat_M10 <- merge(dat_full,dat_M10[c("student_id","p3","PS3B","PS3W","pscore","pscore.j","PSW3B","PSW3W","Site_match","subclass","N","weights")],by=c("student_id"))

# Check balance #
source("Check balance function.r")

balance.M10 <- fun.balance(dat_full = dat_full,
						   dat_match = dat_M10,
						   tx.var = "D",
						   vars = c("Z_k","r_k","W_jk","u_jk","X_jk","X_ijk","U_ijk"))

balance.M10

write.csv(balance.M10,file="test/matching results (M10).csv",row.names=FALSE)

############################
# Prepare data - Method 11 #
############################

# Read in data #
dat_full <- read.csv("test/student data.csv")
dat_M11 <- read.csv("test/matched data (M11).csv")

# Create weights for comparison teachers selected more than once #
for(i in unique(dat_M11$teacher_id)){
  dat_M11[which(dat_M11$teacher_id==i),"N"] <- nrow(dat_M11[which(dat_M11$teacher_id==i),])
}
dat_M11$weights <- 1/dat_M11$N

# Drop the duplicate matched cases #
dat_M11 <- dat_M11[!duplicated(dat_M11[c("teacher_id")]),]

# Merge into student data #
dat_M11 <- merge(dat_full,dat_M11[c("teacher_id","p3","PS3B","PS3W","pscore","pscore.j","PSW3B","PSW3W","Site_match","subclass","N","weights")],by=c("teacher_id"))

# Check balance #
source("Check balance function.r")

balance.M11 <- fun.balance(dat_full = dat_full,
						   dat_match = dat_M11,
						   tx.var = "D",
						   vars = c("Z_k","r_k","W_jk","u_jk","X_jk","X_ijk","U_ijk"))

balance.M11

write.csv(balance.M11,file="test/matching results (M11).csv",row.names=FALSE)
