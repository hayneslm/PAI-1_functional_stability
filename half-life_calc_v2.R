#---
#title: "Calculating Halflives"
#output: html_notebook
#---

#Initializing libraries and set working directory

setwd("~/Documents/R_projects/48h_trimmed") 

# install.packages("tidyverse")
# install.packages("data.table")
# install.packages("nlstools")
# install.packages("minpack.lm")
library(tidyverse)
library(data.table)
library(nlstools)
library(minpack.lm)
# install.packages("devtools")
# install.packages("reticulate")
library(reticulate)
library(devtools)



#Setup time and A0 raw data vectors

#Set up Time vector
time <- c(0, 1, 2, 4, 8, 11, 24, 48)
time<-as.data.frame(time)

#A0 Raw Data
A0RC <- c(1,1.02,0.83,0.87,0.62,0.39,0.28,0.13)





#Import data and format to make it more usable

all_data = NULL
for (i in c(1:12)) {
  temp = read_tsv(paste("trimmed_count/Counts_amp",i,".txt",
                        sep = ""))
  all_data = rbind(all_data,temp)
}
A.dt = as.data.table(all_data)

#A0 series extraction
A0.dt = A.dt[, .(Variants, Input = LH1_lane_1, h0 = LH2_lane_1, h1 = LH3_lane_1,h2 = LH4_lane_1,h4 = LH5_lane_1, h8 = LH6_lane_1, h11 = LH7_lane_1,h24 = LH8_lane_1, h48 = LH9_lane_1)]

#Split Variants column apart into components and Assign type
A0.dt = A0.dt %>% separate(Variants, into = c("AA", "rest"),
                           "(?<=[A-Z])(?=[0-9])") %>% 
  separate(rest, into = c("Pos", "Mut"),
           "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(variant = paste(AA, Pos, Mut, sep = ""))
missense.A0 = A0.dt %>% filter(AA != Mut, Mut != 'X') %>% 
  mutate(type = 'missense')
nonsense.A0 = A0.dt %>% filter(Mut == 'X') %>% mutate(type = 'nonsense')
WT.A0 = A0.dt %>% filter(AA == Mut) %>% mutate(type = "WT")

#want to propagate missense and WT through these calculations
misWT.A0 = rbind(missense.A0, WT.A0)

#eliminate rows with zeros from missense.A0
misWT.A0[misWT.A0==0]<-NA
misWT.A0 = misWT.A0 %>% na.omit()

#Now start doing math

#determine fractions of missense+WT sequenced (aka divide by sums)
A0.frac = misWT.A0[, .(variant,
                       Input = Input/sum(Input), h0 = h0/sum(h0), 
                       h1 = h1/sum(h1), 
                       h2 = h2/sum(h2), 
                       h4 = h4/sum(h4), 
                       h8 = h8/sum(h8), 
                       h11 = h11/sum(h11),
                       h24 = h24/sum(h24), 
                       h48 = h48/sum(h48))]
#I was taking averages of the three here before but I will skip that for now

#convert to data.frame and transpose
A0matrixT = as.data.frame(A0.frac) %>% t()
columns<-ncol(A0matrixT) %>% as.numeric()

p1 = vector()
p2<-vector()
p3 = vector()

for (i in c(1:columns)){ 
  print(i)
  title = A0matrixT[1,i]
  print(title)
  
  x<-time
  x<-data.frame(x)
  y<-A0matrixT[2:9,i]
  y<-as.numeric(y)
  y<-y*A0RC
  y<-y/y[1]
  y<-as.data.frame(y)

  holder<-data.frame(x,y)
  x<-as.matrix(holder[1])
  y<-as.matrix(holder[2])
  
  error = try(nlsLM(y ~ a*(exp(-b*x)+c*(exp(-0.007*x))), 
                start = list(a=0.1, b=0.33 ,c=0.9),
                control = nls.control(maxiter = 10000)))
  
  if(class(error) != "try-error"){
  output = nlsLM(y ~ a*(exp(-b*x)+c*(exp(-0.007*x))), 
                  start = list(a=0.1, b=0.33 ,c=0.9),
                  control = nls.control(maxiter = 10000))
  print(output)
  parameters = coef(output)

  
  }else{
    parameters=c("error","error","error")
  } 
  
  print(parameters)
  p1<-as.vector(parameters)
  p2<-rbind(p2,p1)
}
mutations = A0matrixT %>% t() %>% as.data.table() %>% select(variant) 
p3 = cbind(mutations, p2)

#Now manipulate to determine half-lives
halflife_list = log(2)/as.numeric(p3$V2)
halflives = cbind(mutations,halflife_list)

#Normalize by WT 2h halflife
halflife_called = left_join(halflives, misWT.A0, by = 'variant') %>% 
  select('variant', 'AA', 'Pos', 'Mut', 'type', halflife_list) %>% 
  rename(halflife = halflife_list)

positions = unique(halflife_called$Pos)
HL_adj = NULL

for (pos in positions) {
  temp = halflife_called %>% filter(Pos == pos) 
  WThl = temp %>% filter(type == 'WT') %>% select(halflife) %>% as.numeric()
  fudge = 2/WThl
  temp2 = temp %>% mutate(Adj_HL = (temp$halflife)*fudge) %>% print()
  HL_adj = rbind(HL_adj, temp2)
}
