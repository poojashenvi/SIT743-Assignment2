#========================================================================================================================
#1.5) Generate Bayesian network and perform the d-separation tests
#========================================================================================================================
library(igraph)
library(ggm)
dag <- DAG(R ~ S,T ~ S,N ~ U+S,P ~ F+R+N, C ~ N+T)
drawGraph(dag, adjust = FALSE)
dSep(dag, first="S", second="U", cond=NULL)
dSep(dag, first="F", second="U", cond=c("N", "P"))

#========================================================================================================================
#2) Belief Network
#========================================================================================================================
library("gRain")
source("https://bioconductor.org/biocLite.R")
biocLite("RBGL")
library(RBGL)
library(gRbase)
library(gRain)
biocLite("Rgraphviz")
#define the appropriate network and use the “compileCPT()”function to Compile list of conditional probability tables and create the network.
lh <- c("low","high")
lhn <- c("low","high","normal")
eh <- cptable(~economyhealth, values=c(60,40),levels=lh)
bp.oil <- cptable(~bp|oil, values=c(20,60,20,70,20,10),levels=lhn)
oil.eh <- cptable(~oil|economyhealth, values=c(30,70,15,85),levels=lh)
rt.inf.eh <- cptable(~retailerstock|inflationrate:economyhealth, values=c(40,60,70,30,35,65,80,20,90,10,60,40),levels=lh)
inf.oil.eh <- cptable(~inflationrate|oil:economyhealth, values=c(15,80,5,40,30,30,65,20,15,19,1,80),levels=lhn)
plist <- compileCPT(list(eh,bp.oil,oil.eh,rt.inf.eh,inf.oil.eh))
plist
#plist$oil
#plist$eh 
net1 <- grain(plist)
net1
summary(net1)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")
plot(net1)

summary(plist)

# Probability of Oil prices when BP stock price is low and retailer-stock is high
net12 <- setEvidence(net1,nodes=c("bp", "retailerstock"), states=c("low", "high"))
#The probability of observing this evidence under the model is
#pEvidence( net12 )
querygrain( net12, nodes=c("oil") )

# Probability of BP Stock price when inflation rate is high
net13 <- setEvidence(net1,nodes=c("inflationrate"), states=c("high"))
#The probability of observing this evidence under the model is
#pEvidence( net13 )
querygrain( net13, nodes=c("bp") )

#========================================================================================================================
#4.1) Belief Network # Bayesian Structure Learning 
#========================================================================================================================
library (bnlearn)
# load the ‘insurance’data.
data(insurance)
summary(insurance)

#Build the true network structure
library(bnlearn)
#create and plot the network structure.
modelstring = paste0("[Age][Mileage][SocioEcon|Age][GoodStudent|Age:SocioEcon]",
 "[RiskAversion|Age:SocioEcon][OtherCar|SocioEcon][VehicleYear|SocioEcon:RiskAversion]",
 "[MakeModel|SocioEcon:RiskAversion][SeniorTrain|Age:RiskAversion]",
 "[HomeBase|SocioEcon:RiskAversion][AntiTheft|SocioEcon:RiskAversion]",
 "[RuggedAuto|VehicleYear:MakeModel][Antilock|VehicleYear:MakeModel]",
 "[DrivingSkill|Age:SeniorTrain][CarValue|VehicleYear:MakeModel:Mileage]",
 "[Airbag|VehicleYear:MakeModel][DrivQuality|RiskAversion:DrivingSkill]",
 "[Theft|CarValue:HomeBase:AntiTheft][Cushioning|RuggedAuto:Airbag]",
 "[DrivHist|RiskAversion:DrivingSkill][Accident|DrivQuality:Mileage:Antilock]",
 "[ThisCarDam|RuggedAuto:Accident][OtherCarCost|RuggedAuto:Accident]",
 "[MedCost|Age:Accident:Cushioning][ILiCost|Accident]",
 "[ThisCarCost|ThisCarDam:Theft:CarValue][PropCost|ThisCarCost:OtherCarCost]")
dag = model2network(modelstring)
graphviz.plot(dag)

#a) 100 (first 100 data)
#bic
bic_net1 = hc(insurance[0:100,], score = "bic")
bic_net1
graphviz.plot(bic_net1)
scoreBic100<-score(bic_net1, insurance[0:100,], type = "bic")
print(scoreBic100)

#bde
bde_net1 = hc(insurance[0:100,], score = "bde")
bde_net1
graphviz.plot(bde_net1)
scoreBde100<-score(bde_net1, insurance[0:100,], type = "bde")
print(scoreBde100)

#b) 500 (first 500 data)
#bic
bic_net2 = hc(insurance[0:500,], score = "bic")
bic_net2
graphviz.plot(bic_net2)
scoreBic500<-score(bic_net2, insurance[0:500,], type = "bic")
print(scoreBic500)

#bde
bde_net2 = hc(insurance[0:500,], score = "bde")
bde_net2
graphviz.plot(bde_net2)
scoreBde500<-score(bde_net2, insurance[0:500,], type = "bde")
print(scoreBde500)

#c) 1000 (first 1000 data)
#bic
bic_net3 = hc(insurance[0:1000,], score = "bic")
bic_net3
graphviz.plot(bic_net3)
scoreBic1000<-score(bic_net3, insurance[0:1000,], type = "bic")
print(scoreBic1000)

#bde
bde_net3 = hc(insurance[0:1000,], score = "bde")
bde_net3
graphviz.plot(bde_net3)
scoreBde1000<-score(bde_net3, insurance[0:1000,], type = "bde")
print(scoreBde1000)

#d) 5000 (first 5000 data)
#bic
bic_net4 = hc(insurance[0:5000,], score = "bic")
bic_net4
graphviz.plot(bic_net4)
scoreBic5000<-score(bic_net4, insurance[0:5000,], type = "bic")
print(scoreBic5000)

#bde
bde_net4 = hc(insurance[0:5000,], score = "bde")
bde_net4
graphviz.plot(bde_net4)
scoreBde5000<-score(bde_net4, insurance[0:5000,], type = "bde")
print(scoreBde5000)

#e) 15000 (first 15000 data)
#bic
bic_net5 = hc(insurance[0:15000,], score = "bic")
bic_net5
graphviz.plot(bic_net5)
scoreBic15000 <-score(bic_net5, insurance[0:15000,], type = "bic")
print(scoreBic15000)

#bde
bde_net5 = hc(insurance[0:15000,], score = "bde")
bde_net5
graphviz.plot(bde_net5)
scoreBde15000 <-score(bde_net5, insurance[0:15000,], type = "bde")
print(scoreBde15000)

#e) All data
#bic
bic_net6 = hc(insurance, score = "bic")
bic_net6
graphviz.plot(bic_net6)
scoreBicAll <-score(bic_net6, insurance, type = "bic")
print(scoreBicAll)

#bde
bde_net6 = hc(insurance, score = "bde")
bde_net6
graphviz.plot(bde_net6)
scoreBdeAll <-score(bde_net6, insurance, type = "bde")
print(scoreBdeAll)

#========================================================================================================================
#3) R code for question 3 (extra)
#========================================================================================================================
library("gRain")
source("https://bioconductor.org/biocLite.R")
#biocLite("RBGL")
library(RBGL)
library(gRbase)
library(gRain)
#biocLite("Rgraphviz")
#define the appropriate network and use the “compileCPT()”function to Compile list of conditional probability tables and create the network.
lh <- c("0","1")
A <- cptable(~A, values=c(20,80),levels=lh)
B <- cptable(~B, values=c(15,85),levels=lh)
C.A.B <- cptable(~C|A:B, values=c(10,90,20,80,30,70,80,20),levels=lh)
D.C <- cptable(~D|C, values=c(5,95,30,70),levels=lh)
plist <- compileCPT(list(A,B,C.A.B,D.C))
plist
net1 <- grain(plist)
net1
summary(net1)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")
plot(net1)
summary(plist)


net12 <- setEvidence(net1,nodes=c("B"), states=c("0"))
#The probability of observing this evidence under the model is
#pEvidence( net12 )
querygrain( net12, nodes=c("D"))