library(MuMIn)
data<-read.csv("data.csv")

temp.mod<-lm(temp~time+dfw+time*dfw,data=data)
sal.mod<-lm(sal~time+dfw+time*dfw,data=data)
sd.mod<-lm(sdens~time+dfw+time*dfw,data=data)
lai.mod<-lm(lai~time+dfw+time*dfw,data=data)
el.mod<-lm(epload~time+dfw+time*dfw,data=data)