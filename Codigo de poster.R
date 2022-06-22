library(MASS)
library(rstanarm)
library(insight)
library(MCMCpack)
library(compiler)
library(parallel)
library(tidyr)
library(foreach)
library(doParallel)
library(xtable)
library(ggplot2)
library(gridExtra)

f<-function(m,n,sigma,B0,B1){
  b01=numeric()
  b11=numeric()
  b02=numeric()
  b12=numeric()
  b03=numeric()
  b13=numeric()
  con1=numeric()
  con2=numeric()
  con3=numeric()
  x=runif(n,0,4)
  
  for (i in 1:m) {
    e=rnorm(n,0,sigma)
    y=B0+B1*x+e #Parametro B0=3,B1=3
    mod<-lm(y~x)
    b01[i]<-mod$coefficients[1]
    b11[i]<-mod$coefficients[2]
    mod1<-rlm(y~x)
    b02[i]<-mod1$coefficients[1]
    b12[i]<-mod1$coefficients[2]
    mod2<-summary(MCMCregress(y~x))$statistics
    b03[i]<-mod2[1]
    b13[i]<-mod2[2]
    
    a<-min(abs(b01[i]+b11[i]-(B0+B1)),
           abs(b02[i]+b12[i]-(B0+B1)),
           abs(b03[i]+b13[i]-(B0+B1)))
    
    ifelse(a==abs(b01[i]+b11[i]-(B0+B1)),con1[i]<-1,con1[i]<-0)
    ifelse(a==abs(b02[i]+b12[i]-(B0+B1)),con2[i]<-1,con2[i]<-0)
    ifelse(a==abs(b03[i]+b13[i]-(B0+B1)),con3[i]<-1,con3[i]<-0)
    
  }
  
  b0_pos =c(mean(tail(sort(b01),m*0.005)),mean(tail(sort(b02),m*0.005)),
            mean(tail(sort(b03),m*0.005)))
  b0_neg =c(mean(head(sort(b01),ceiling(m*0.005))),
            mean(head(sort(b02),ceiling(m*0.005))),
            mean(head(sort(b03),ceiling(m*0.005))))
  b1_pos=c(mean(tail(sort(b11),m*0.005)),mean(tail(sort(b12),m*0.005)),
           mean(tail(sort(b13),m*0.005)))
  b1_neg=c(mean(head(sort(b11),ceiling(m*0.005))),
           mean(head(sort(b12),ceiling(m*0.005))),
           mean(head(sort(b13),ceiling(m*0.005))))
  
  e01=mean(tail(sort(abs(b01-B0)),m*0.01))
  e11=mean(tail(sort(abs(b11-B1)),m*0.01))
  e02=mean(tail(sort(abs(b02-B0)),m*0.01))
  e12=mean(tail(sort(abs(b12-B1)),m*0.01))
  e03=mean(tail(sort(abs(b03-B0)),m*0.01))
  e13=mean(tail(sort(abs(b13-B1)),m*0.01))
  
  e0=c(e01,e02,e03)
  e1=c(e11,e12,e13)
  et=e0+e1
  

  
  prop=c(mean(con1),mean(con2),mean(con3))
  
  
  Resultados<-data.frame(Metodo=c("lm","rlm","MCMC"),Error_b0=e0,
                         Error_b1=e1,Error_total=et,b0_neg,b0_pos,
                         b1_neg,b1_pos)
  
  Resultados1<-data.frame(n,sigma,Metodo=c("lm","rlm","MCMC"),Error_b0=e0,
                         Error_b1=e1,Error_total=et,prop)

  return(Resultados1)  
  
}
vs_model<-cmpfun(f)


vs_model(m=1000,n=100,sigma=1,B0=1,B1=1)  #una simulación (puede tardar algunos minutos)



registerDoParallel(makeCluster(detectCores() - 1)) #computo paralelo

res<-foreach(m=rep(100,15),n=c(5,5,5,10,10,10,20,20,20,40,40,40,100,100,100),
        sigma=c(1,3,5,1,3,5,1,3,5,1,3,5,1,3,5),
        B0=rep(1,15),B1=rep(1,15))%do%vs_model(m,n,sigma,B0,B1)

stopImplicitCluster()

res

resultados<-do.call(rbind.data.frame, res) 
resultados<-resultados[,-c(3,4)]

resultados[1]<-as.integer(resultados$n)
resultados[2]<-as.integer(resultados$sigma)

resultados  ####resultados de la simulacion


######Graficas 
##Este proceso puede tardar algunos minutos, incluso horas


graf<-foreach(m=rep(100,288),n=rep(seq(5,100),3),
              sigma=c(rep(1,96),rep(3,96),rep(5,96)),
              B0=rep(1,288),B1=rep(1,288))%do%f(m,n,sigma,B0,B1)

graf
grafica<-do.call(rbind.data.frame, graf) 

grafica

grafica<-grafica[,-c(4,5,7)]
graf1<-grafica[grafica$sigma==1,-2]
graf1
graf2<-grafica[grafica$sigma==3,-2]
graf2
graf3<-grafica[grafica$sigma==5,-2]
graf3


p1<-graf1 %>%
  ggplot() +
  aes(x = n, y = Error_total,label=Metodo, color=Metodo)+
  geom_line()+
  xlab("n") +
  ylab(" Error total promedio del 10% \n de las peores estimaciones") +
  ggtitle("Datos para sigma = 1") +
  geom_hline(yintercept = 0)


p2<-graf2 %>%
  ggplot() +
  aes(x = n, y = Error_total,label=Metodo, color=Metodo)+
  geom_line()+
  xlab("n") +
  ylab(" Error total promedio del 10% \n de las peores estimaciones") +
  ggtitle("Datos para sigma = 3") +
  geom_hline(yintercept = 0)


p3<-graf3 %>%
  ggplot() +
  aes(x = n, y = Error_total,label=Metodo, color=Metodo)+
  geom_line(size=0.5)+
  xlab("n") +
  ylab(" Error total promedio del 10% \n de las peores estimaciones") +
  ggtitle("Datos para sigma = 5") +
  geom_hline(yintercept = 0)+
  theme (plot.title = element_text(size=rel(1.2)))+
  theme(axis.title.x = element_text( size=rel(1))) +
  theme(axis.title.y = element_text( size=rel(1)))
  

f1<-grid.arrange(p1,p2,p3,nrow=3)


 

