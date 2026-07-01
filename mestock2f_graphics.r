rm(list=ls())

#current_path = rstudioapi::getActiveDocumentContext()$path 
#setwd(dirname(current_path ))



source('read.admb.R')
source('read.admbFit.R')
source('por_recluta_r.R')
library(areaplot)


# Nombre de las flotas

nF1="Trampas"
nF2="Buceo"
name="Langosta_base12"

 #   shell('erase for_R.rep ')  # for model running
  #  system(paste('mestock2f -ind ',name,'.dat -nox',sep=""))  # for model running
  #  shell(paste('copy for_R.rep ',name,'.rep',sep=""))  # for model running


data <-read.rep(paste(name,'.rep',sep=""))

attach(data)


target=0.4

# Ajuste desembarques--------------------------------------------------------------------------
par(mfrow = c(2, 2))
sim=20

plot(Yrs,Y_obs_pred_F1[1,],cex.lab=1.5, xlab="Año",ylab="Toneladas",ylim = c(0,max(Y_obs_pred_F1[1,])*1.1), 
        main=nF1, cex.main = 1.5, major.ticks = "years",pch=sim,type="b",cex=1.5)
lines(Yrs,Y_obs_pred_F1[2,],col="red")

plot(Yrs,Y_obs_pred_F2[1,], cex.lab=1.5, xlab="Año",ylab="Toneladas",ylim = c(0,max(Y_obs_pred_F2[1,])*1.1), 
        main=nF2, cex.main = 1.5, major.ticks = "years",pch=sim,type="b",cex=1.5)
lines(Yrs,Y_obs_pred_F2[2,],col="red")


Ytot=Y_obs_pred_F1[1,]+Y_obs_pred_F2[1,]

plot(Yrs,Ytot, cex.lab=1.5, xlab="Año",ylab="Toneladas",ylim = c(0,max(Ytot)*1.1), 
     main="Total", cex.main = 1.5, major.ticks = "years",pch=sim,type="b",cex=1.5)
lines(Yrs,Y_obs_pred_F1[2,]+Y_obs_pred_F2[2,],col="red")
lines(Yrs,Y_obs_pred_F1[1,],col="blue",lwd=2)
lines(Yrs,Y_obs_pred_F2[1,],col="green",lwd=2)
legend("topleft",c(nF1,nF2,"total"),lty=c(1,1),col =c("blue","green","brown"),lwd=2,bty="n")


#Indices--------------------------------------------------------------------------

par(mfrow = c(2, 2))
sim=20

ubi=which(CPUE_obs_pred_F1[1,]>0)
max=max(CPUE_obs_pred_F1[,ubi]);
suave = smooth.spline(Yrs[ubi],CPUE_obs_pred_F1[1,ubi], spar=0.6)
cv=sd(log(CPUE_obs_pred_F1[2,ubi])-log(CPUE_obs_pred_F1[1,ubi]))
plot(Yrs[ubi],(CPUE_obs_pred_F1[1,ubi]),ylab="Indice",xlab="Año",main=paste(nF1," (cv=",round(cv,2),")"),
     pch=sim,type="b",cex=1.5,ylim=c(0,max))
lines(Yrs[ubi],(CPUE_obs_pred_F1[2,ubi]),col="red",lwd=2)
lines(suave,col="green",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

ubi=which(CPUE_obs_pred_F2[1,]>0)
max=max(CPUE_obs_pred_F2[,ubi]);
suave = smooth.spline(Yrs[ubi],CPUE_obs_pred_F2[1,ubi], spar=0.6)
cv=sd(log(CPUE_obs_pred_F2[2,ubi])-log(CPUE_obs_pred_F2[1,ubi]))
plot(Yrs[ubi],(CPUE_obs_pred_F2[1,ubi]),ylab="Indice",xlab="Año",main=paste(nF2," (cv=",round(cv,2),")"),
     pch=sim,type="b",cex=1.5,ylim=c(0,max))
lines(Yrs[ubi],(CPUE_obs_pred_F2[2,ubi]),col="red",lwd=2)
lines(suave,col="green",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)



#Talla media-------------------------------------------------------------------------

#par(mfrow = c(1, 1))
sim=20
cv_max=0.1


suave = smooth.spline(Lm_F1[1,],Lm_F1[2,], spar=0.6)
cv=sd(log(Lm_F1[3,])-log(Lm_F1[2,]))
mu=mean(Lm_F1[2,])
n=(1.96*cv*mu/(mu*cv))^2
plot(Lm_F1[1,],Lm_F1[2,],ylab="Talla promedio",xlab="Año",main=paste(nF1,"(n=",round(n,1),")"),
     pch=sim,type="b",cex=1.5,ylim=c(0.95*min(Lm_F1[2,],Lm_F1[3,]),1.05*max(Lm_F1[2,],Lm_F1[3,])))
lines(Lm_F1[1,],Lm_F1[3,],lwd=2,col="red")
lines(suave,col="green",lwd=2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)



#Comps_tallas_F1----------------------------------------------------------------------
n=length(Tallas)
filas=4
cols=3

par(mfcol = c(filas, cols))
for (i in 1:length(Nanos_frec_F1))
{
  areaplot(Tallas,pobs_F1[i,],main=paste(nF1,Nanos_frec_F1[i]),lwd=0.5, col="lightgray",ylab="", xlab="Talla",
               ylim=c(0,max(c(pobs_F1[i,],ppred_F1[i,]))))
  lines(Tallas,ppred_F1[i,],col="red",lwd=2)
}

#Comps_tallas_F2----------------------------------------------------------------------

if(length(Nanos_frec_F2)>1){
par(mfcol = c(3, 2))
for (i in 1:length(Nanos_frec_F2))
{
  areaplot(Tallas,pobs_F2[i,],main=paste(nF2,Nanos_frec_F2[i]),lwd=0.5, col="lightblue",ylab="", xlab="Talla",
           ylim=c(0,max(c(pobs_F2[i,],ppred_F2[i,]))))
  lines(Tallas,ppred_F2[i,],col="red",lwd=2)
}
} else{
  
  areaplot(Tallas,pobs_F2,main=paste(nF2,Nanos_frec_F2),lwd=0.5, col="lightblue",ylab="", xlab="Talla",
           ylim=c(0,max(c(pobs_F2,ppred_F2))))
  lines(Tallas,ppred_F2,col="red",lwd=2)

}

#Marginales---------------------
par(mfrow = c(1, 2))

sd1=sd(Frec_marg_F1[1,]-Frec_marg_F1[2,])
sd2=sd(Frec_marg_F2[1,]-Frec_marg_F2[2,])


areaplot(Tallas,Frec_marg_F1[1,],main=paste(nF1),lwd=0.5, col="gray",ylab="", xlab="Talla",
         ylim=c(0,max(Frec_marg_F1)))
#lines(Tallas,Frec_marg_F1[2,],col="red",lwd=2)

areaplot(Tallas,Frec_marg_F2[1,],main=paste(nF2),lwd=0.5, col="lightblue",ylab="", xlab="Talla",
         ylim=c(0,max(Frec_marg_F2)))
#lines(Tallas,Frec_marg_F2[2,],col="red",lwd=2)

  
  
#Crecimiento------------
par(mfrow = c(2, 1))
nedades=length(N_tot[1,])

plot(Tallas,No[1]*Prob_talla[1,],type="l",xlab="Talla",ylab="Densidad",main="Componentes de edades")

for (i in 2:nedades)
{ 
  lines(Tallas,No[i]*Prob_talla[i,])
}

abline(v=Lmed_edad,lty=2, col="gray")
lines(Tallas,No[1]*Prob_talla[1,],col="red",lwd=2)
text(round(Lmed_edad[1],2),100,paste(round(Lmed_edad[1],2)),col="blue")


edades=seq(min_edad,min_edad+length(Lmed_edad)-1)
plot(edades,Lmed_edad,type="b", xlab="Edad relativa",ylab="Talla",main="Crecimiento individual",
     ylim=c(0,max(Lmed_edad))*1.05)
abline(v=Lmed_edad[1:nedades-1],lty=2, col="gray")
abline(h=Lmed_edad[2:nedades],lty=2, col="gray")


#Selectividad----------------------------------------------------------------------
par(mfrow = c(2, 2))
persp(Yrs,Lmed_edad,(Sel_F1),xlab="Año",ylab="Talla",zlab="Proporcion",main=nF1,phi = 30, theta=60, expand=0.5,col="lightblue")
persp(Yrs,Lmed_edad,(Sel_F2),xlab="Año",ylab="Talla",zlab="Proporcion",main=nF2,phi = 30, theta=60, expand=0.5, col="lightblue")

matplot(Lmed_edad,t(Sel_F1),xlab="Talla",ylab="Talla",main=nF1,type="l",lty=1,ylim=c(0,1),col=1)
lines(Lmed_edad,Msex_edad,col="red")
matplot(Lmed_edad,t(Sel_F2),xlab="Talla",ylab="Talla",main=nF2,type="l",lty=1,col=1,ylim=c(0,1))
lines(Lmed_edad,Msex_edad,col="red")


#Poblacionales----------------------------------------------------------------------
par(mfrow = c(2, 2))

plot(Yrs,BT,type="l",xlab="Años",ylab="Biomasa",ylim=c(0, max(BT)),lwd=2,col="black")#, cex.lab=1.2,cex.axis=1.2)
lines(Yrs,SSB,col="green",lwd=3)
legend("topright",c("Total", "Desovante"),lty =c(1,1),col=c("black","green"),lwd=2,cex=1.0,bty = "n")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)


sdR=R_sd
li=Reclutas[1,]-1.96*sdR
ls=Reclutas[1,]+1.96*sdR
nyrs=length(SSB)
plot(Yrs,Reclutas[1,],ylim = c(0,max(Reclutas[1,]+1.96*sdR)*1.01),type="l",ylab="Reclutamientos",xlab="Año",pch = 16,cex=1,lwd=2)#,cex.lab=1.2,cex.main=2,cex.axis=1.2)
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])
polygon(x,y,col="#DCDCDC",border="#DCDCDC")
lines(Yrs,Reclutas[1,],lwd="2")
lines(Yrs,Reclutas[2,],lwd="2",col="red")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

plot(Yrs,Reclutas[3,],type="l",xlab="Años",ylab="Desvio log_R",lwd=2,col="black")#,cex.lab=1.2,cex.main=2,cex.axis=1.2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
abline(h=0,col="red",lwd = 2)

#plot(Yrs,Ftot,type="l",xlab="Años",ylab="Mortalidad por pesca",lwd=2,col="black",cex.lab=1.2,cex.axis=1.2)
#grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)


#Stock-recluta--------------------------------------
alfa=4*h*R0/(5*h-1)
beta=(1-h)*B0/(5*h-1)

ssb=seq(0,max(c(B0,max(SSB))),B0/50)
rec=alfa*ssb/(beta+ssb)

plot(SSB[seq(1,length(Yrs)-min_edad)],Reclutas[1,seq(min_edad+1,length(Yrs))],
     ylim=c(0,max(Reclutas[1,])),type="b",
     xlim=c(0,max(c(B0,max(SSB)))),
     xlab="Biomasa desovante",
     ylab="Reclutamiento")
lines(ssb,rec,col="blue",lwd=2)
text(SSB[seq(1,length(Yrs)-min_edad)],Reclutas[1,seq(min_edad+1,length(Yrs))],
     paste(Yrs[seq(min_edad+1,length(Yrs))]),cex=0.8)

text(SSB[length(Yrs)-min_edad],Reclutas[1,length(Yrs)],
     paste(Yrs[length(Yrs)]),cex=0.8,col="red")


#Analisis por recluta------------

tmax=length(Msex_edad)
target=0.4
ypr_out<-por_recluta_r(tmax,Sel_tot[nyrs,],Msex_edad,Peso_edad,target,h,M,dts)
attach(ypr_out)
Ftar=Ftar
BPRtar=BPRtar
YPRtar=YPRtar
RMS=YPRtar*R0

ubi=max(which(B>0))

par(mfrow = c(1, 1))
plot(Fcr,Y/max(Y),type="l", col="green", lwd=2, main="Analisis por recluta", xlab="Mortalidad por pesca", ylab="BPR, YPR relativos",
     cex.lab=1.1,cex.main=1.5,ylim = c(0,1),xlim=c(0,Fcr[ubi]))
lines(Fcr,B/max(B), col="red", lwd=2)
lines(Ftar, BPRtar/max(B),type="p", lwd=5)
text(Ftar, 0,Ftar,cex=1.1)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
abline(h = target, lty = 2,lwd=1)
abline(v = Ftar, lty = 2,lwd=1)


Ftar_vec=rep(0,length(Yrs))
for (i in 1:length(Yrs))
{ 
  ypr_out<-por_recluta_r(tmax,Sel_tot[i,],Msex_edad,Peso_edad,target,h,M,dts)
  Ftar_vec[i]=ypr_out$Ftar
}

#Biomasa y F con IC-----------------------------------------------------------
par(mfrow = c(2, 2))

li=SSB-1.96*BD_sd
ls=SSB+1.96*BD_sd

plot(Yrs,SSB,ylim = c(0,max(ls)*1.01),type="l",
     ylab="Biomasa",xlab="Año",pch = 16,cex=1,lwd=2,main="Biomasa") #, cex.lab=1.5,cex.main=1.5)
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])

polygon(x,y,col="#DCDCDC",border="#DCDCDC")  

lines(Yrs,SSB,lwd=2)
abline(h = target*B0, col = "red",lty = 2,lwd=2)



li=Ftot-1.96*F_sd
ls=Ftot+1.96*F_sd

plot(Yrs,Ftot,ylim = c(0,max(ls)*1.01),type="l",
     ylab="F",xlab="Año",pch = 16,cex=1,lwd=2,main="Mortalidad por pesca")#,cex.lab=1.5,cex.main=1.5)
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])

polygon(x,y,col="#DCDCDC",border="#DCDCDC")  

lines(Yrs,Ftot,lwd=2)
lines(Yrs,Ftar_vec, col = "red",lty = 2,lwd=2)

for (i in 1:length(Yrs))
{ 
  F_F1[i]=max(Fa_F1[i,])
  F_F2[i]=max(Fa_F2[i,])
}



lines(Yrs,F_F1, col = "blue")
lines(Yrs,F_F2, col = "green")
legend("topleft",c("Total",nF1,nF2),lty =c(1,1,1),col=c("black","blue","green","brown"),lwd=c(2,1,1,1),cex=1.0,bty = "n")



li=B_B0-1.96*SPR_sd
ls=B_B0+1.96*SPR_sd

plot(Yrs,B_B0,ylim = c(0,max(ls)*1.01),type="l",
     ylab="B/B0",xlab="Año",pch = 16,cex=1,lwd=2,main="B/B0")#,cex.lab=1.5,cex.main=1.5)
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])

polygon(x,y,col="#DCDCDC",border="#DCDCDC")  

lines(Yrs,B_B0,lwd=2)
abline(h = target, col = "red",lty = 2,lwd=2)



#Kobe---------------------------------------------------------------------
par(mfrow = c(1, 1))
BRMS=B0*target
FRMS=Ftar
B_Btar=SSB/BRMS
F_Ftar=Ftot/Ftar_vec

SPR=B_B0
p_low=pnorm(target,SPR[nyrs],SPR_sd[nyrs])
p_high=1-pnorm(Ftar,Ftot[nyrs],F_sd[nyrs])


plot(B_Btar,F_Ftar,pch = 16,ylab="F/Frms",xlab="B/Brms",xlim = c(0,max(B_Btar)), ylim = c(0,max(F_Ftar)*1.5), 
     type="o",col="black",lty="dashed",main=paste("B/Brms=",round(B_Btar[nyrs],2),"(risk=",round(p_low,2),")",
                                                  " F/Frms=",round(F_Ftar[nyrs],2),"(risk=",round(p_high,2),")"))

polygon(c(0,1,1,0),c(0,0,1,1),col="yellow1") #amarillo
polygon(c(1,1.1*max(B_Btar),1.1*max(B_Btar),1),c(0,0,1,1),col="green") #verde
polygon(c(1,1.1*max(B_Btar),1.1*max(B_Btar),1),c(1,1,1.5*max(F_Ftar),1.5*max(F_Ftar)),col="orange") #amarillo
polygon(c(0,1,1,0),c(1,1,1.5*max(F_Ftar),1.5*max(F_Ftar)),col="tomato1") #rojo

lines(B_Btar,F_Ftar,pch = 16, type="o",col="black",lty="dashed")
lines(B_Btar[nyrs],F_Ftar[nyrs],type="p",col="blue",pch = 16,cex=2)
text(B_Btar*.95,F_Ftar,paste(Yrs),cex=0.8)


X0=B_Btar[nyrs]
Y0=F_Ftar[nyrs]
cvF=F_sd[nyrs]/Ftot[nyrs];
cvB=BD_sd[nyrs]/SSB[nyrs];

arrows(X0,Y0-1.96*cvF*Y0,X0,Y0+1.96*cvF*Y0,
       length = 0.05, code = 3, angle = 90, lwd=2, col="blue")

arrows(X0-1.96*cvB*X0,Y0,X0+1.96*cvB*X0,Y0,
       length = 0.05, code = 3, angle = 90,lwd=2, col="blue")


box()



#Proyecciones---------------------------------
par(mfcol = c(2, 2))

nanos=length(Yrs)
nysim=length(Bio_proy[,1])
yproy=seq(Yrs[length(Yrs)]+1,Yrs[length(Yrs)]+nysim)
vecto=seq(nyrs-9,nyrs)

plot(Yrs[vecto],SSB[vecto]/B0,type="l",xlim=c(min(Yrs[vecto]),max(yproy)),
     ylim=c(0,max(Bio_proy/B0)),lwd="2",lty=1,xlab="Año", ylab="B/B0",
     main="Biomasa")
abline(h=target,  col = "red",lty = 2,lwd=2)
abline(v=Yrs[length(Yrs)]+1,lty = 2,lwd=1)
matlines(yproy,Bio_proy/B0,lwd="2",lty=1)
x=yproy[nysim]
y=Bio_proy[nysim,]*0.98/B0
text(x,y,paste(round(Mult_F,2)),cex=1.1)

Desembarques=Y_obs_pred_F1[1,]+Y_obs_pred_F2[1,]
plot(Yrs[vecto],Desembarques[vecto],type="l",xlim=c(min(Yrs[vecto]),max(yproy)),
     ylim=c(0,max(Desembarques)),lwd="2",lty=1,xlab="Año", ylab="Capturas",
     main="Capturas")
matlines(yproy,Capt_proy,lwd=2,lty=1)
x=yproy[nysim]
y=Capt_proy[nysim,]*0.98
text(x,y,paste(round(Mult_F,2)),cex=1.1)
abline(v=Yrs[length(Yrs)]+1,lty = 2,lwd=1)
abline(h=RMS,col = "red",lty = 2,lwd=2)

plot(Yrs[vecto],Ftot[vecto],type="l",xlim=c(min(Yrs[vecto]),max(yproy)),
     ylim=c(0,c(max(max(Ftot),max(F_proy)))),lwd="2",lty=1,xlab="Año", ylab="F",
     main="Mort por pesca")

matlines(yproy, F_proy,type="l",lwd=2,lty=1)
abline(h=Ftar_vec[length(Yrs)],  col = "red",lty = 2,lwd=2)
abline(v=Yrs[length(Yrs)]+1,lty = 2,lwd=1)

y=F_proy[nysim,]*0.98
text(x,y,paste(round(Mult_F,2)),cex=1.1)



#Riesgos----------------------------------------------------------------------------
par(mfrow = c(1, 2))
barplot(Red_stock[1,]~Mult_F,ylim=c(0,1.1),xlab="Mult del esfuerzo", ylab="B/B0", main="Reduccion del stock")
abline(h=target,  col = "red",lty = 2,lwd=2)
box()


riesgo=1-pnorm(Red_stock[1,],target,Red_stock[2,])
riesgo_crash=1-pnorm(Red_stock[1,],0.5*target,Red_stock[2,])

plot(Mult_F,riesgo,type = "p",ylim=c(0,1.1),xlab="Mult del esfuerzo", ylab="p(B<Brms)", main="Riesgo largo plazo")
lines(Mult_F,riesgo,col="red",lwd=2)
text(Mult_F*1.02,riesgo,round(riesgo,3))
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

box()


#Genera excel--------------------------------------------------------------------------------------

ubi=which(Yrs>2000)

Variables=data.frame(Año=Yrs[ubi],Biomasa=SSB[ubi],R_R0=Reclutas[2,ubi]/R0,
                     Fcr=Ftot[ubi], F_Fmrs=Ftot[ubi]/Ftar_vec[ubi],B_Brms=SSB[ubi]/BRMS,
                     B_B0=SPR[ubi])
write.csv(round(Variables,2), paste('Var_Pobl_',name,'.csv'))


Variables2=data.frame(Mult_Eff=Mult_F,B_B0=Red_stock[1,],Riesgo=riesgo,Riesgo_colapso=riesgo_crash,Captura_cp=Capt_proy[1,], 
                      Captura_lp=colMeans(Capt_proy[seq(nysim,nysim-4,-1),]))


write.csv(Variables2, paste('Manejo_',name,'.csv'))






