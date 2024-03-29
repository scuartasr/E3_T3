# Eliminaci?n de todos los objetos anteriores
rm(list=ls(all=TRUE))

# _____________________________________________________________________________
# _____________________________________________________________________________
# Paquetes y funciones necesarias =============================================


# Paquetes

library(forecast)
library(TSA)
library(car)
library(lmtest)
library(FitAR)
library(uroot)
library(pdR)

# Funciones de usuario

source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SuavizamientoEstacional.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Descomp.Loess.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-BP.LB.test-pruebaDW1.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexpo.ErrorARMA.R")
##
##
##

# _____________________________________________________________________________
# _____________________________________________________________________________
# 1. Lectura de datos =========================================================

# Lectura de la base de datos
#enlace <- "./anexos-emmet-noviembre-2021-1-total industria-modif.csv"

datos <- read.table(file = file.choose(),
                    header=T,
                    sep=";",
                    skip=15,
                    dec=",")
# Selecci?n de la columna necesaria
datos <- as.data.frame(datos[, 7])
colnames(datos) <- c("ventas.nominal")

# Creaci?n del objeto ts
datos <- ts(datos,
            freq=12,
            start=c(2001,1))

##
##
##

# _____________________________________________________________________________
# _____________________________________________________________________________
# 2. An?lisis descriptivo =====================================================

# Gr?fica de la serie de tiempo lineal

par(adj = 0.5, col = 'black')
plot(datos,
     xlab = "A?o\n",
     ylab = "?ndice de ventas nominal",
     cex.main = 1,
     lwd = 1)
grid(col = 'gray', lwd = 1)
par(adj = 1,
    col = 'black')

# Gr?fica de la serie de tiempo en escala logar?tmica

par(adj = 0.5, col = 'black')
plot(log(datos),
     xlab = "A?o\n",
     ylab = "Logaritmo del ?ndice de ventas nominal",
     cex.main = 1,
     lwd = 1)
grid(col = 'gray', lwd = 1)
par(adj = 1,
    col = 'black')

# Barplot para la estacionalidad

par(adj = 0.5, col = 'black')
boxplot(log(datos)~cycle(log(datos)),
        names = month.abb,
        cex.axis = 0.8,
        main = "Valor del ?ndice mensual de ventas nominal
        en Colombia por mes",
        xlab = "Mes",
        ylab = "?ndice mensual de ventas nominal")
par(adj = 1)
title(sub = "Fuente. DANE")

# Periodograma

par(adj = 0.5)
log.desest <- diff(log(datos))
periodogram(log.desest,
            lwd = 2,
            main = "Periodograma para el ?ndice de ventas nominal
            en Colombia por mes",
            xlab = "Frecuencias",
            ylab = "Asociaci?n")
abline(v = c(1:6)/12,
       col=2, lty=2)

# Gr?fico de la componente de tendencia de la gr?fica

par(adj = 0.5)
plot(decompose(log(datos))$trend,
     xlab = "A?o",
     ylab = "Logaritmo del ?ndice mensual de ventas nominales",
     main = "Tendencia del logaritmo ?ndice mensual de ventas
     nominales en Colombia entre 2001 y 2021",
     lwd = 2)
par(adj = 1,
    col = 'black')
title(sub = "Fuente. DANE. Estimaci?n de la tendencia hecha con R.")
grid(col = 'gray', lwd = 1)

# Gr?fico de la funci?n de autocorrelaci?n muestral (ACF) del logaritmo

par(adj = 0.5)
acf(as.numeric(log(datos)),
    xlab = 'Rezago',
    lag.max = 36,
    ci.type = "ma",
    col = 1,
    ci.col = 2,
    main = '')

# Gr?fico de la componente estacional de la serie

descomposicion <- decompose(datos, type = "multiplicative")

par(adj = 0.5)
plot(descomposicion$seasonal,
     main = "Componente estacional del ?ndice mensual de
     ventas nominales en Colombia entre 2001 y 2021",
     xlab = "A?o",
     ylab = "?ndice mensual de ventas nominales")
par(adj = 1,
    col = 'black')
title(sub = "Fuente. DANE. Estimaci?n de la estacionalidad hecha con R.")

par(adj = 1,
    col = 'black')

##
##
##

# _____________________________________________________________________________
# _____________________________________________________________________________
# 3. Preparac?n para el modelo exponencial polinomial estacional ==============

# Definici?n de variables necesarias

m <- 12                # N?mero de periodos a pronosticar dentro de la muestra
n <- length(datos) - m # Tama?o de la muestra para el ajuste
t <- 1:n               # ?ndice de tiempo en los periodos de ajuste
t2 <- t^2
t3 <- t^3
t4 <- t^4
t5 <- t^5
t6 <- t^6

# Funciones trigonom?tricas para la componente estacional

sen1 <- sin(pi*t/6)
cos1 <- cos(pi*t/6)
sen2 <- sin(pi*t/3)
cos2 <- cos(pi*t/3)
sen3 <- sin(pi*t/2)
cos3 <- cos(pi*t/2)
sen4 <- sin(2*pi*t/3)
cos4 <- cos(2*pi*t/3)
sen5 <- sin(5*pi*t/6)
cos5 <- cos(5*pi*t/6)

# Serie de tiempo con los n datos para construir el modelo

yt <- ts(datos[t],
         frequency = 12,
         start=c(2001,1))

# Matriz de dise?o

X1 <- data.frame(t, t2, t3, t4, t5, t6,
                 sen1, cos1, sen2, cos2,
                 sen3, cos3, sen4, cos4,
                 sen5, cos5)

# Valores de las variables en la validaci?n cruzada

tnuevo <- (n+1):length(datos) # ?ndice de tiempo en los pron?sticos
t2nuevo <- tnuevo^2
t3nuevo <- tnuevo^3
t4nuevo <-  tnuevo^4
t5nuevo <-  tnuevo^5
t6nuevo <-  tnuevo^6

#Funciones trigonom?tricas en los periodos de pron?stico

sen1n <- sin(pi*tnuevo/6)
cos1n <- cos(pi*tnuevo/6)
sen2n <- sin(pi*tnuevo/3)
cos2n <- cos(pi*tnuevo/3)
sen3n <- sin(pi*tnuevo/2)
cos3n <- cos(pi*tnuevo/2)
sen4n <- sin(2*pi*tnuevo/3)
cos4n <- cos(2*pi*tnuevo/3)
sen5n <- sin(5*pi*tnuevo/6)
cos5n <- cos(5*pi*tnuevo/6)

# Serie de tiempo para los valores que van a ser usados para el ajuste

ytf <- ts(datos[tnuevo],
          freq = 12,
          start = c(2020, 12))

# Matriz de dise?o con las potencias de t en polinomio de grado seis
# y trigonom?tricas para la estacionalidad

X1nuevo <- data.frame(t = tnuevo, t2 = t2nuevo, t3 = t3nuevo,
                      t4 = t4nuevo, t5 = t5nuevo, t6 = t6nuevo,
                      sen1 = sen1n, cos1 = cos1n, sen2 = sen2n,
                      cos2 = cos2n, sen3 = sen3n, cos3 = cos3n,
                      sen4 = sen4n, cos4 = cos4n, sen5 = sen5n,
                      cos5=cos5n)

##
##
##

# _____________________________________________________________________________
# _____________________________________________________________________________
# 4. Gr?ficos de las series diferenciadas y sus ACFs ==========================

# Logserie

par(adj = 0.5)
plot(log(yt),
     xlab = "A?o\n",
     ylab = "Logaritmo del ?ndice de ventas nominal",
     cex.main = 1,
     lwd = 1)
grid(col = 'gray', lwd = 1)
par(adj = 1,
    col = 'black')

par(adj = 0.5)
acf(as.numeric(log(yt)),
    xlab = 'Rezago',
    lag.max = 36,
    ci.type = "ma",
    col = 1,
    ci.col = 2,
    main=expression(paste("ACF",sep=" ","de",sep=" ",log(Y[t]))))

## __________________________________________
## 4.1. Diferencias sobre la serie =========
## __________________________________________

d = 1
D = 1

# Serie logar?tmica y su ACF

plot(diff(log(yt)),
     ylab = expression(paste(nabla,sep="",log(Y[t]))),
     xlab = "A?o")
abline(h=mean(diff(log(yt))),
       col = 2)
grid(col = 'gray', lwd = 1)

### Primera diferencia regular (???logY???)

lny=ts(log(datos[t]),freq=12,start=c(2001,1)) # Log. de la serie recortada a
# ajustar y su primera diferencia
diff1lny=diff(lny)                            # Primera diferencia regular

plot(diff1lny,
     ylab = expression(paste(nabla,sep="",log(Y[t]))),
     xlab = "A?o")
abline(h=mean(diff(log(yt))),
       col = 2)
grid(col = 'gray', lwd = 1)

acf(as.numeric(diff1lny),
    xlab = 'Rezago',
    lag.max=36,
    ci.type="ma",
    ci.col="red",
    main=expression(paste("ACF",sep=" ","de",sep=" ",nabla,sep="",log(Y[t]))))

# Primera diferencia estacional (?????????logY???)

diff12lny = diff(lny,
                 lag = 12)

plot(diff12lny,
     ylab = expression(paste(nabla[12],sep="",log(Y[t]))),
     xlab = "A?o")
abline(h=mean(diff12lny),
       col = 2)
grid(col = 'gray', lwd = 1)

acf(as.numeric(diff12lny),
    xlab = 'Rezago',
    lag.max=36,
    ci.type="ma",
    ci.col="red",
    main=expression(paste("ACF",sep=" ","de",
                          sep=" ",nabla[12],sep="",log(Y[t]))))
grid(col = 'gray', lwd = 1)

# Diferencia mixta. Se agrega su PACF.

diffmixta = diff(diff(lny, lag = 12))

plot(diffmixta,
     ylab = expression(paste(nabla,nabla[12],sep="",log(Y[t]))),
     xlab = "A?o")
abline(h=mean(diffmixta),
       col = 2)
grid(col = 'gray', lwd = 1)

acf(as.numeric(diffmixta),
    xlab = 'Rezago',
    lag.max=36,
    ci.type="ma",
    ci.col="red",
    main=expression(paste("ACF",sep=" ","de",sep=" ",
                          nabla,nabla[12],sep="",log(Y[t]))))
grid(col = 'gray', lwd = 1)

pacf(as.numeric(diffmixta),
     xlab = 'Rezago',
     lag.max=36,
     ci.col="red",
     main=expression(paste("ACF",sep=" ","de",sep=" ",
                           nabla,nabla[12],sep="",log(Y[t]))))

##
##
##

# _____________________________________________________________________________
# _____________________________________________________________________________
# 5. Test HEGY  ---------------------------------------------------------------

HEGY.test(wts = lny,
          itsd = c(0, 0, c(0)),
          selectlags = list(mode = "aic", Pmax = 12))$stats
##
##
##

# _____________________________________________________________________________
# _____________________________________________________________________________
# 6. Identificaci?n de m?todos SARIMA =========================================

## __________________________________________
## 6.1. Usando auto.arima ===================
## __________________________________________

auto.arima(lny,ic="aic",seasonal.test="ocsb")
auto.arima(lny,ic="aic",seasonal.test="ch")
auto.arima(lny,ic="aic",seasonal.test="seas")
auto.arima(lny,ic="bic",seasonal.test="ocsb")
auto.arima(lny,ic="bic",seasonal.test="ch")
auto.arima(lny,ic="bic",seasonal.test="seas")

## __________________________________________
## 6.2. Usando armasubsets ==================
## __________________________________________

win.graph(heigh=5,width=9) 
plot(armasubsets(diffmixta,nar=12,nma=12,y.name='AR',ar.method='ols'))

win.graph(heigh=5,width=9) 
plot(armasubsets(diffmixta,nar=18,nma=18,y.name='AR',ar.method="ols"))

##
##
##
##

# _____________________________________________________________________________
# _____________________________________________________________________________
# 7. Ajuste de modelos ========================================================

# Modelo uno. ARIMA(2, 1, 0)(0, 1, 2)[12]

modelo1 <- Arima(lny,
                 order = c(2, 1, 0),
                 seasonal = list(order = c(0, 1, 2)),
                 method = "ML")
modelo1
coeftest(modelo1)

# Modelo dos. ARMA(4, 1, 0)(1, 1, 2)[12]

modelo2 <- Arima(lny,
                 order = c(4, 1, 0),
                 seasonal = list(order = c(1, 1, 2)),
                 method = 'ML')
modelo2
coeftest(modelo2)

# Modelo tres. ARMA(6, 1, 10)(0, 1, 1)[12]

p3 <-  c(NA, NA, 0, 0, 0, NA)
q3 <-  c(rep(0, 9), NA)
Q3 <- c(NA)
arma3 <- c(p3, q3, Q3)
modelo3 <- Arima(lny,
                 order = c(6, 1, 10),
                 seasonal = list(order = c(0, 1, 1)),
                 method = 'ML',
                 fixed = arma3)
modelo3
coeftest(modelo3)

# Modelo cuatro. 

p4 <- c(NA, NA, rep(0, 6), NA)
q4 <- c(rep(0, 9), NA)
Q4 <- c(NA)
arma4 <- c(p4, q4, Q4)

modelo4 <- Arima(lny,
                 order = c(9, 1, 10),
                 seasonal = list(order = c(0, 1, 1)),
                 method = 'ML',
                 fixed = arma4)
modelo4
coeftest(modelo4)


#C?lculo de AIC y BIC versi?n exp(Cn*(p))
yhat1=exp(modelo1$fitted)*exp(modelo1$sigma2/2)
yhat2=exp(modelo2$fitted)*exp(modelo2$sigma2/2)
yhat3=exp(modelo3$fitted)*exp(modelo3$sigma2/2)
yhat4=exp(modelo4$fitted)*exp(modelo4$sigma2/2)

#seudo residuos
resorig1=yt-yhat1
resorig2=yt-yhat2
resorig3=yt-yhat3
resorig4=yt-yhat4


k1=length(coef(modelo1)[coef(modelo1)!=0]);k1 #n?mero de par?metros del modelo1
k2=length(coef(modelo2)[coef(modelo2)!=0]);k2 #n?mero de par?metros del modelo2
k3=length(coef(modelo3)[coef(modelo3)!=0]);k3 #n?mero de par?metros del modelo3
k4=length(coef(modelo4)[coef(modelo4)!=0]);k4 #n?mero de par?metros del modelo4


Criterios1=exp.crit.inf.resid(resorig1,n.par=k1);Criterios1
Criterios2=exp.crit.inf.resid(resorig2,n.par=k2);Criterios2
Criterios3=exp.crit.inf.resid(resorig3,n.par=k3);Criterios3
Criterios4=exp.crit.inf.resid(resorig4,n.par=k4);Criterios4


#Tabla n?mero de par?metros y criterios de informaci?n
tablacriterios=rbind(Criterios1,Criterios2,Criterios3,Criterios4)
rownames(tablacriterios)=paste0("Modelo",1:4)
tablacriterios


#Tabla pronosticos y los IP

predmod1=exp(as.data.frame(forecast(modelo1,h=12,level=95)))*exp(modelo1$sigma2/2) 
predmod1=ts(predmod1,freq=12,start=c(2020,12)); predmod1
Amplcobmodelo1=amplitud.cobertura(real=ytf,LIP=predmod1[,2],LSP=predmod1[,3])

predmod2=exp(as.data.frame(forecast(modelo2,h=12,level=95)))*exp(modelo2$sigma2/2) 
predmod2=ts(predmod2,freq=12,start=c(2020,12))
Amplcobmodelo2=amplitud.cobertura(real=ytf,LIP=predmod2[,2],LSP=predmod2[,3])

predmod3=exp(as.data.frame(forecast(modelo3,h=12,level=95)))*exp(modelo3$sigma2/2) 
predmod3=ts(predmod3,freq=12,start=c(2020,12))
Amplcobmodelo3=amplitud.cobertura(real=ytf,LIP=predmod3[,2],LSP=predmod3[,3])

predmod4=exp(as.data.frame(forecast(modelo4,h=12,level=95)))*exp(modelo4$sigma2/2) 
predmod4=ts(predmod4,freq=12,start=c(2020,12))
Amplcobmodelo4=amplitud.cobertura(real=ytf,LIP=predmod4[,2],LSP=predmod4[,3])


win.graph()

vector.auxiliar <- c(ytf, predmod1[,1], predmod2[,1],
                     predmod3[,1], predmod4[,1])
par(adj = 0.5)
plot(ytf, type = "b", pch = 19, lty = 1, col = 1, lwd = 2,
     ylab = "?ndice de ventas nominales",
     xlab = "Periodo [mmm - yy]",
     ylim = c(min(vector.auxiliar), max(vector.auxiliar)),
     xaxt = "n")
lines(predmod1[,1], col = 2, pch = 2, lty = 2, type = "b", lwd = 2)
lines(predmod2[,1], col = 3, pch = 3, lty = 3, type = "b", lwd = 2)
lines(predmod3[,1], col = 4, pch = 4, lty = 4, type = "b", lwd = 2)
lines(predmod4[,1], col = 5, pch = 5, lty = 5, type = "b", lwd = 2)
legend("bottomright",
       legend = c("Real", "Modelo1", "Modelo 2",
                  "Modelo 3", "Modelo 4" ), 
       col = 1:5,
       pch = c(19, 2:5),
       lty = 1:5,
       lwd = 1,
       cex= 0.7)
axis(1,at = time(ytf), 
     labels = c("dic-20", "ene-21", "feb-21", "mar-21",
                "abr-21", "may-21", "jun-21", "jul-21",
                "ago-21", "sep-21", "oct-21", "nov-21"))

#GRAFICO COMPARATIVA PRON PUNTUAL 

#mejor modelo global

mmgl <- lm(log(yt)~., data = X1)
mmgl <- exp(predict(mmgl,
                           newdata = X1nuevo,
                           interval = "prediction",
                           level = 0.95)) * exp(summary(modelo1)$sigma^2/2)
mmgl <- ts(mmgl,
                  freq = 12,
                  start = start(ytf))
ytpronmmgl <- mmgl[,1]


#mejor modelo local

mml <- Descomp.Loess(serie.ajuste = yt,
                         h = m,
                         tipo.descomp = "multiplicative",
                         grado = 1,
                         criterio = "gcv")
ytpronmml <- mml$ytpron

#Mejor modelo errores arma ARMA(12,10), con phi7 y theta 10
param2 <- c(paste0("beta",0:6),
            "alfa1", "gamma1", "alfa2",
            "gamma2", "alfa3", "gamma3",
            "alfa4", "gamma4", "alfa5",
            "gamma5")
mma = regexpo.ErrorARMA(respuesta=yt,names.param=param2,data=X1,
                            newdata=X1nuevo,order=c(12,0,10),
                            fixed= c(NA,NA,NA,rep(0,3),NA,rep(0,4),NA,rep(0,3),NA,rep(0,4),NA,NA),
                            method="ML")
ytpronmma = mma$forecast
win.graph(width=18,height=18)
vector.auxiliar <- c(ytf,ytpronmmgl, ytpronmml, ytpronmma,predmod3[,1])
par(adj = 0.5)
plot(ytf, type = "b", pch = 19, lty = 1, col = 1, lwd = 2,
     ylab = "�ndice de ventas nominales",
     xlab = "Periodo [mmm - yy]",
     ylim = c(min(vector.auxiliar), max(vector.auxiliar)),
     xaxt = "n")
lines(ytpronmmgl, col = 2, pch = 2, lty = 2, type = "b", lwd = 2)
lines(ytpronmml, col = 3, pch = 3, lty = 3, type = "b", lwd = 2)
lines(ytpronmma, col = 4, pch = 4, lty = 4, type = "b", lwd = 2)
lines(predmod3[,1], col = 5, pch = 5, lty = 5, type = "b", lwd = 2)
legend("bottomright",
       legend = c("Real","Global", "Local", "Errores ARMA", "Modelo 3"
                   ), 
       col = 1:5,
       pch = c(19, 2:5),
       lty = 1:5,
       lwd = 1,
       cex= 0.7)
axis(1,at = time(ytf), 
     labels = c("dic-20", "ene-21", "feb-21", "mar-21",
                "abr-21", "may-21", "jun-21", "jul-21",
                "ago-21", "sep-21", "oct-21", "nov-21"))

accuracy(ytpronmma,ytf)

acf(as.numeric(residuals(mml)),ci.type="ma",lag.max=36,main="ACF modelo local",ci.col=2)

