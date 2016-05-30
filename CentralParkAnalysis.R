# Kevin O'Connor
# 4/28/16

# Analysis of Central Park temperature data.
## Reading and cleaning data.
setwd("/Users/kevinoconnor/Documents/School/Research/")
temp_data = read.table("CentralPark.txt", col.names = c("Date", "High", "Low"))
temp_data$High = temp_data$High/10
temp_data$Low = temp_data$Low/10
temp_data[which(temp_data$High == -999.9),2] = NaN
temp_data[which(temp_data$Low == -999.9),3] = NaN

## Plotting the data
pdf("MaxTemps1876to1886.pdf", width=8, height=6)
plot(temp_data$High[1:which(temp_data$Date == 18851231)], main="Daily Maximum Temperature for 1876 to 1886", ylab="Maximum Temperature (degrees Celsius)")
dev.off()
pdf("MaxTemps2005to2015.pdf", width=8, height=6)
plot(temp_data$High[which(temp_data$Date == 20050101):which(temp_data$Date == 20151230)], main="Daily Maximum Temperature for 2005 to 2015", ylab="Maximum Temperature (degrees Celsius)")
dev.off()
pdf("MinTemps1876to1886.pdf", width=8, height=6)
plot(temp_data$Low[1:which(temp_data$Date == 18851231)], main="Daily Minimum Temperature for 1876 to 1886", ylab="Minimum Temperature (degrees Celsius)")
dev.off()
pdf("MinTemps2005to2015.pdf", width=8, height=6)
plot(temp_data$Low[which(temp_data$Date == 20050101):which(temp_data$Date == 20151230)], main="Daily Minimum Temperature for 2005 to 2015", ylab="Minimum Temperature (degrees Celsius)")
dev.off()
pdf("MaxTemps1994.pdf", width=8, height=6)
plot(temp_data$High[which(temp_data$Date == 19940101):which(temp_data$Date == 19941231)], main="Daily Maximum Temperature for 1994", ylab="Maximum Temperature (degrees Celsius)")
dev.off()
pdf("MaxTemps1984.pdf", width=8, height=6)
plot(temp_data$High[which(temp_data$Date == 19840101):which(temp_data$Date == 19841231)], main="Daily Maximum Temperature for 1984", ylab="Maximum Temperature (degrees Celsius)")
dev.off()

## Carrying out triangle wave, uniform distribution simulation.
u_0 = 22.5; u_1 = 12.5
l_0 = 5; l_1 = 15
a_0 = 365.25/2
vals = runif(3650)
for(i in 1:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
}
pdf("TriangleBoundsUniformDistributionSimulation.pdf", width=13, height=8)
plot(vals, main="Triangle Bounds with Uniform Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
## Carrying out triangle wave, AR(1) process with uniform distribution simulation.
vals = runif(3650)
u_in = u_1*triangle_wave(1,a_0)+u_0
l_in = l_1*triangle_wave(1,a_0)+l_0
vals[1] = (u_in - l_in)*vals[1] + l_in
for(i in 2:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
	vals[i] = mean(vals[i:(i-1)])
}
pdf("AR1UniformSimulation.pdf", width=13, height=8)
plot(vals, main="Triangle Bounds and AR(1) Process with Uniform Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
### Repeating this example with different parametrizations.  
#### AR(2)
vals = runif(3650)
u_in = u_1*triangle_wave(1,a_0)+u_0
l_in = l_1*triangle_wave(1,a_0)+l_0
vals[1] = (u_in - l_in)*vals[1] + l_in
u_in = u_1*triangle_wave(2,a_0)+u_0
l_in = l_1*triangle_wave(2,a_0)+l_0
vals[2] = (u_in - l_in)*vals[2] + l_in
for(i in 3:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
	vals[i] = mean(vals[(i-2):i])
}
pdf("AR2UniformSimulation.pdf", width=8, height=6)
plot(vals, main="Triangle Bounds and AR(2) Process with Uniform Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
pdf("AR2UniformSimulation1Year.pdf", width=13, height=8)
plot(vals[1000:1365], main="Triangle Bounds and AR(2) Process with Uniform Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
#### AR(5)
vals = runif(3650)
for(i in 1:5){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
}
for(i in 6:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
	vals[i] = mean(vals[i:(i-5)])
}
pdf("AR5UniformSimulation.pdf", width=8, height=6)
plot(vals, main="Triangle Bounds and AR(5) Process with Uniform Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
#### AR(10)
vals = runif(3650)
for(i in 1:10){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
}
for(i in 11:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
	vals[i] = mean(vals[i:(i-10)])
}
pdf("AR10UniformSimulation.pdf", width=8, height=6)
plot(vals, main="Triangle Bounds and AR(10) Process with Uniform Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
#### AR(20)
vals = runif(3650)
for(i in 1:20){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
}
for(i in 21:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
	vals[i] = mean(vals[i:(i-20)])
}
pdf("AR20UniformSimulation.pdf", width=8, height=6)
plot(vals, main="Triangle Bounds and AR(20) Process with Uniform Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
pdf("AR20UniformSimulation1Year.pdf", width=13, height=8)
plot(vals[250:615], main="Triangle Bounds and AR(20) Process with Uniform Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
## Triangle wave, beta distribution simulation.
u_0 = 60; l_0 = -40
alpha = 10; beta = 10
hist(rbeta(1000, alpha, beta), breaks = 20)
### AR(0)
vals = rbeta(3650, alpha, beta)
for(i in 1:length(vals)){
	u_in = u_1*triangle_wave(i,a_0)+u_0
	l_in = l_1*triangle_wave(i,a_0)+l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
}
pdf("AR0BetaSimulation.pdf", width=8, height=6)
plot(vals, main="Triangle Bounds and AR(0) Process with Beta Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
### AR(1)
vals = rbeta(3650, alpha, beta)
u_in = u_1*triangle_wave(1,a_0) + u_0
l_in = l_1*triangle_wave(1,a_0) + l_0
vals[1] = (u_in - l_in)*vals[1] + l_in
for(i in 2:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
	vals[i] = mean(vals[i:(i-1)])
}
ar1_sim_data = vals
pdf("AR1BetaSimulation.pdf", width=8, height=6)
plot(vals, main="Triangle Bounds and AR(1) Process with Beta Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
### AR(2)
vals = rbeta(3650, alpha, beta)
for(i in 1:2){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
}
for(i in 3:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
	vals[i] = mean(vals[i:(i-2)])
}
ar2_sim_data = vals
pdf("AR2BetaSimulation.pdf", width=8, height=6)
plot(vals, main="Triangle Bounds and AR(2) Process with Beta Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
### AR(3)
vals = rbeta(3650, alpha, beta)
for(i in 1:3){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
}
for(i in 4:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
	vals[i] = mean(vals[i:(i-3)])
}
ar3_sim_data = vals
pdf("AR3BetaSimulation.pdf", width=8, height=6)
plot(vals, main="Triangle Bounds and AR(3) Process with Beta Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
### AR(5)
vals = rbeta(3650, alpha, beta)
for(i in 1:5){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
}
for(i in 6:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
	vals[i] = mean(vals[i:(i-5)])
}
ar5_sim_data = vals
pdf("AR5BetaSimulation.pdf", width=8, height=6)
plot(vals, main="Triangle Bounds and AR(5) Process with Beta Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
### AR(10)
vals = rbeta(3650, alpha, beta)
for(i in 1:10){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
}
for(i in 11:length(vals)){
	u_in = u_1*triangle_wave(i,a_0) + u_0
	l_in = l_1*triangle_wave(i,a_0) + l_0
	vals[i] = (u_in - l_in)*vals[i] + l_in
	vals[i] = mean(vals[i:(i-10)])
}
ar10_sim_data = vals
pdf("AR10BetaSimulation.pdf", width=8, height=6)
plot(vals, main="Triangle Bounds and AR(10) Process with Beta Distribution Simulation", ylab="Simulated Maximum Temperature (degrees Celsius)")
dev.off()
# Going to try and automate the simulation.
sim_tbbeta <- function(p, t){
	vals = rbeta(t, alpha, beta)
	for(i in 1:p){
		u_in = u_1*triangle_wave(i,a_0)+u_0
		l_in = l_1*triangle_wave(i,a_0)+l_0
		vals[i] = (u_in-l_in)*vals[i]+l_in
	}
	for(i in (p+1):length(vals)){
		u_in = u_1*triangle_wave(i,a_0)+u_0
		l_in = l_1*triangle_wave(i,a_0)+l_0
		vals[i] = (u_in-l_in)*vals[i]+l_in
		vals[i] = mean(vals[i:(i-p)])
	}
	return(vals)
}

# Comparing dependence to that of the data.
pdf("ACFSimData.pdf", width=8, height=6)
acf(vals, lag=1000, main="ACF Simulated Data")
dev.off()
pdf("ACFCentralParkData.pdf", width=8, height=6)
acf(temp_data$High, lag=1000, na.action=na.pass, main="ACF Central Park Data")
dev.off()
pdf("PACFSimData.pdf", width=8, height=6)
pacf(vals, lag=100, main="PACF Simulated Data")
dev.off()
pdf("PACFCentralParkData.pdf", width=8, height=6)
pacf(temp_data$High, lag=100, na.action=na.pass, main="PACF Central Park Data")
dev.off()

## 5/26/16 
## Deseasonalizing and comparing deseasonalized data to deseasonalized simulation data
library(stats)
library(zoo)
temp_data$High = na.approx(temp_data$High)
temp_data_hts = ts(temp_data$High, start=c(1876, 1), deltat=1/365.25)
dec_hts = stl(as.ts(temp_data_hts), "periodic", na.action=na.pass, s.window="period")$time.series
deseas_hts = ts(dec_hts[,2]+dec_hts[,3], start=c(1876, 1), deltat=1/365.25)
plot(deseas_hts, main="Deseasonalized Central Park Max Temperatures", ylab="Deseasonalized Max Temperature")
acf(deseas_hts, lag=30)
pacf(deseas_hts, lag=30)
### Deseasonalizing simulated data from triangle-bounded Beta model
#### ar1
ar1_sim_data = ts(ar1_sim_data, start=c(1876, 1), deltat=1/365)
dec_ar1 = stl(as.ts(ar1_sim_data), "periodic", s.window="period")$time.series
deseas_ar1 = ts(dec_ar1[,2]+dec_ar1[,3], deltat=1/365, start=c(1876, 1))
plot(deseas_ar1)
acf(deseas_ar1, lag=30)
pacf(deseas_ar1, lag=30)
#### ar5
ar5_sim_data = ts(ar5_sim_data, start=c(1876,1), deltat=1/365)
dec_ar5 = stl(as.ts(ar5_sim_data), "periodic", s.window="period")$time.series
deseas_ar5 = ts(dec_ar5[,2]+dec_ar5[,3], deltat=1/365, start=c(1876,1))
plot(deseas_ar5)
acf(deseas_ar5, lag=30)
pacf(deseas_ar5, lag=30)
#### ar15
alpha=3; beta=3
ar15_sim_data = ts(sim_tbbeta(15,10*365), start=c(1876,1), deltat=1/365)
dec_ar15 = stl(as.ts(ar15_sim_data), "periodic", s.window="period")$time.series
deseas_ar15 = ts(dec_ar15[,2]+dec_ar15[,3], deltat=1/365, start=c(1876,1))
plot(ar15_sim_data, ylim=c(-50, 100))
points(time(ar15_sim_data), triangle_wave(1:3650, a_0)*u_1+u_0)
points(time(ar15_sim_data), triangle_wave(1:3650, a_0)*l_1+l_0)
plot(deseas_ar15)
acf(deseas_ar15, lag=30)
pacf(deseas_ar15, lag=30)

## 5/29/16
### Independent Beta simulation with triangle bounds.
sim.length = 3650
u_0 = 30
vals = rbeta(sim.length, 10, 10)
for(i in 1:sim.length){
	u = u_1*triangle_wave(i,a_0)+u_0
	l = l_1*triangle_wave(i,a_0)+l_0
	vals[i] = (u-l)*vals[i] + l
}
pdf("IndependentBetaSimWithTriangleBounds.pdf", width=13, height=8)
plot(vals, main="Simulation of Independent Beta with Triangle Bounds", ylab="Simulated Temperature")
dev.off()
### Fitting bounds to the data.
plot(temp_data$High)
u_0 = 45; u_1 = 0
l_0 = -15; l_1 = 0
delta = 0.1
is_valid_bound <- function(u_0, u_1, l_0, l_1){
	valid = T
	for (i in 1:length(temp_data$High)){
		if (!is.na(temp_data$High[i])){
			if ((temp_data$High[i] > u_1*triangle_wave(i,a_0)+u_0) || (temp_data$High[i] < l_1*triangle_wave(i,a_0)+l_0) ){
				valid = F
				break
			}
		}
	}
	return(valid)
}
new_u_0 = 100
new_u_1 = 0
while(new_u_0 != u_0 || new_u_1 != u_1){
	if(is_valid_bound(u_0-delta, u_1+delta, l_0, l_1)){
		new_u_0 = u_0 - delta
		new_u_1 = u_1 + delta	
	} else if(is_valid_bound(u_0-delta, u_1, l_0, l_1)) {
		new_u_0 = u_0 - delta
	} else if(is_valid_bound(u_0, u_1+delta, l_0, l_1)) {
		new_u_1 = u_1 + delta
	} 
	u_0 = new_u_0
	u_1 = new_u_1
}

## Simulations
### Setting parameters
a_0 = 365.25/2 # half period



## Kernel Functions
### Triangle kernel
tri_kernel <- function(i, p){
	return((2*p - 2*i + 	1)/((p+1)^2))
}
### Epanechnikov kernel
epa_kernel <- function(i, p){
	return((3*p^2+6*p-3*i^2-3*i+2)/(2*(p+1)^3))
}

## Bound function
### Triangle wave
triangle_wave <- function(t,a){
	return((2/a)*(t-a*floor(t/a+1/2))*(-1)^(floor(t/a+1/2)))
}

## TESTING ##
### Testing kernels
sum(tri_kernel(seq(0,15), 15)) == 1
sum(epa_kernel(seq(0,15), 15)) == 1
tri_kernel(seq(0,15), 15) > 0
epa_kernel(seq(0,15), 15) > 0
