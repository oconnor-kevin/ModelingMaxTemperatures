# Kevin O'Connor
# 4/28/16

### TABLE OF CONTENTS ###
# README
# Reading and cleaning data
# Plots of Central Park data
# Deseasonalization
# Parameter estimation
# Simulations
# Kernel (weighting) functions
# Bound functions
# Simulation routines
# Bound fitting routines
# Standardization function
# Dependence removal function
# Alpha/Beta estimation functions
# TESTING

# README
# This file was created by Kevin O'Connor at the University of Chicago as part of a research project on Extreme Temperature Models under Michael Stein.  The paper for which this code was written can be found on my GitHub (https://github.com/oconnor-kevin).  Please email me at oconnor.kevin.ant@gmail.com with any questions or comments.


# Analysis of Central Park temperature data.
## Reading and cleaning data.
setwd("/Users/kevinoconnor/Documents/School/Research/")
temp_data = read.table("CentralPark.txt", col.names = c("Date", "High", "Low"))
temp_data$High = temp_data$High/10
temp_data$Low = temp_data$Low/10
temp_data[which(temp_data$High == -999.9),2] = NaN
temp_data[which(temp_data$Low == -999.9),3] = NaN


## Plots
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
pdf("DeseasonalizedMax.pdf", width=13, height=8)
plot(deseas_hts, main="Deseasonalized Central Park Max Temperatures", ylab="Deseasonalized Max Temperature")
dev.off()
pdf("DeseasonalizedMaxACF.pdf", width=8, height=6)
acf(deseas_hts, lag=30)
dev.off()
pdf("DeseasonalizedMaxPACF.pdf", width=8, height=6)
pacf(deseas_hts, lag=30)
dev.off()
pdf("EstimatedBoundsMax.pdf", width=13, height=8)
plot(temp_data$High[1:1000], main="Daily Maximum Temperature with Fitted Triangle Wave Bounds", ylim=c(-25,60), ylab="Daily Maximum Temperature")
points(upper_params[[1]][2]*triangle_wave(seq(1,1000)-upper_params[[1]][3],a_0)+upper_params[[1]][1])
points(lower_params[[1]][2]*triangle_wave(seq(1,1000)-lower_params[[1]][3],a_0)+lower_params[[1]][1])
dev.off()
pdf("MaxTempsWithoutEstimatedBounds.pdf", width=13, height=8)
plot(unseason_data, main="Daily Maximum Temperature Deseasonalized Based on Estimated Bounds", ylab="Standardized Daily Maximum Temperature")
dev.off()
pdf("MaxTempsWithoutEstimatedBoundsPACF.pdf", width=13, height=8)
pacf(unseason_data)
dev.off()
pdf("MaxTempsAlphaBetaEstimationUniform.pdf", width=8, height=6)
plot(unif_alphas, pch=20, main="Alpha and Beta vs. Autoregressive Order for Uniform Weighting", xlab="Order", ylab="")
points(unif_betas, pch=1)
legend("topright", legend=c("alpha", "beta"), pch=c(20,1))
dev.off()
pdf("MaxTempsAlphaBetaEstimationTri.pdf", width=8, height=6)
plot(tri_alphas, pch=20, main="Alpha and Beta vs. Autoregressive Order for Triangle Weighting", xlab="Order", ylab="")
points(tri_betas, pch=1)
legend("topright", legend=c("alpha", "beta"), pch=c(20,1))
dev.off()
pdf("MaxTempsAlphaBetaEstimationEpa.pdf", width=8, height=6)
plot(epa_alphas, pch=20, main="Alpha and Beta vs. Autoregressive Order for Epanechnikov Weighting", xlab="Order", ylab="")
points(epa_betas, pch=1)
legend("topright", legend=c("alpha", "beta"), pch=c(20,1))
dev.off()
pdf("UnifSimAR1.pdf", width=8, height=6)
plot(unif_sim_data[,1], main="Uniform Weighting AR(1) Simulation", ylab="Simulated Temperature", ylim=c(-20,50))
dev.off()
pdf("UnifSimAR2.pdf", width=8, height=6)
plot(unif_sim_data[,2], main="Uniform Weighting AR(2) Simulation", ylab="Simulated Temperature", ylim=c(-20,50))
dev.off()
pdf("TriSimAR1.pdf", width=8, height=6)
plot(tri_sim_data[,1], main="Triangle Weighting AR(1) Simulation", ylab="Simulated Temperature", ylim=c(-20,50))
dev.off()
pdf("TriSimAR2.pdf", width=8, height=6)
plot(tri_sim_data[,2], main="Triangle Weighting AR(2) Simulation", ylab="Simulated Temperature", ylim=c(-20,50))
dev.off()
pdf("EpaSimAR1.pdf", width=8, height=6)
plot(epa_sim_data[,1], main="Epanechnikov Weighting AR(1) Simulation", ylab="Simulated Temperature", ylim=c(-20,50))
dev.off()
pdf("EpaSimAR2.pdf", width=8, height=6)
plot(epa_sim_data[,2], main="Epanechnikov Weighting AR(2) Simulation", ylab="Simulated Temperature", ylim=c(-20,50))
dev.off()
pdf("TriSimAR2Max.pdf", width=8, height=6)
hist(tri_maxes[,2], main="Simulated Annual Maxima for Triangle Weighted AR(2)", xlab="Simulated Annual Max Temperature")
dev.off()
pdf("TriSimAR2MaxInds.pdf", width=8, height=6)
hist(tri_max_inds[,2], main="Indices of Simulated Annual Maxima for Triangle Weighted AR(2)", xlab="Indices of Simulated Annual Max Temperature")
dev.off()
pdf("EpaSimAR2Max.pdf", width=8, height=6)
hist(epa_maxes[,2], main="Simulated Annual Maxima for Epanechnikov Weighted AR(2)", xlab="Simulated Annual Max Temperature")
dev.off()
pdf("EpaSimAR2MaxInds.pdf", width=8, height=6)
hist(epa_max_inds[,2], main="Indices of Simulated Annual Maxima for Epanechnikov Weighted AR(2)", xlab="Indices of Simulated Annual Max Temperature")
dev.off()
pdf("TriSimAR2MaxGEVDiag.pdf", width=13, height=8)
gev.diag(fit_tri_max)
dev.off()
pdf("EpaSimAR2MaxGEVDiag.pdf", width=13, height=8)
gev.diag(fit_epa_max)
dev.off()


## Deseasonalizing and comparing deseasonalized data to deseasonalized simulation data
library(stats)
library(zoo)
temp_data$High = na.approx(temp_data$High)
temp_data_hts = ts(temp_data$High, start=c(1876, 1), deltat=1/365.25)
dec_hts = stl(as.ts(temp_data_hts), "periodic", na.action=na.pass, s.window="period")$time.series
deseas_hts = ts(dec_hts[,2]+dec_hts[,3], start=c(1876, 1), deltat=1/365.25)


## Parameter Estimation
### Fitting bounds to the data.
#### To save on computation time, we will fit the model using a random subset of the data
inds = sample(seq(1,length(temp_data$High)), 1000) # Generate random indices
fit_data = temp_data$High[inds] # Store selected data
delta = 0.01; a_0 = 365.25/2; delta_s = 0.5; shift = -241.5 # Initialize for bound optimization
u_0 = 40; u_1 = 0; l_0 = -15; l_1 = 0
upper_params = optimize_upper_bound(u_0,u_1,temp_data$High[1:1000]) # Generate upper bounds
lower_params = optimize_lower_bound(l_0,l_1,temp_data$High[1:1000]) # Generate lower bounds
### Deseasonalizing based on the estimated bounds.
unseason_data = standardize_data(upper_params[[1]][1], upper_params[[1]][2], lower_params[[1]][1], lower_params[[1]][2], upper_params[[1]][3], temp_data$High)
### Estimating alpha and beta.  
indep_unif_data = c(); unif_alphas = c(); unif_betas = c()
indep_tri_data = c(); tri_alphas = c(); tri_betas = c()
indep_epa_data = c(); epa_alphas = c(); epa_betas = c()
for(i in 1:10){
	indep_unif_data = cbind(indep_unif_data, na.approx(remove_dependence(i, rep(1,i+1)/(i+1), unseason_data)))
	indep_tri_data = cbind(indep_tri_data, na.approx(remove_dependence(i, tri_kernel(seq(0,i),i), unseason_data)))
	indep_epa_data = cbind(indep_epa_data, na.approx(remove_dependence(i, epa_kernel(seq(0,i),i), unseason_data)))
	unif_alphas = c(unif_alphas, estimate_alpha(mean(indep_unif_data[,i]), var(indep_unif_data[,i])))
	tri_alphas = c(tri_alphas, estimate_alpha(mean(indep_tri_data[,i]), var(indep_tri_data[,i])))
	epa_alphas = c(epa_alphas, estimate_alpha(mean(indep_epa_data[,i]), var(indep_epa_data[,i])))
	unif_betas = c(unif_betas, estimate_beta(mean(indep_unif_data[,i]), var(indep_unif_data[,i])))
	tri_betas = c(tri_betas, estimate_beta(mean(indep_tri_data[,i]), var(indep_tri_data[,i])))
	epa_betas = c(epa_betas, estimate_beta(mean(indep_epa_data[,i]), var(indep_epa_data[,i])))
}


## Simulations
### Setting parameters
a_0 = 365.25/2 # half period
p = 2         # order
alpha = 5      # beta shape parameter 1
beta = 5       # beta shape parameter 2
l_0 = lower_params[[1]][1]      # parameters for bounds
l_1 = lower_params[[1]][2]
u_0 = upper_params[[1]][1]
u_1 = upper_params[[1]][2]
years = 100
### Generating data for all orders up to p=10.  
unif_sim_data = c()
tri_sim_data = c()
epa_sim_data = c()
for (i in 1:2){
	p = i # Setting order
	alpha = unif_alphas[i]; beta = unif_betas[i]
	unif_sim_data = cbind(unif_sim_data, sim_unif_beta(365*years)) # Generating data
	alpha = tri_alphas[i]; beta = tri_betas[i]
	tri_sim_data = cbind(tri_sim_data, sim_tri_beta(365*years))
	alpha = epa_alphas[i]; beta = epa_betas[i]
	epa_sim_data = cbind(epa_sim_data, sim_epa_beta(365*years))
}
### Finding annual maxima.
unif_maxes = c(); unif_max_inds = c()
tri_maxes = c(); tri_max_inds = c()
epa_maxes = c(); epa_max_inds = c()
for (i in 1:2){
	unif_maxes_0 = c(); unif_max_ind = c()
	tri_maxes_0 = c(); tri_max_ind = c()
	epa_maxes_0 = c(); epa_max_ind = c()
	for (j in 1:years){
		unif_max_ind = c(unif_max_ind, which(unif_sim_data[((j-1)*365+1):(j*365),i]==max(unif_sim_data[((j-1)*365+1):(j*365),i])))
		tri_max_ind = c(tri_max_ind, which(tri_sim_data[((j-1)*365+1):(j*365),i]==max(tri_sim_data[((j-1)*365+1):(j*365),i])))
		epa_max_ind = c(epa_max_ind, which(epa_sim_data[((j-1)*365+1):(j*365),i]==max(epa_sim_data[((j-1)*365+1):(j*365),i])))
		unif_maxes_0 = c(unif_maxes_0, unif_sim_data[unif_max_ind[j]])
		tri_maxes_0 = c(tri_maxes_0, tri_sim_data[tri_max_ind[j]])
		epa_maxes_0 = c(epa_maxes_0, epa_sim_data[epa_max_ind[j]])
		
	}
	unif_maxes = cbind(unif_maxes, unif_maxes_0)
	tri_maxes = cbind(tri_maxes, tri_maxes_0)
	epa_maxes = cbind(epa_maxes, epa_maxes_0)
	unif_max_inds = cbind(unif_max_inds, unif_max_ind)
	tri_max_inds = cbind(tri_max_inds, tri_max_ind)
	epa_max_inds = cbind(epa_max_inds, epa_max_ind)
}
### Fitting annual maxima to a GEV distribution.
library(ismev)
fit_tri_max = gev.fit(tri_maxes[,2])
fit_epa_max = gev.fit(epa_maxes[,2])
gev.diag(fit_tri_max)
gev.diag(fit_epa_max)


## Kernel Functions
### Triangle kernel
tri_kernel <- function(i, p){
	return((2*p - 2*i + 	1)/((p+1)^2))
}
### Epanechnikov kernel
epa_kernel <- function(i, p){
	return((3*p^2+6*p-3*i^2-3*i+2)/(2*(p+1)^3))
}


## Bound functions
### Triangle wave
triangle_wave <- function(t,a){
	return((2/a)*(t-a*floor(t/a+1/2))*(-1)^(floor(t/a+1/2)))
}
### Triangle wave with shift
shift_triangle_wave <- function(t,a,s){
	return(triangle_wave((t-s), a))
}


## Simulation routines
### AR Beta Simulation with Arbitrary Weighting Vector
#### w is a weight vector of length p+1 which should sum to 1 and be in order from i=0 to p
#### Initialize p, alpha, beta, u_1, l_1, u_0, l_0, a_0
sim_beta <- function(t, w){
	pre_vals = rbeta(p, alpha, beta) # Beta variables used just for initialization
	vals = rbeta(t, alpha, beta) # Beta variables
	ss_vals = c() # Shifted beta variables
	for(i in (p+1):(length(vals)+p)){
		vals[i-p] = sum(w*(c(pre_vals, vals)[i:(i-p)])) # Applying weight vector
		u_in = u_1*triangle_wave(i,a_0)+u_0 # Creating bounds
		l_in = l_1*triangle_wave(i,a_0)+l_0
		ss_vals = c(ss_vals, (u_in - l_in)*vals[i-p]+l_in) # Shifting and scaling beta variable
	}
	return(ss_vals)
}
### Uniform Weighting AR Beta Simulation
#### Initialize p, alpha, beta, u_1, l_1, u_0, l_0, a_0
sim_unif_beta <- function(t){
	return(sim_beta(t, rep(1,p+1)/(p+1)))
}
### Triangle Weighting AR Beta Simulation
sim_tri_beta <- function(t){
	return(sim_beta(t, tri_kernel(seq(0,p),p)))	
}
### Epanechnikov Weighting AR Beta Simulation
sim_epa_beta <- function(t){
	return(sim_beta(t, epa_kernel(seq(0,p),p)))
}


## Bound fitting routines
### Checks that bounds are valid
#### Initialize a_0, shift
is_valid_bound <- function(u_0, u_1, l_0, l_1, shift, x){
	valid = T
	for (i in 1:length(x)){
		if (!is.na(x[i])){
			if ((x[i] > u_1*triangle_wave(i-shift,a_0)+u_0) || (x[i] < l_1*triangle_wave(i-shift,a_0)+l_0) ){
				valid = F
				break
			}
		}
	}
	return(valid)
}
### Optimize upper bound
#### Initialize delta, a_0, delta_s
optimize_upper_bound <- function(h_0, h_1, x){
	done = F
	while(!done){
		if(is_valid_bound(h_0-delta, h_1+delta, -Inf, 0, shift, x)){
			h_0 = h_0 - delta
			h_1 = h_1 + delta	
		} else if (is_valid_bound(h_0-delta, h_1, -Inf, 0, shift, x)){
			h_0 = h_0 - delta
		} else if (is_valid_bound(h_0, h_1+delta, -Inf, 0, shift, x)){
			h_1 = h_1 + delta
		} else if (is_valid_bound(h_0, h_1, -Inf, 0, shift-delta_s, x)){
			shift = shift-delta_s
		} else {
			done = T
		}
	}
	return(list(c(h_0, h_1, shift)))
}
### Optimize lower bound
#### Initialize delta, a_0, delta_s
optimize_lower_bound <- function(h_0, h_1, x){
	done = F
	while(!done){
		if(is_valid_bound(Inf, 0, h_0+delta, h_1+delta, shift, x)){
			h_0 = h_0 + delta
			h_1 = h_1 + delta	
		} else if (is_valid_bound(Inf, 0, h_0+delta, h_1, shift, x)){
			h_0 = h_0 + delta
		} else if (is_valid_bound(Inf, 0, h_0, h_1+delta, shift, x)){
			h_1 = h_1 + delta
		} else if (is_valid_bound(Inf, 0, h_0, h_1, shift-delta_s, x)){
			shift = shift-delta_s
		} else {
			done = T
		}
	}
	return(list(c(h_0, h_1, shift)))
}


## Standardization function
### Initialize a_0
standardize_data <- function(u_0_0, u_1_0, l_0_0, l_1_0, sh, x){
	std_vals = x
	for (i in 1:length(x)){
		u = u_1_0*shift_triangle_wave(i, a_0, sh)+u_0_0
		l = l_1_0*shift_triangle_wave(i, a_0, sh)+l_0_0
		std_vals[i] = (x[i] - l)/(u-l)
	}
	return(std_vals)
}


## Dependence Removal Function
remove_dependence <- function(order, weights, x){
	new_x = x
	for (i in (order+1):length(x)){
		new_x[i] = (x[i] - sum(weights[-1]*x[(i-1):(i-order)]))/weights[1]
	}
	return(new_x)
}

## Alpha/Beta Estimation Functions
estimate_alpha <- function(mu, sigma_sq){
	return(mu*((mu*(1-mu))/(sigma_sq)-1))
}
estimate_beta <- function(mu, sigma_sq){
	return((1-mu)*((mu*(1-mu))/(sigma_sq)-1))	
}


## TESTING ##

### Testing kernels
sum(tri_kernel(seq(0,15), 15)) == 1
sum(epa_kernel(seq(0,15), 15)) == 1
tri_kernel(seq(0,15), 15) > 0
epa_kernel(seq(0,15), 15) > 0

### Testing simulation routines
p=5; alpha=3; beta=3; a_0 = 365.25
u_1 = 10; u_0 = 70; l_1 = 10; l_0 = -10
plot(sim_unif_beta(3650))
plot(sim_tri_beta(3650))
plot(sim_epa_beta(3650))

### Testing bound optimization
delta = 0.01; a_0 = 365.25/2; delta_s = 0.5
shift = a_0/2 - (which(triangle_wave(seq(1,1000),365.25/2)==min(triangle_wave(seq(1,1000), 365.25/2))) - which(temp_data$High[1:1000]==min(temp_data$High[1:1000])))
u_0 = 40; u_1 = 0; l_0 = -15; l_1 = 0
upper_params = optimize_upper_bound(u_0,u_1,temp_data$High[1:1000])
lower_params = optimize_lower_bound(l_0,l_1,temp_data$High[1:1000])
plot(temp_data$High[1:1000], ylim=c(-20,60))
points(upper_params[[1]][2]*triangle_wave(seq(1,1000)-upper_params[[1]][3],a_0)+upper_params[[1]][1])
points(lower_params[[1]][2]*triangle_wave(seq(1,1000)-lower_params[[1]][3],a_0)+lower_params[[1]][1])
