# Kevin O'Connor
# Modeling temperature as a beta distributed random variable

## Independent, constant bounds.
# Generating data
maxvec1 = c()
l = 45; u = 105
for (i in 1:10000){
	data = l + (u-l)*rbeta(10000, 100, 100)
	maxvec1 = c(maxvec1, max(data))
}
# Histogram of data
pdf("/Users/kevinoconnor/Documents/School/Research/BetaIndConstantBounds1.pdf", width=13, height=8)
hist(maxvec1, main="Max of 10000 iid beta(100,100) with constant bounds [45, 105]")
dev.off()
# Fitting to GEV distribution
gevcdf <- function(x, mu, sigma, xi){
	return(exp(-(1 + xi*((x-mu)/sigma))^(-1/xi)))
}
fit1 = gev.fit(maxvec1)
# Plotting empirical CDF with fitted CDF
pdf("/Users/kevinoconnor/Documents/School/Research/BetaPlot7.pdf", width=13, height=8)
plot(ecdf(maxvec1), main="Empirical CDF with Fit Points")
points(maxvec1, gevcdf(maxvec1, fit1$mle[1], fit1$mle[2], fit1$mle[3]))
dev.off()


## Independent, sinusoidal bounds.
l <- function(i,n){
	return(l_0+l_1*sin(2*pi*i/n))
}
u <- function(i,n){
	return(u_0+u_1*sin(2*pi*i/n))
}
# u_1 = 20, l_1 = 20
maxvec2 = c(); maxindvec2 = c()
l_0 = 10; l_1 = 20; u_0 = 115; u_1 = 20
n = 365
for (j in 1:1000){
	data = c()
	for (i in 1:365){
		data = c(data, l(i,n) + (u(i,n)-l(i,n))*rbeta(1, 100, 100))
	}	
		maxvec2 = c(maxvec2, max(data))
		maxindvec2 = c(maxindvec2, which.max(data))
}
pdf("/Users/kevinoconnor/Documents/School/Research/BetaPlot3.pdf", width=13, height=8)
hist(maxvec2, main="Simulated Annual Maxima with Sinusoidal Bounds")
dev.off()
pdf("/Users/kevinoconnor/Documents/School/Research/BetaPlot4.pdf", width=13, height=8)
hist(maxindvec2, breaks=20, main="Day at Occurrence of Annual Maximum", xlim=c(1,365))
dev.off()
# u_1 = 40, l_1 = 40
maxvec3 = c(); maxindvec3 = c()
l_0 = 10; l_1 = 40; u_0 = 115; u_1 = 40
n = 365
for (j in 1:1000){
	data = c()
	for (i in 1:365){
		data = c(data, l(i,n) + (u(i,n)-l(i,n))*rbeta(1, 100, 100))
	}	
		maxvec3 = c(maxvec3, max(data))
		maxindvec3 = c(maxindvec3, which.max(data))
}
pdf("/Users/kevinoconnor/Documents/School/Research/BetaPlot5.pdf", width=13, height=8)
hist(maxvec3, main="Simulated Annual Maxima with Sinusoidal Bounds")
dev.off()
pdf("/Users/kevinoconnor/Documents/School/Research/BetaPlot6.pdf", width=13, height=8)
hist(maxindvec3, breaks=20, main="Day at Occurrence of Annual Maximum", xlim=c(1,365))
dev.off()
# u_1 = 5, l_1 = 5
maxvec4 = c(); maxindvec4 = c()
l_0 = 10; l_1 =5; u_0 = 115; u_1 = 5
n = 365
for (j in 1:1000){
	data = c()
	for (i in 1:n){
		data = c(data, l(i,n) + (u(i,n)-l(i,n))*rbeta(1, 100, 100))
	}	
		maxvec4 = c(maxvec4, max(data))
		maxindvec4 = c(maxindvec4, which.max(data))
}
pdf("/Users/kevinoconnor/Documents/School/Research/BetaPlot8.pdf", width=13, height=8)
hist(maxvec4, main="Simulated Annual Maxima with Sinusoidal Bounds")
dev.off()
pdf("/Users/kevinoconnor/Documents/School/Research/BetaPlot9.pdf", width=13, height=8)
hist(maxindvec4, breaks=20, main="Day at Occurrence of Annual Maximum", xlim=c(1,365))
dev.off()
# Sinusoidal bounds, u_0-l_0 large.
maxvec5 = c(); maxindvec5 = c()
l_0 = -10^6; l_1 =20; u_0 = 115; u_1 = 20
n = 365
l <- function(i,n){
	return(l_0+l_1*sin(2*pi*i/n))
}
u <- function(i,n){
	return(u_0+u_1*sin(2*pi*i/n))
}
for (j in 1:1000){
	data = c()
	for (i in 1:n){
		data = c(data, l(i,n) + (u(i,n)-l(i,n))*rbeta(1, 100, 100))
	}	
		maxvec5 = c(maxvec5, max(data))
		maxindvec5 = c(maxindvec5, which.max(data))
}
pdf("/Users/kevinoconnor/Documents/School/Research/BetaPlot10.pdf", width=13, height=8)
hist(maxvec5, main="Simulated Annual Maxima with Sinusoidal Bounds")
dev.off()
pdf("/Users/kevinoconnor/Documents/School/Research/BetaPlot11.pdf", width=13, height=8)
hist(maxindvec5, breaks=20, main="Day at Occurrence of Annual Maximum", xlim=c(1,365))
dev.off()
