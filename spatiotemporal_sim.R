##  spatiotemporal_sim.R: script to simulate data reflecting two species
##  percent cover across a spatially and temporally varying covariate.



rm(list=ls(all.names = TRUE))

####
####  LIBRARIES ----
####
library(mvtnorm)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggthemes)
library(rjags)
library(coda)
library(synchrony)



####
####  MY PLOTTING THEME ----
####
my_theme <- theme_bw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color="white"),
        panel.background   = element_rect(fill = "#EFEFEF"),
        axis.text          = element_text(size=10, color="grey35", family = "Arial Narrow"),
        axis.title         = element_text(size=12, family = "Arial Narrow", face = "bold"),
        panel.border       = element_blank(),
        axis.line.x        = element_line(color="black"),
        axis.line.y        = element_line(color="black"),
        strip.background   = element_blank(),
        strip.text         = element_text(size=10, color="grey35", family = "Arial Narrow"),
        legend.title       = element_text(size=10, family = "Arial Narrow"),
        legend.text        = element_text(size=8, color="grey35", family = "Arial Narrow"))




####
#### SIMULATE ENVRIONMENTAL DATA ----
####
yrs <- 10
lat <- c(-11:10)
lon <- c(-11:10)
saddle_length <- 0.5*length(lat)
snowdepth <- c(seq(-5,5, length.out = saddle_length),
               seq(5,-5, length.out = saddle_length))
environment <- rbeta(n = yrs, shape1 = 1, shape2 = 1)

snowgrid_time <- list()
for(t in 1:yrs){
  snowgrid <- matrix(data = NA, ncol = length(lon), nrow = length(lat))
  for(i in 1:length(lat)){
    snowgrid[i,] <- rnorm(ncol(snowgrid),snowdepth[i],sd=1) 
    snowgrid[i,] <- snowgrid[i,] * environment[t]
  }
  snowgrid_time[[t]] <- snowgrid
}

### Grab first year grid for plotting example
snowgrid_df <- as.data.frame(snowgrid_time[[1]])
colnames(snowgrid_df) <- lon
snowgrid_df$lat <- lat
snowgrid_long <- melt(snowgrid_df, id.vars = "lat")
colnames(snowgrid_long) <- c("lat", "lon", "snowdepth")
snowgrid_long$lon <- as.numeric(as.character(snowgrid_long$lon))

### Set up evenly spaced sample grid
lon_ids <- seq.int(1L, length(lon), 3L)
lat_ids <- seq.int(1L, length(lat), 3L)
sample_lons <- lon[lon_ids]
sample_lats <- lat[lat_ids]
sample_grid <- expand.grid(sample_lats, sample_lons)
colnames(sample_grid) <- c("sample_lat", "sample_lon")

### Extract snow depth values at each sample location for each year
snow_sample <- list()
for(t in 1:yrs){
  snow_sample[[t]] <- snowgrid_time[[t]][lon_ids, lat_ids]
}


### Make example plot for snow depth gradient and sample points
ggplot(snowgrid_long, aes(x=lon, y=lat))+
  geom_raster(data = snowgrid_long, aes(x=lon, y=lat, fill=snowdepth))+
  geom_point(data=sample_grid, aes(x=sample_lon, y=sample_lat), color="coral")+
  geom_point(data=sample_grid, aes(x=sample_lon, y=sample_lat), shape=1)+
  scale_fill_viridis(direction = -1, "Snow depth")+
  xlab("Easting")+
  ylab("Northing")+
  my_theme
ggsave(filename = "snowdepth_grid_example.png", width = 4, height=2.5, units = "in", dpi=300)



####
####  GENERATE SPECIES COVER DATA ----
####
spp1_var <- 0.1
spp2_var <- 0.2
beta_int <- 0
beta_dd1 <- 0.8
beta_dd2 <- 0.8
beta_mu <- -1.5
beta_spp_mu <- rnorm(2, beta_mu, sd = 0.5)
sigma <- 1
rho <- 0
Sigma <- matrix(data = c(sigma, sigma*rho,
                         sigma*rho, sigma), 2,2)
Sigma <- as.matrix(Sigma)
beta_spp_yrs <- rmvnorm(yrs, mean = beta_spp_mu, sigma = Sigma)

beta_spp_int <- rnorm(2, beta_int, sd=2)
sigma_int <- 2
rho_int <- 0.5
Sigma_int <- matrix(data = c(sigma_int, sigma_int*rho_int,
                         sigma_int*rho_int, sigma_int), 2,2)
Sigma_int <- as.matrix(Sigma_int)
beta_int_yrs <- rmvnorm(yrs, mean = beta_spp_int, sigma = Sigma_int)

mu_spp1 <- mu_spp2 <- list()
mu_spp1[[1]] <- rnorm(64, snow_sample[[1]],0.01)
mu_spp2[[1]] <- rnorm(64, snow_sample[[1]],0.01)
for(t in 2:yrs){
  mu_spp1[[t]] <- rnorm(64,beta_int_yrs[t,1] + beta_dd1*mu_spp1[[t-1]] + beta_spp_mu[1]*snow_sample[[t]], 2)
  mu_spp2[[t]] <- rnorm(64,beta_int_yrs[t,2] + beta_dd2*mu_spp2[[t-1]] + beta_spp_mu[2]*snow_sample[[t]], 2) 
}

### Convert to really long dataframe
sim_data <- list()
for(t in 1:yrs){
  tmp1 <- as.numeric(mu_spp1[[t]])
  tmp2 <- as.numeric(mu_spp2[[t]])
  tmp_data <- data.frame(year = t,
                         lat = sample_grid$sample_lat,
                         lon = sample_grid$sample_lon,
                         plot_id = c(1:length(tmp1)),
                         spp1_value = rnorm(length(tmp1),tmp1,spp1_var),
                         spp2_value = rnorm(length(tmp2),tmp2,spp2_var),
                         snow_depth = as.numeric(snow_sample[[t]]))
  sim_data <- rbind(sim_data, tmp_data)
}

### Melt to make even longer dataframe
sim_data_long <- melt(sim_data, id.vars = c("year","lat","lon","plot_id","snow_depth"))
sim_data_long2 <- sim_data_long
sim_data_long2$year <- sim_data_long2$year+1
colnames(sim_data_long2)[7] <- "lag_value"
sim_data_long3 <- merge(sim_data_long, sim_data_long2[,c("year", "plot_id", "variable", "lag_value")])
sim_data_long <- sim_data_long3
sim_data_long$year <- sim_data_long$year-1


####
####  JAGS STATE-SPACE MODEL ----
####
my_model <- "  
  model{

    #### Variance Priors
    tau_obs ~ dgamma(0.0001, 0.0001)
    tau_proc ~ dgamma(0.0001, 0.0001)
    sigma_proc <- 1/sqrt(tau_proc)
    for(i in 1:2) { 
      betayr_tau[i] ~ dgamma(0.0001,0.0001) 
      intyr_tau[i] ~ dgamma(0.0001,0.0001) 
    }

    #### Fixed Effects Priors
    beta0 ~ dnorm(0,0.00001)
    beta1 ~ dnorm(0,0.00001)
    beta2_mu ~ dnorm(0,0.00001)
    for(i in 1:nspp){
      beta2[i] ~ dnorm(0, 0.00001)
      int[i] ~ dnorm(0,0.00001)
      for(t in 1:nyrs){
        int_yr[t,i] ~ dnorm(int[i], intyr_tau[i])
      }
    }

    #### Process Model
    for(t in 1:n){
      Nmed[t] <- int_yr[yr[t],spp[t]] + beta1*lag_value[t] + beta2[spp[t]]*snow[t]
    }

  #### Data Model
  for(i in 1:n){
    Nobs[i] ~ dnorm(Nmed[i], tau_obs)
  }

}"



####
####  FIT STATE-SPACE MODEL ----
####
###  Prepare data list
snow_by_year <- ddply(sim_data_long, .(year), summarise,
                      avg_snowdepth = mean(snow_depth))
mydat         <- list(Nobs = sim_data_long$value, 
                      lag_value = sim_data_long$lag_value,
                      snow = sim_data_long$snow_depth,
                      n = nrow(sim_data_long),
                      yr = sim_data_long$year,
                      nyrs = length(unique(sim_data_long$year)),
                      npreds = length(unique(sim_data_long$year)),
                      spp = as.numeric(sim_data_long$variable),
                      nspp = length(unique(sim_data_long$variable)))
out_variables <- c("int_yr")

###  Send to JAGS
mc3     <- jags.model(file=textConnection(my_model), data=mydat, n.chains=1)
update(mc3, n.iter = 10000)
mc3.out <- coda.samples(model=mc3, variable.names=out_variables, n.iter=5000)
# plot(mc3.out)
mfit <- as.data.frame(as.matrix(mc3.out,chains=TRUE))
mfit$iteration <- 1:nrow(mfit)



####
####  CALCULATE POSTERIOR CORRELATIONS OF RANDOM YEAR EFFECTS ----
####
mfit_long <- melt(mfit, id.vars = c("CHAIN","iteration")) %>%
  separate(variable, c("variable","species"), sep = ",") %>%
  separate(variable, c("variable", "year"), sep="\\[") %>%
  separate(species, c("species", "extra"), sep=-2)
mfit_long <- mfit_long[,-which(colnames(mfit_long)=="extra")]

corr_snoweff <- mfit_long%>%
  spread(key=species, value=value, fill=0) %>%
  group_by(iteration) %>%
  summarise(correlation = cor(`1`,`2`))

mycol <- viridis(1, begin = 0.5)
ggplot(corr_snoweff, aes(x = correlation))+
  geom_histogram(color="#EFEFEF", fill=mycol, bins = 30)+
  xlab("Correlation of Random Year Effects")+
  ylab("Frequency")+
  my_theme
ggsave("post_correlations.png", width = 3, height = 2.5, units="in", dpi=300)



####
####  PLOT POSTERIOR DISTRIBUTIONS OF SNOW EFFECT ----
####
out_variables <- c("beta2")
mc3     <- jags.model(file=textConnection(my_model), data=mydat, n.chains=1)
update(mc3, n.iter = 10000)
mc3.out <- coda.samples(model=mc3, variable.names=out_variables, n.iter=5000)
# plot(mc3.out)
mfit <- as.data.frame(as.matrix(mc3.out,chains=TRUE))
mfit$iteration <- 1:nrow(mfit)

mfit_long <- melt(mfit, id.vars = c("CHAIN","iteration")) %>%
  separate(variable, c("variable","species"), sep = "\\[") %>%
  separate(species, c("species", "extra"), sep=-2)
mfit_long <- mfit_long[,-which(colnames(mfit_long)=="extra")]

ggplot(mfit_long, aes(x = value, fill=species))+
  geom_histogram(color="#EFEFEF", bins=40)+
  xlab("Snow Depth Effect")+
  ylab("Frequency")+
  scale_fill_viridis(begin=0.2, end=0.7, discrete=TRUE)+
  guides(fill=FALSE)+
  my_theme
ggsave("post_snoweffect.png", width = 3, height = 2.5, units="in", dpi=300)



####
####  SIMULATE COMMUNITY WITH AND WITHOUT SNOW EFFECT ----
####

### Sample model, tracing all parameters
out_variables <- c("int_yr", "beta2", "beta1")
mc3     <- jags.model(file=textConnection(my_model), data=mydat, n.chains=1)
update(mc3, n.iter = 10000)
mc3.out <- coda.samples(model=mc3, variable.names=out_variables, n.iter=5000)
mfit <- as.data.frame(as.matrix(mc3.out,chains=TRUE))
mfit$iteration <- 1:nrow(mfit)

mfit_long <- melt(mfit, id.vars = c("CHAIN","iteration")) 
mfit_long$variable <- as.character(mfit_long$variable)
mfit_long[which(mfit_long$variable=="beta1"),"variable"] <- "beta1[0,0]"
mfit_long[which(mfit_long$variable=="beta2[1]"),"variable"] <- "beta2[0,1]"
mfit_long[which(mfit_long$variable=="beta2[2]"),"variable"] <- "beta2[0,2]"


mfit_decomp <- mfit_long %>%
  separate(variable, c("variable","species"), sep = ",") %>%
  separate(variable, c("variable", "year"), sep="\\[") %>%
  separate(species, c("species", "extra"), sep=-2)
mfit_decomp <- mfit_decomp[,-which(colnames(mfit_decomp)=="extra")]

mfit_means <- mfit_decomp %>%
  group_by(variable, year, species) %>%
  summarise(mean_value = mean(value))
mfit_means$year <- as.numeric(mfit_means$year)
mfit_means$species <- as.numeric(mfit_means$species)

### Simulate with climate and random years
nsims <- 100
nspp <- 2
yrs <- sample(c(1:9), size = nsims, replace = TRUE)
yhat <- matrix(nrow = nsims, ncol=nspp)
yhat[1,] <- 1
for(j in 1:nspp){
  for(i in 2:nsims){
    beta_dd <- as.numeric(subset(mfit_means, variable=="beta1")["mean_value"])
    snow_effs <- as.numeric(as.matrix(subset(mfit_means, variable=="beta2")["mean_value"]))
    yeardat <- subset(mfit_means, year==yrs[i])
    intercepts <- as.numeric(as.matrix(subset(yeardat, variable=="int_yr")["mean_value"]))
    
    year_snow <- snow_sample[[yrs[i]]]
    snow_row <- year_snow[sample(1:8, 1)]
  
    yhat[i,j] <- intercepts[j] + beta_dd*yhat[i-1,j] + snow_effs[j]*snow_row
  }
}

synch_snow_yrs <- community.sync(yhat)[[1]]

### Simulate with climate only
yhat <- matrix(nrow = nsims, ncol=nspp)
yhat[1,] <- 1
for(j in 1:nspp){
  for(i in 2:nsims){
    beta_dd <- as.numeric(subset(mfit_means, variable=="beta1")["mean_value"])
    snow_effs <- as.numeric(as.matrix(subset(mfit_means, variable=="beta2")["mean_value"]))
    yeardat <- subset(mfit_means, year==yrs[i])
    intercepts <- mean(as.numeric(as.matrix(subset(yeardat, variable=="int_yr")["mean_value"])))

    year_snow <- snow_sample[[yrs[i]]]
    snow_row <- year_snow[sample(1:8, 1)]

    yhat[i,j] <- intercepts + beta_dd*yhat[i-1,j] + snow_effs[j]*snow_row
  }
}
synch_snow <- community.sync(yhat)[[1]]

### Simulate with random years only
yhat <- matrix(nrow = nsims, ncol=nspp)
yhat[1,] <- 1
for(j in 1:nspp){
  for(i in 2:nsims){
    beta_dd <- as.numeric(subset(mfit_means, variable=="beta1")["mean_value"])
    snow_effs <- mean(as.numeric(as.matrix(subset(mfit_means, variable=="beta2")["mean_value"])))
    yeardat <- subset(mfit_means, year==yrs[i])
    intercepts <- as.numeric(as.matrix(subset(yeardat, variable=="int_yr")["mean_value"]))
    
    year_snow <- snow_sample[[yrs[i]]]
    snow_row <- year_snow[sample(1:8, 1)]
    
    yhat[i,j] <- intercepts[j] + beta_dd*yhat[i-1,j] + snow_effs*snow_row
  }
}
synch_yrs <- community.sync(yhat)[[1]]


plot_dat <- data.frame(synchrony = c(synch_snow_yrs, synch_snow, synch_yrs),
                       simulation = c("ASnowYrs","BSnow","CYears"))

ggplot(plot_dat, aes(x=simulation, y=synchrony, fill=simulation))+
  geom_bar(stat="identity", width=0.5)+
  scale_y_continuous(limits=c(0,0.75))+
  scale_fill_viridis(end=0.8, discrete=TRUE)+
  ylab("Synchrony")+
  xlab("Simulation")+
  scale_x_discrete(labels = c("","",""))+
  guides(fill=FALSE)+
  my_theme
ggsave("synchrony_sims.png", width = 3, height = 2.5, units="in", dpi=300)
