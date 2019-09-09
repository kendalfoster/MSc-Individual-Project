### Plotting functions
source("Bayesian-Hierarchical-Model/fitting_plots.r")

### Rough Timing
library('tictoc')

########## Stan- Write and Compile Model #######################################
library('rstan')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

exp_log_model_stanc <- stanc(file = 'Bayesian-Hierarchical-Model/exp_logistic_model_h.stan',
                             model_name = "exponential_logistic_model")
exp_log_model <- stan_model(stanc_ret = exp_log_model_stanc)
saveRDS(exp_log_model, 'Bayesian-Hierarchical-Model/exp_logistic_model_h.Rds')
exp_log_model <- readRDS('Bayesian-Hierarchical-Model/exp_logistic_model_h.Rds')

########## Define Functions ####################################################
exp_log <- function(rt, V, E, A, D, L, H) {
  return( (V + E*exp(-A*rt)) * (1 - D/(1 + exp(-L*(rt-H)))) )
}
sim_exp_log <- function(rt, V, E, A, D, L, H, var) {
  return( exp_log(rt=rt, E=E, A=A, V=V, D=D, L=L, H=H)
          + rnorm(length(rt), 0, var) )
}


########## Test Fitting (with SIMULATED data) ##################################
# Prepare Simulated Data
K <- 2
N <- 100
NTI <- rep(N, K)
NK <- N*K
Y0 <- rep(NA, NK)
Y1 <- rep(NA, NK)
V <- c(1.1, 1.5, .75)
E <- c(.25, .25, .25)
A <- c(.25, .25, .25)
D <- c(.4, 0, .5)
L <- c(.4, .4, .4)
H <- c(50, 60, 30)
st <- 0
for (i in seq_len(K)) {
  x <- seq_len(NTI[i])
  Y0[(st+1):(st+NTI[i])] <- sim_exp_log(x, V=V[i], E=E[i], A=A[i], D=0, L=L[i], H=H[i], var=0.15)
  Y1[(st+1):(st+NTI[i])] <- sim_exp_log(x, V=V[i], E=E[i], A=A[i], D=D[i], L=L[i], H=H[i], var=0.15)
  st = st + NTI[i]
}

# Fit Simulated Data
tic()
sim_fit <- summary(sampling(
  exp_log_model,
  data = list('K'=K, 'NTI'=NTI, 'NK'=NK, 'Y0'=Y0, 'Y1'=Y1),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10)
))$summary
toc()
saveRDS(sim_fit, "Bayesian-Hierarchical-Model/SampleFits/sim_fit_h_two.RDS")

# Plot Fits
plot_fits_sim_gg(K, NTI, V, E, A, D, L, H, Y0, Y1, fits=sim_fit[,1], save=TRUE)




########## Test Fitting (with REAL data) #######################################
# Visually examine data (optional)
data <- read.csv(file = "Bayesian-Hierarchical-Model/simplified-exp-1.csv",
                 header = TRUE, sep = ",")
sub_ids <- unique(data$subject_id)
par(mfrow=c(length(sub_ids),1))
for (i in 1:length(sub_ids)) {
  s0 <- subset(data, subject_id == sub_ids[i] & is_predictable==0)$rt/1000
  s1 <- subset(data, subject_id == sub_ids[i] & is_predictable==1)$rt/1000
  n <- min(length(s0), length(s1))
  s0 <- s0[seq_len(n)]
  s1 <- s1[seq_len(n)]
  x <- seq_len(n)
  plot(x, s1, pch=20, col='black', ylab = paste("Subject", i, sep=" "))
  points(x, s0, pch=20, col='gray50')
}




# Sample subject who learned
data <- read.csv(file = "Bayesian-Hierarchical-Model/simplified-exp-1.csv",
                 header = TRUE, sep = ",")
sub_ids <- unique(data$subject_id)[78]
df <- subset(data, subject_id==sub_ids)
label <- vector()
for (i in 1:nrow(df)) {
  if (df$is_predictable[i] == 0) {
    label <- c(label, "Non-Learned Data")
  } else {
    label <- c(label, "Learned Data")
  }
}
df$label <- label
ggplot(df) +
  geom_point(aes(x=t, y=rt/1000, color=label), shape=20) +
  scale_color_manual(values=c("black", "gray60")) +
  labs(title = "Example Learning Data",
       x = "Trial Number", y = "Response Time (sec)",
       subtitle = "Subject 78") +
   theme(legend.title=element_blank())
ggsave("Paper/Images/RTLearn.png")

# Sample subject who learned with fit
data <- read.csv(file = "Bayesian-Hierarchical-Model/simplified-exp-1.csv",
                 header = TRUE, sep = ",")
sub_ids <- unique(data$subject_id)[c(78, 11, 8)]
KK <- length(sub_ids)
gs_idx <- vector()
NTI <- rep(0, KK)
Y0 <- vector()
Y1 <- vector()

for (i in seq_len(KK)) {
  temp <- subset(data, subject_id==sub_ids[i])
  temp0 <- subset(temp, is_predictable==0)$rt/1000
  temp1 <- subset(temp, is_predictable==1)$rt/1000
  mm <- min(length(temp0), length(temp1))
  if (mm >= 50) {
    NTI[i] <- mm
    Y0 <- c(Y0, temp0[seq_len(mm)])
    Y1 <- c(Y1, temp1[seq_len(mm)])
    gs_idx <- c(gs_idx, i)
  }
}
good_sub_ids <- sub_ids[gs_idx]
NTI <- NTI[NTI > 0]
NK <- sum(NTI)
K <- length(good_sub_ids)
K == length(NTI)
NK == length(Y0)


# Fit the Data
tic()
real_fit <- summary(sampling(
  exp_log_model,
  data = list('K'=K, 'NTI'=NTI, 'NK'=NK, 'Y0'=Y0, 'Y1'=Y1),
  refresh = FALSE, chains = 1, iter = 250, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10)
))$summary
toc()
saveRDS(real_fit, "Bayesian-Hierarchical-Model/SampleFits/real_fit_paper_three.Rds")

# Plot Fits
plot_fits_gg(K, NTI, Y0, Y1, fits=real_fit[,1], save=TRUE)



### Learning Occurs
# Prepare the data
data <- read.csv(file = "Bayesian-Hierarchical-Model/simplified-exp-1.csv",
                 header = TRUE, sep = ",")
sub_ids <- unique(data$subject_id)
sub_ids <- sub_ids[c(8, 78, 13)]
KK <- length(sub_ids)
gs_idx <- vector()
NTI <- rep(0, KK)
Y0 <- vector()
Y1 <- vector()

for (i in seq_len(KK)) {
  temp <- subset(data, subject_id==sub_ids[i])
  temp0 <- subset(temp, is_predictable==0)$rt/1000
  temp1 <- subset(temp, is_predictable==1)$rt/1000
  mm <- min(length(temp0), length(temp1))
  if (mm >= 50) {
    NTI[i] <- mm
    Y0 <- c(Y0, temp0[seq_len(mm)])
    Y1 <- c(Y1, temp1[seq_len(mm)])
    gs_idx <- c(gs_idx, i)
  }
}
good_sub_ids <- sub_ids[gs_idx]
NTI <- NTI[NTI > 0]
NK <- sum(NTI)
K <- length(good_sub_ids)
K == length(NTI)
NK == length(Y0)


# Fit the Data
tic()
real_fit <- summary(sampling(
  exp_log_model,
  data = list('K'=K, 'NTI'=NTI, 'NK'=NK, 'Y0'=Y0, 'Y1'=Y1),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10)
))$summary
toc()
saveRDS(real_fit, "Bayesian-Hierarchical-Model/SampleFits/real_fit_three.RDS")

# Plot Fits
subV <- c(.9, .75, .78)
subE <- c(.5, .5, .25)
subA <- c(.35, .5, .25)
subD <- c(.4, .45, .4)
subL <- c(.5, .9, .4)
subH <- c(15, 49, 40)
plot_fits_h(K, NTI, Y0, Y1, fits=real_fit[,1], P=TRUE)





################################################################################
########## Fitting ALL of the REAL data ########################################

##### Prepare the data
data <- read.csv(file = "Bayesian-Hierarchical-Model/simplified-exp-1.csv",
                 header = TRUE, sep = ",")
sub_ids <- unique(data$subject_id)
KK <- length(sub_ids)
gs_idx <- vector()
dif <- vector()
NTI <- rep(0, KK)
Y0 <- vector()
Y1 <- vector()

for (i in seq_len(KK)) {
  temp <- subset(data, subject_id==sub_ids[i])
  temp0 <- subset(temp, is_predictable==0)$rt/1000
  temp1 <- subset(temp, is_predictable==1)$rt/1000
  mm <- min(length(temp0), length(temp1))
  if (mm >= 50) {
    NTI[i] <- mm
    Y0 <- c(Y0, temp0[seq_len(mm)])
    Y1 <- c(Y1, temp1[seq_len(mm)])
    gs_idx <- c(gs_idx, i)
    dif <- c(dif, temp$subject_condition[1])
  }
}
good_sub_ids <- sub_ids[gs_idx]
NTI <- NTI[NTI > 0]
NK <- sum(NTI)
K <- length(good_sub_ids)

dif[dif==1] <- "Easy"
dif[dif==2] <- "Medium"
dif[dif==3] <- "Hard"
dif <- factor(dif, levels=c("Easy", "Medium", "Hard"))

K == length(NTI)
K == length(dif)
NK == length(Y0)


##### Fit the Data
tic()
real_fit <- summary(sampling(
  exp_log_model,
  data = list('K'=K, 'NTI'=NTI, 'NK'=NK, 'Y0'=Y0, 'Y1'=Y1),
  refresh = FALSE, chains = 2, iter = 250, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10)
))$summary
toc()
saveRDS(real_fit, "Bayesian-Hierarchical-Model/SampleFits/real_fit_ALL_by_par.RDS")
real_fit <- readRDS("Bayesian-Hierarchical-Model/SampleFits/real_fit_ALL_by_par.RDS")

# Reorganize the fitted data
nind <- K # number of individuals in the file
all_fits <- data.frame(matrix(NA, nrow=16*nind, ncol=10))
colnames(all_fits) <- c('mean', 'se_mean', 'sd', '2.5%', '25%', '50%', '75%', '97.5%', 'n_eff', 'Rhat')
varnames <- c('V', 'v_alpha', 'v_beta', 'E', 'A', 'P', 'D', 'd_alpha', 'd_beta', 'L', 'l_alpha', 'l_beta', 'H', 'h_loc', 'h_scale', 'sigma')
rnames <- c()
ctr <- 1
for (j in seq_len(nind)) {
  all_fits[ctr:(ctr+15),] <- real_fit[c(j, j+nind, j+2*nind, j+3*nind, j+4*nind, j+5*nind, j+6*nind, j+7*nind, j+8*nind, j+9*nind, j+10*nind, j+11*nind, j+12*nind, j+13*nind, j+14*nind, j+15*nind),]
  rnames <- c(rnames, paste0(varnames, '[', j, ']'))
  ctr = ctr + 15
}
rownames(all_fits) <- rnames
saveRDS(all_fits, "Bayesian-Hierarchical-Model/SampleFits/real_fit_ALL_by_sub.RDS")

### Plot Fits
plot_fits(k=K, NTI=NTI, y0=Y0, y1=Y1, fits=real_fit[,1])

### Plot distributions of fitted parameters grouped by difficulty
library("ggplot2")
library("latex2exp")
P_vec <- real_fit[(5*K+1):(6*K),1]

# V, mean non-learning response time
V_df <- data.frame(var=real_fit[(0*K+1):(1*K),1], dif=dif)
V_plot <- ggplot(V_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         scale_color_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = "V, mean non-learning response time") +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/V_plot.png")

V_alpha_df <- data.frame(var=real_fit[(1*K+1):(2*K),1], dif=dif)
V_alpha_plot <- ggplot(V_alpha_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         scale_color_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = TeX('$v_\\alpha$, hyperparameter for V')) +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/V_alpha_plot.png")

V_beta_df <- data.frame(var=real_fit[(2*K+1):(3*K),1], dif=dif)
V_beta_plot <- ggplot(V_beta_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         scale_color_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = TeX('$v_\\beta$, hyperparameter for V')) +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/V_beta_plot.png")

# E, scale of adaptation
E_df <- data.frame(var=real_fit[(3*K+1):(4*K),1], dif=dif)
E_plot <- ggplot(E_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         scale_color_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = "E, scale of adaptation") +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/E_plot.png")

# A, rate of adaptation
A_df <- data.frame(var=real_fit[(4*K+1):(5*K),1], dif=dif)
A_plot <- ggplot(A_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         scale_color_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = "A, rate of adaptation") +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/A_plot.png")

# P, probablity of learning
P_df <- data.frame(var=real_fit[(5*K+1):(6*K),1], dif=dif, )
P_prior <- data.frame(prior=rep(rbeta(1000, 0.05, 0.05),3),
                      lev=factor(c(rep("Easy",1000), rep("Medium",1000), rep("Hard",1000)), levels=c("Easy", "Medium", "Hard")))
P_plot <- ggplot() +
         geom_violin(data=P_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif), trim=TRUE) +
         scale_fill_manual(values=c("#ace8bd", "#48b066", "#005e1b")) +
         scale_color_manual(values=c("#ace8bd", "#48b066", "#005e1b")) +
         geom_boxplot(data=P_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif), width=0.1, fill="white") +
         # geom_violin(data=P_prior, aes(x=lev, y=prior, group=lev, fill=NA), size=0.5) +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = "P, probability of learning") +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/P_plot.png")

# D, amount of learning
D_df <- data.frame(var=real_fit[(6*K+1):(7*K),1][P_vec>0.5], dif=dif[P_vec>0.5])
D_plot <- ggplot(D_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#c4d7ff", "#6790e0", "#001c59")) +
         scale_color_manual(values=c("#c4d7ff", "#6790e0", "#001c59")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = "D, amount fr learning") +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/D_plot.png")

D_alpha_df <- data.frame(var=real_fit[(7*K+1):(8*K),1][P_vec>0.5], dif=dif[P_vec>0.5])
D_alpha_plot <- ggplot(D_alpha_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#c4d7ff", "#6790e0", "#001c59")) +
         scale_color_manual(values=c("#c4d7ff", "#6790e0", "#001c59")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = TeX('$d_\\alpha$, hyperparameter for D')) +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/D_alpha_plot.png")

D_beta_df <- data.frame(var=real_fit[(8*K+1):(9*K),1][P_vec>0.5], dif=dif[P_vec>0.5])
D_beta_plot <- ggplot(D_beta_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
        geom_violin(trim=TRUE) +
        scale_fill_manual(values=c("#c4d7ff", "#6790e0", "#001c59")) +
        scale_color_manual(values=c("#c4d7ff", "#6790e0", "#001c59")) +
        geom_boxplot(width=0.1, fill="white") +
        labs(title="Densities of Fitting Results",
             x="Difficulty Level", y="Fitted Value",
             subtitle = TeX('$d_\\beta$, hyperparameter for D')) +
        theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/D_beta_plot.png")

# L, learning rate
L_df <- data.frame(var=real_fit[(9*K+1):(10*K),1][P_vec>0.5], dif=dif[P_vec>0.5])
L_plot <- ggplot(L_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         scale_color_manual(values=c("#ebe7b0", "#bdb43c", "#756d00")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = "L, learning rate") +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/L_plot.png")

# H, onset of learning
H_df <- data.frame(var=real_fit[(12*K+1):(13*K),1][P_vec>0.5],dif=dif[P_vec>0.5])
H_plot <- ggplot(H_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#ffcfe9", "#ff78c1", "#bf084e")) +
         scale_color_manual(values=c("#ffcfe9", "#ff78c1", "#bf084e")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = "H, onset of learning") +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/H_plot.png")

H_loc_df <- data.frame(var=real_fit[(13*K+1):(14*K),1][P_vec>0.5],dif=dif[P_vec>0.5])
H_loc_plot <- ggplot(H_loc_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#ffcfe9", "#ff78c1", "#bf084e")) +
         scale_color_manual(values=c("#ffcfe9", "#ff78c1", "#bf084e")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = TeX('$h_{loc}$, hyperparameter for H')) +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/H_loc_plot.png")

H_scale_df <- data.frame(var=real_fit[(14*K+1):(15*K),1][P_vec>0.5],dif=dif[P_vec>0.5])
H_scale_plot <- ggplot(H_scale_df, aes(x=dif, y=var, group=dif, color=dif, fill=dif)) +
         geom_violin(trim=TRUE) +
         scale_fill_manual(values=c("#ffcfe9", "#ff78c1", "#bf084e")) +
         scale_color_manual(values=c("#ffcfe9", "#ff78c1", "#bf084e")) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Densities of Fitting Results",
              x="Difficulty Level", y="Fitted Value",
              subtitle = TeX('$h_{scale}$, hyperparameter for H')) +
         theme(legend.position="none")
ggsave("Bayesian-Hierarchical-Model/SampleFits/ViolinPlots/H_scale_plot.png")
