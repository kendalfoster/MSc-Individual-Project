# Testing for implementing constant drift rate response time algorithms
    # tests for accuracy (consistency) and speed

source("PDF-Approximation/benchmark.r") # this contains other necessary imports





################################################################################
####################### Accuracy (Consistency among methods)

# Sample parameters
rt = seq(0.1, 3, by=0.01)
v = rnorm(1, 0, 2)
a = runif(1, 0.5, 3)
w = runif(1, 0.4, 0.6) # w = z/a, so needs to be rescaled for rtdists
t0 = 0.0001#runif(1, 0, 0.5)
eps = sqrt(.Machine$double.eps)


# Calculate the various densities
RTDists_R = ddiffusion(rt, rep("lower", length(rt)), v=v, a=a, z=w*a, t0=t0) # rtdists
BGK2014_R = fs14_R(t=rt-t0, v=v, a=a, w=w, eps=eps) # BGK's 2014 R Code
RWiener_R = dwiener(rt, resp=rep("lower", length(rt)), delta=v, alpha=a, beta=w, tau=t0)
fs_eps_14 = fs_eps_2014(rt=rt, v=v, a=a, w=w, t0=t0, eps=eps)
fs_eps_17 = fs_eps_2017(rt=rt, v=v, a=a, w=w, t0=t0, eps=eps)
fs_Nav_14 = fs_Nav_2014(rt=rt, v=v, a=a, w=w, t0=t0, eps=eps)
fs_Nav_17 = fs_Nav_2017(rt=rt, v=v, a=a, w=w, t0=t0, eps=eps)
fs_BGK_14 = fs_BGK_2014(rt=rt, v=v, a=a, w=w, t0=t0, eps=eps)
fs_BGK_17 = fs_BGK_2017(rt=rt, v=v, a=a, w=w, t0=t0, eps=eps)
fl_eps_00 = fl_eps(     rt=rt, v=v, a=a, w=w, t0=t0, eps=eps)
fl_Nav_00 = fl_Nav(     rt=rt, v=v, a=a, w=w, t0=t0, eps=eps)
fb_Nav_Nav= f_Nav_Nav(  rt=rt, v=v, a=a, w=w, t0=t0, eps=eps)
fb_BGK_Nav= f_BGK_Nav(  rt=rt, v=v, a=a, w=w, t0=t0, eps=eps)


# Visual inspection
plot(  rt, RTDists_R, pch=8, cex=.15,  col='black')
points(rt, BGK2014_R, pch=16, col='green')
points(rt, RWiener_R, pch=4, col='orange')
points(rt, fs_eps_14, pch=1, cex=.25,  col='black')
points(rt, fs_eps_17, pch=1,  col='blue')
points(rt, fs_Nav_14, pch=17, cex=.25, col='red')
points(rt, fs_Nav_17, pch=17, col='blue')
points(rt, fs_BGK_14, pch=16, col='red')
points(rt, fs_BGK_17, pch=16, col='blue')
points(rt, fl_eps_00, pch=5, cex=.25, col='pink')
# points(rt, fl_Nav_00, pch=17, cex=.5, col='pink')
# points(rt, fb_Nav_Nav,pch=20, cex=.5, col='green')
# points(rt, fb_BGK_Nav,pch=20, cex=.5, col='blue')


# Check actual consistency; from Code/Almost_Equal.cpp
almost_equal(RTDists_R, BGK2014_R, thresh=1e-7)
almost_equal(RTDists_R, RWiener_R, thresh=1e-7)
almost_equal(RTDists_R, fs_eps_14, thresh=1e-7)
almost_equal(RTDists_R, fs_eps_17, thresh=1e-7)
almost_equal(RTDists_R, fs_Nav_14, thresh=1e-7)
almost_equal(RTDists_R, fs_Nav_17, thresh=1e-7)
almost_equal(RTDists_R, fs_BGK_14, thresh=1e-7)
almost_equal(RTDists_R, fs_BGK_17, thresh=1e-7)
almost_equal(RTDists_R, fl_eps_00, thresh=1e-6)





################################################################################
####################### Microbenchmark

### Run benchmark tests and save output
# Set parameter exploration vectors
RT = c(0.1, seq(0.5, 3, by = 0.5))
V = seq(-6, 6, by = 0.5)
A = seq(0.5, 5, by = 0.5)
W = seq(0.2, 0.8, by = 0.1)
t0 = 0.0001
eps = sqrt(.Machine$double.eps)

# Small time
cnst_sm <- rt_benchmark(RT=RT, V=V, A=A, W=W, t0=t0, size="small",
                        RTvec=FALSE, times=1000, unit="us")
saveRDS(cnst_sm, file = "PDF-Approximation/Benchmark-Results/cnst_sm.Rds")
cnst_sm <- readRDS("PDF-Approximation/Benchmark-Results/cnst_sm.Rds")

cnst_sm_vec <- rt_benchmark(RT=RT, V=V, A=A, W=W, t0=t0, size="small",
                        RTvec=TRUE, times=1000, unit="us")
saveRDS(cnst_sm_vec, file = "PDF-Approximation/Benchmark-Results/cnst_sm_vec.Rds")
cnst_sm_vec <- readRDS("PDF-Approximation/Benchmark-Results/cnst_sm_vec.Rds")

# Large time
cnst_lg <- rt_benchmark(RT=RT, V=V, A=A, W=W, t0=t0,
                        size="large", times=1000, unit="us")
saveRDS(cnst_lg, file = "PDF-Approximation/Benchmark-Results/cnst_lg.Rds")
cnst_lg <- readRDS("PDF-Approximation/Benchmark-Results/cnst_lg.Rds")

# Small and large time
cnst <- rbind(subset(cnst_sm, FuncName %in% c("RTDists_R", "BGK2014_R",
                                              "RWiener_R",
                                              "fs_eps_14", "fs_eps_17",
                                              "fs_Nav_14", "fs_Nav_17",
                                              "fs_BGK_14", "fs_BGK_17")),
              subset(cnst_lg, FuncName %in% c("large_eps", "large_Nav",
                                              "both_Nav_Nav", "both_BGK_Nav")))
saveRDS(cnst, file = "PDF-Approximation/Benchmark-Results/cnst.Rds")
cnst <- readRDS("PDF-Approximation/Benchmark-Results/cnst.Rds")

# Small and large time C++ code
cnst_c <- rbind(subset(cnst_sm, FuncName %in% c("fs_eps_14", "fs_eps_17",
                                              "fs_Nav_14", "fs_Nav_17",
                                              "fs_BGK_14", "fs_BGK_17")),
              subset(cnst_lg, FuncName %in% c("large_eps", "large_Nav",
                                              "both_Nav_Nav", "both_BGK_Nav")))
saveRDS(cnst_c, file = "PDF-Approximation/Benchmark-Results/cnst_c.Rds")
cnst_c <- readRDS("PDF-Approximation/Benchmark-Results/cnst_c.Rds")

### Visualizations
anim_bm(cnst_sm, filepath="PDF-Approximation/Benchmark-Results/",
        filename="cnst_sm", t_trunc=c(125,125,125,125), fac=2)
anim_bm(cnst_sm_vec, filepath="PDF-Approximation/Benchmark-Results/",
        filename="cnst_sm_vec", t_trunc=c(125,125,125,125), fac=2)
anim_bm(cnst_lg, filepath="PDF-Approximation/Benchmark-Results/",
        filename="cnst_lg", t_trunc=c(125,125,125,125))
anim_bm(cnst_c, filepath="PDF-Approximation/Benchmark-Results/",
        filename="cnst_c", t_trunc=c(125,125,125,125))



cnst_c <- readRDS("PDF-Approximation/Benchmark-Results/cnst_c.Rds")
filepath="PDF-Approximation/Benchmark-Results/"
filename="cnst_c"
fac <- 2
bm <- cnst
t_idx <- match("FuncName", colnames(bm))
temp <- bm[,-seq_len(t_idx)]
if (fac != 1) {
  niter <- ncol(temp)
  nf <- niter/fac
  bm[,(seq_len(niter/fac)+t_idx)] <- temp[,seq_len(niter/fac)]/(fac*1000)
  for (i in seq_len(fac-1)) {
    bm[,(seq_len(niter/fac)+t_idx)] <- bm[,(seq_len(niter/fac)+t_idx)] +
                                       temp[,(seq_len(niter/fac)+(i*nf))]/(fac*1000)
  }
  bm <- bm[,seq_len(t_idx+niter/fac)]
}
else {
  bm[,-seq_len(t_idx)] <- bm[, -seq_len(t_idx)]/1000
}
mbm <- melt(bm, measure.vars = -seq_len(t_idx),
                 variable.name="iter", value.name="time")

#### make an animated gif for each varied parameter
# RT
anim <- ggplot(subset(mbm, time<15), aes(x=time, color=FuncName, fill=FuncName)) +
          geom_density(alpha=0.4) +
          transition_states(W,
                            transition_length = 2,
                            state_length = 1) +
          ease_aes('cubic-in-out') + # Slow start and end for a smoother look
          labs(title="Density of Microbenchmark Results",
               x="Time (ms)", y="Density",
               subtitle = "w: {closest_state}")
anim_save(paste0(filepath, filename, "_w.gif"), animation=anim)

#### make a still png for each varied parameter
# RT
ggplot(subset(mbm, time<15 & RT==RT[7]), aes(x=time, color=FuncName, fill=FuncName)) +
          geom_density(alpha=0.4) +
          labs(title="Density of Microbenchmark Results",
               x="Time (ms)", y="Density",
               subtitle = "rt: 3")
ggsave("PDF-Approximation/Images/rt/6.png")

#### Box/Violin Plot to show densities
mbm$FuncName2 <- factor(mbm$FuncName, levels=unique(mbm$FuncName))
ggplot(subset(mbm, time < 200), aes(x=FuncName2, y=time, color=FuncName2, fill=FuncName2)) +
          geom_violin(trim=FALSE) +
          geom_boxplot(width=0.1, fill="white") +
          labs(title="Density of Microbenchmark Results",
               x="Method", y="Time (ms)") +
          theme(axis.text.x = element_text(angle = 45, hjust=1)) +
          theme(legend.position="none")
ggsave("PDF-Approximation/Images/violin_master.png")


#### Determine significant differences between algorithms

fs_eps_14 <- melt(subset(cnst, FuncName == "fs_eps_14")[,-seq_len(5)])$value
fs_eps_17 <- melt(subset(cnst, FuncName == "fs_eps_17")[,-seq_len(5)])$value
fs_Nav_14 <- melt(subset(cnst, FuncName == "fs_Nav_14")[,-seq_len(5)])$value
fs_Nav_17 <- melt(subset(cnst, FuncName == "fs_Nav_17")[,-seq_len(5)])$value
fs_BGK_14 <- melt(subset(cnst, FuncName == "fs_BGK_14")[,-seq_len(5)])$value
fs_BGK_17 <- melt(subset(cnst, FuncName == "fs_BGK_17")[,-seq_len(5)])$value
large_eps <- melt(subset(cnst, FuncName == "large_eps")[,-seq_len(5)])$value
large_Nav <- melt(subset(cnst, FuncName == "large_Nav")[,-seq_len(5)])$value
both_Nav_Nav <- melt(subset(cnst, FuncName == "both_Nav_Nav")[,-seq_len(5)])$value
both_BGK_Nav <- melt(subset(cnst, FuncName == "both_BGK_Nav")[,-seq_len(5)])$value
RWiener_R <- melt(subset(cnst, FuncName == "RWiener_R")[,-seq_len(5)])$value
BGK2014_R <- melt(subset(cnst, FuncName == "BGK2014_R")[,-seq_len(5)])$value
RTDists_R <- melt(subset(cnst, FuncName == "RTDists_R")[,-seq_len(5)])$value

all_times <- data.frame(fs_eps_14, fs_eps_17, fs_Nav_14, fs_Nav_17, fs_BGK_14, fs_BGK_17, large_eps, large_Nav, both_Nav_Nav, both_BGK_Nav, RWiener_R, BGK2014_R, RTDists_R)
saveRDS(all_times, file = "PDF-Approximation/Benchmark-Results/all_times.Rds")
all_times <- readRDS("PDF-Approximation/Benchmark-Results/all_times.Rds")

FuncNames <- c("fs_eps_14", "fs_eps_17", "fs_Nav_14", "fs_Nav_17", "fs_BGK_14", "fs_BGK_17", "large_eps", "large_Nav", "both_Nav_Nav", "both_BGK_Nav", "RWiener_R", "BGK2014_R", "RTDists_R")

diffs <- data.frame(matrix(0, ncol=13, nrow=13))
colnames(diffs) <- FuncNames
rownames(diffs) <- FuncNames

for (i in 1:(length(FuncNames)-1)) {
  f1 <- mean(subset(all_times, select=FuncNames[i])[[1]])
  for (j in (i+1):length(FuncNames)) {
    f2 <- mean(subset(all_times, select=FuncNames[j])[[1]])
    diffs[i,j] <- 200*abs(f1-f2)/(f1+f2)
    diffs[j,i] <- diffs[i,j]
  }
}

library("tidyverse")
dt2 <- diffs %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

heatmap <- ggplot(dt2, aes(x = factor(rowname, levels=unique(dt2$rowname)),
                           y = factor(colname, levels=rev(unique(dt2$colname))),
                           fill = value)) +
             geom_tile() +
             labs(title="Percent Differences in Microbenchmark Results",
                  x=NULL, y=NULL, fill="Percent Difference") +
             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("PDF-Approximation/Images/heatmap.png")
