# Testing for implementing constant drift rate response time algorithms
    # tests for accuracy (consistency) and speed

source("PDF-Approximation/benchmark.r") # this contains other necessary imports
library("readr")


# Import Data
med <- read_csv("PDF-Approximation/medical_dm.csv")
RT <- subset(med, response == "non-blast")$rt


################################################################################
####################### Microbenchmark

### Run benchmark tests and save output
# Set parameter exploration vectors
V = seq(-6, 6, by = 0.5)
A = seq(0.5, 5, by = 0.5)
W = seq(0.2, 0.8, by = 0.1)
t0 = 0.0001
eps = sqrt(.Machine$double.eps)

# Small time
real_sm <- rt_benchmark(RT=RT, V=V, A=A, W=W, t0=t0, size="small",
                        RTvec=TRUE, times=1000, unit="us")
saveRDS(real_sm, file = "PDF-Approximation/Benchmark-Results/real_sm.Rds")
real_sm <- readRDS("PDF-Approximation/Benchmark-Results/real_sm.Rds")

# Large time
real_lg <- rt_benchmark(RT=RT, V=V, A=A, W=W, t0=t0, size="large",
                        RTvec=TRUE, times=1000, unit="us")
saveRDS(real_lg, file = "PDF-Approximation/Benchmark-Results/real_lg.Rds")
real_lg <- readRDS("PDF-Approximation/Benchmark-Results/real_lg.Rds")

# Small and large time
real <- rbind(subset(real_sm, FuncName %in% c("RTDists_R", "BGK2014_R",
                                              "RWiener_R",
                                              "fs_eps_14", "fs_eps_17",
                                              "fs_Nav_14", "fs_Nav_17",
                                              "fs_BGK_14", "fs_BGK_17")),
              subset(real_lg, FuncName %in% c("large_eps", "large_Nav",
                                              "both_Nav_Nav", "both_BGK_Nav")))
saveRDS(real, file = "PDF-Approximation/Benchmark-Results/real.Rds")
real <- readRDS("PDF-Approximation/Benchmark-Results/real.Rds")

real_cp <- rbind(subset(real, FuncName %in% c("BGK2014_R",
                                              "fs_eps_14", "fs_eps_17",
                                              "fs_Nav_14", "fs_Nav_17",
                                              "fs_BGK_14", "fs_BGK_17",
                                              "large_eps", "large_Nav",
                                              "both_Nav_Nav", "both_BGK_Nav")))
saveRDS(real_cp, file = "PDF-Approximation/Benchmark-Results/real_cp.Rds")
real_cp <- readRDS("PDF-Approximation/Benchmark-Results/real_cp.Rds")





####################### Visualizations
cnst <- readRDS("PDF-Approximationcnst.Rds")
filepath="PDF-Approximation/"
filename="cnst"
fac <- 2
bm <- real
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
} else {
  bm[,-seq_len(t_idx)] <- bm[, -seq_len(t_idx)]/1000
}
mbm <- melt(bm, measure.vars = -seq_len(t_idx),
                 variable.name="iter", value.name="time")

ggplot(subset(mbm, time<6 & A==A[9]), aes(x=FuncName, y=time, color=FuncName, fill=FuncName)) +
         geom_violin(trim=FALSE) +
         geom_boxplot(width=0.1, fill="white") +
         labs(title="Density of Microbenchmark Results",
              x="method", y="Time (ms)",
              subtitle = "a: 5") +
         theme(axis.text.x = element_text(angle = 45, hjust=1)) +
         theme(legend.position="none")
ggsave("PDF-Approximation/Images/a/9.png")



ggplot(subset(mbm, time<80000), aes(x=time, color=FuncName, fill=FuncName)) +
         geom_line(aes(y = 1 - ..y..), stat='ecdf') +
         labs(title="Complementary ECDF of Microbenchmark Results",
              x="Time (ms)", y="Complementary ECDF",
              color="Method",
              subtitle = "subtitle") +
         theme(legend.position="right")
ggsave("PDF-Approximation/Images/ecdf1.png")
