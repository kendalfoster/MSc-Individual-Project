# R functions for running Microbenchmark tests on the response time models

# Utility packages and related local code
library("Rcpp")
library("microbenchmark")
sourceCpp("PDF-Approximation/Almost_Equal.cpp")

# Distribution packages and local code
library("rtdists")
library("RWiener")
source("PDF-Approximation/BGK2014.r")
sourceCpp("PDF-Approximation/constant_dr.cpp")
sourceCpp("PDF-Approximation/variable_dr.cpp")

# Plotting packages
library("reshape2")
library("ggplot2")
library("gganimate")





####################### Constant Drift Rate, Small Time
rt_bm_cnst_sv_small <- function(RT, V, A, W, t0=0.0001, times=1000, unit="us"){
  nf <- 9 # number of functions being benchmarked
  nRT <- length(RT) # number of response times
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  fnames <- c('RTDists_R', 'RWiener_R', 'BGK2014_R', 'fs_eps_14', 'fs_eps_17',
              'fs_Nav_14', 'fs_Nav_17', 'fs_BGK_14', 'fs_BGK_17')

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol=times+5, nrow=nf*nRT*nV*nA*nW))
  colnames(mbm_res) <- c('RT', 'V', 'A', 'W', 'FuncName',
                         paste0("bm", seq_len(times)))
  start <- 1
  stop <- nf
  sl_time <- seq_len(times)+5

  # Loop through each combination of parameters and record microbenchmark results
  for (rt in 1:nRT) {
    # cat(paste("progress: ", rt, "/", nRT, " = ", round(rt*100/nRT), "%\n", sep=""))
    for (v in 1:nV) {
      for (a in 1:nA) {
        for (w in 1:nW) {
          mbm <- microbenchmark(
            RTDists_R = ddiffusion(RT[rt], "lower", v=V[v], a=A[a], z=W[w]*A[a], t0=t0),
            RWiener_R = dwiener(RT[rt], resp="lower", delta=V[v], alpha=A[a], beta=W[w], tau=t0),
            BGK2014_R = fs14_R(t=RT[rt]-t0, v=V[v], a=A[a], w=W[w], eps=eps),
            fs_eps_14 = fs_eps_2014(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fs_eps_17 = fs_eps_2017(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fs_Nav_14 = fs_Nav_2014(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fs_Nav_17 = fs_Nav_2017(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fs_BGK_14 = fs_BGK_2014(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            fs_BGK_17 = fs_BGK_2017(rt=RT[rt], v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
            times = times, unit = unit)
          # add the rt, v, a, w values and function names to the dataframe
          mbm_res[start:stop, 1] <- rep(RT[rt], nf)
          mbm_res[start:stop, 2] <- rep(V[v]  , nf)
          mbm_res[start:stop, 3] <- rep(A[a]  , nf)
          mbm_res[start:stop, 4] <- rep(W[w]  , nf)
          mbm_res[start:stop, 5] <- fnames
          # add the microbenchmark results to the dataframe
          for (i in 1:nf) {
            bm <- subset(mbm, expr==fnames[i])
            mbm_res[start+i-1, sl_time] <- bm$time
          }
          # iterate start and stop values
          start = start + nf
          stop = stop + nf
        }
      }
    }
  }
  return(mbm_res)
}

rt_bm_cnst_sv_small_vec <- function(RT, V, A, W, t0=0.0001, times=1000, unit="us"){
  nf <- 9 # number of functions being benchmarked
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  resp <- rep("lower", length(RT)) # responses for rtdists
  fnames <- c('RTDists_R', 'RWiener_R', 'BGK2014_R', 'fs_eps_14', 'fs_eps_17',
              'fs_Nav_14', 'fs_Nav_17', 'fs_BGK_14', 'fs_BGK_17')

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol=times+4, nrow=nf*nV*nA*nW))
  colnames(mbm_res) <- c('V', 'A', 'W', 'FuncName',
                         paste0("bm", seq_len(times)))
  start <- 1
  stop <- nf
  sl_time <- seq_len(times)+4

  # Loop through each combination of parameters and record microbenchmark results
  for (v in 1:nV) {
    for (a in 1:nA) {
      for (w in 1:nW) {
        mbm <- microbenchmark(
          RTDists_R = ddiffusion(RT, resp, v=V[v], a=A[a], z=W[w]*A[a], t0=t0),
          RWiener_R = dwiener(RT, resp=resp, delta=V[v], alpha=A[a], beta=W[w], tau=t0),
          BGK2014_R = fs14_R(t=RT-t0, v=V[v], a=A[a], w=W[w], eps=eps),
          fs_eps_14 = fs_eps_2014(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fs_eps_17 = fs_eps_2017(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fs_Nav_14 = fs_Nav_2014(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fs_Nav_17 = fs_Nav_2017(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fs_BGK_14 = fs_BGK_2014(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          fs_BGK_17 = fs_BGK_2017(rt=RT, v=V[v], a=A[a], w=W[w], t0=t0, eps=eps),
          times = times, unit = unit)
        # add the v, a, w values and function names to the dataframe
        mbm_res[start:stop, 1] <- rep(V[v]  , nf)
        mbm_res[start:stop, 2] <- rep(A[a]  , nf)
        mbm_res[start:stop, 3] <- rep(W[w]  , nf)
        mbm_res[start:stop, 4] <- fnames
        # add the microbenchmark results to the dataframe
        for (i in 1:nf) {
          bm <- subset(mbm, expr==fnames[i])
          mbm_res[start+i-1, sl_time] <- bm$time
        }
        # iterate start and stop values
        start = start + nf
        stop = stop + nf
      }
    }
  }
  return(mbm_res)
}



####################### Constant Drift Rate, Large Time
rt_bm_cnst_sv_large <- function(RT, V, A, W, t0=0.0001, times=1000, unit="us"){
  nf <- 5 # number of functions being benchmarked
  nRT <- length(RT) # number of response times
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  resp <- rep("lower", nRT) # necessary for rtdists
  fnames <- c('RTDists_R', 'large_Nav', 'large_eps', 'both_Nav_Nav', 'both_BGK_Nav')

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol=times+5, nrow=nf*nRT*nV*nA*nW))
  colnames(mbm_res) <- c('RT', 'V', 'A', 'W', 'FuncName',
                         paste0("bm", seq_len(times)))
  start <- 1
  stop <- nf
  sl_time <- seq_len(times)+5

  # Loop through each combination of parameters and record microbenchmark results
  for (rt in 1:nRT) {
    # cat(paste("progress: ", rt, "/", nRT, " = ", round(rt*100/nRT), "%\n", sep=""))
    for (v in 1:nV) {
      for (a in 1:nA) {
        for (w in 1:nW) {
          mbm <- microbenchmark(
            RTDists_R = ddiffusion(RT[rt], resp, v=V[v], a=A[a], z=W[w]*A[a], t0=0),
            large_Nav = fl_Nav(RT[rt], v=V[v], a=A[a], w=W[w], t0=0),
            large_eps = fl_eps(RT[rt], v=V[v], a=A[a], w=W[w], t0=0),
            both_Nav_Nav = f_Nav_Nav(RT[rt], v=V[v], a=A[a], w=W[w], t0=0),
            both_BGK_Nav = f_BGK_Nav(RT[rt], v=V[v], a=A[a], w=W[w], t0=0),
            times = times, unit = unit)
          # add the rt, v, a, w values and function names to the dataframe
          mbm_res[start:stop, 1] <- rep(RT[rt], nf)
          mbm_res[start:stop, 2] <- rep(V[v]  , nf)
          mbm_res[start:stop, 3] <- rep(A[a]  , nf)
          mbm_res[start:stop, 4] <- rep(W[w]  , nf)
          mbm_res[start:stop, 5] <- fnames
          # add the microbenchmark results to the dataframe
          for (i in 1:nf) {
            bm <- subset(mbm, expr==fnames[i])
            mbm_res[start+i-1, sl_time] <- bm$time
          }
          # iterate start and stop values
          start = start + nf
          stop = stop + nf
        }
      }
    }
  }
  return(mbm_res)
}

rt_bm_cnst_sv_large_vec <- function(RT, V, A, W, t0=0.0001, times=1000, unit="us"){
  nf <- 5 # number of functions being benchmarked
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  resp <- rep("lower", length(RT)) # necessary for rtdists
  fnames <- c('RTDists_R', 'large_Nav', 'large_eps', 'both_Nav_Nav', 'both_BGK_Nav')

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol=times+4, nrow=nf*nV*nA*nW))
  colnames(mbm_res) <- c('V', 'A', 'W', 'FuncName',
                         paste0("bm", seq_len(times)))
  start <- 1
  stop <- nf
  sl_time <- seq_len(times)+4

  # Loop through each combination of parameters and record microbenchmark results
  for (v in 1:nV) {
    for (a in 1:nA) {
      for (w in 1:nW) {
        mbm <- microbenchmark(
          RTDists_R = ddiffusion(RT, resp, v=V[v], a=A[a], z=W[w]*A[a], t0=0),
          large_Nav = fl_Nav(RT, v=V[v], a=A[a], w=W[w], t0=0),
          large_eps = fl_eps(RT, v=V[v], a=A[a], w=W[w], t0=0),
          both_Nav_Nav = f_Nav_Nav(RT, v=V[v], a=A[a], w=W[w], t0=0),
          both_BGK_Nav = f_BGK_Nav(RT, v=V[v], a=A[a], w=W[w], t0=0),
          times = times, unit = unit)
        # add the rt, v, a, w values and function names to the dataframe
        mbm_res[start:stop, 1] <- rep(V[v], nf)
        mbm_res[start:stop, 2] <- rep(A[a], nf)
        mbm_res[start:stop, 3] <- rep(W[w], nf)
        mbm_res[start:stop, 4] <- fnames
        # add the microbenchmark results to the dataframe
        for (i in 1:nf) {
          bm <- subset(mbm, expr==fnames[i])
          mbm_res[start+i-1, sl_time] <- bm$time
        }
        # iterate start and stop values
        start = start + nf
        stop = stop + nf
      }
    }
  }
  return(mbm_res)
}



####################### Variable Drift Rate
rt_bm_vary_sv <- function(RT, V, A, W, SV, t0=0.0001, times=1000, unit="us") {
  nf <- 2 # number of functions being benchmarked
  nRT <- length(RT) # number of response times
  nV <- length(V) # number of drift rates
  nA <- length(A) # number of boundary separations
  nW <- length(W) # number of starting points
  nSV <- length(SV) # number of drift rate variances
  resp <- rep("lower", nRT) # necessary for rtdists
  fnames <- c('RTDists_R', 'Vary_sv')

  # Initialize the dataframe to contain the microbenchmark results
  mbm_res <- data.frame(matrix(ncol=times+6, nrow=nf*nRT*nV*nA*nW*nSV))
  colnames(mbm_res) <- c('RT', 'V', 'A', 'W', 'SV', 'FuncName',
                         paste0("bm", as.character(seq_len(times))))
  start <- 1
  stop <- nf
  sl_time <- seq_len(times)+6

  # Loop through each combination of parameters and record microbenchmark results
  for (rt in 1:nRT) {
    # cat(paste("progress: ", rt, "/", nRT, " = ", round(rt*100/nRT), "%\n", sep=""))
    for (v in 1:nV) {
      for (a in 1:nA) {
        for (w in 1:nW) {
          for (sv in 1:nSV) {
            mbm <- microbenchmark(
              RTDists_R = ddiffusion(RT[rt], resp, v=V[v], a=A[a], z=W[w]*A[a],
                                     sv=SV[sv], t0=0),
              Vary_sv = f_sv(rt=RT[rt], v=V[v], a=A[a], w=W[w],
                             sv=SV[sv], eps=eps),
              times = times, unit = unit)
            # add the rt, v, a, w values and function names to the dataframe
            mbm_res[start:stop, 1] <- rep(RT[rt], nf)
            mbm_res[start:stop, 2] <- rep(V[v]  , nf)
            mbm_res[start:stop, 3] <- rep(A[a]  , nf)
            mbm_res[start:stop, 4] <- rep(W[w]  , nf)
            mbm_res[start:stop, 5] <- rep(SV[sv]  , nf)
            mbm_res[start:stop, 6] <- fnames
            # add the microbenchmark results to the dataframe
            for (i in 1:nf) {
              bm <- subset(mbm, expr==fnames[i])
              mbm_res[start+i-1, sl_time] <- bm$time
            }
            # iterate start and stop values
            start = start + nf
            stop = stop + nf
          }
        }
      }
    }
  }
  return(mbm_res)
}


####################### Wrapper for Benchmarks #################################
rt_benchmark <- function(RT, V, A, W, SV=0, t0=0.0001, size="sm", RTvec=FALSE,
                         times=1000, unit="us") {
  if (all(SV == 0)) { # constant drift rate
    if (size %in% c("lg", "large", "Large", "LARGE")) { # large time appx
      if (RTvec) {
        return(rt_bm_cnst_sv_large_vec(RT=RT, V=V, A=A, W=W, t0=t0,
               times=times, unit=unit))
      } else {
        return(rt_bm_cnst_sv_large(RT=RT, V=V, A=A, W=W, t0=t0,
               times=times, unit=unit))
      }
    } else { # small time appx
      if (RTvec) { # benchmark all rt's as a vector
        return(rt_bm_cnst_sv_small_vec(RT=RT, V=V, A=A, W=W, t0=t0,
               times=times, unit=unit))
      } else { # benchmark each rt individually
        return(rt_bm_cnst_sv_small(RT=RT, V=V, A=A, W=W, t0=t0,
               times=times, unit=unit))
      }
    }
  } else { # variable drift rate
    return(rt_bm_vary_sv(RT=RT, V=V, A=A, W=W, SV=SV, t0=t0,
           times=times, unit=unit))
  }
}




####################### Animate Benchmark Distributions ########################
anim_bm <- function(bm, filepath, filename, t_trunc=rep(1e9,5), fac=1) {
  # convert to long format for ggplot
  t_idx <- match("FuncName", colnames(bm))
  temp <- bm[,-seq_len(t_idx)]
  niter <- ncol(temp)
  if (fac != 1) {
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
  anim <- ggplot(subset(mbm, time<t_trunc[1]), aes(x=time, color=FuncName, fill=FuncName)) +
            geom_density(alpha=0.4) +
            transition_states(RT,
                              transition_length = 2,
                              state_length = 1) +
            ease_aes('cubic-in-out') + # Slow start and end for a smoother look
            labs(title="Density of Microbenchmark Results",
                 x="Time (ms)", y="Density",
                 subtitle = "rt: {closest_state}")
  anim_save(paste0(filepath, filename, "_rt.gif"), animation=anim)
  # V
  anim <- ggplot(subset(mbm, time<t_trunc[2]), aes(x=time, color=FuncName, fill=FuncName)) +
            geom_density(alpha=0.4) +
            transition_states(V,
                              transition_length = 2,
                              state_length = 1) +
            ease_aes('cubic-in-out') + # Slow start and end for a smoother look
            labs(title="Density of Microbenchmark Results",
                 x="Time (ms)", y="Density",
                 subtitle = "v: {closest_state}")
  anim_save(paste0(filepath, filename, "_v.gif"), animation=anim)
  # A
  anim <- ggplot(subset(mbm, time<t_trunc[3]), aes(x=time, color=FuncName, fill=FuncName)) +
            geom_density(alpha=0.4) +
            transition_states(A,
                              transition_length = 2,
                              state_length = 1) +
            ease_aes('cubic-in-out') + # Slow start and end for a smoother look
            labs(title="Density of Microbenchmark Results",
                 x="Time (ms)", y="Density",
                 subtitle = "a: {closest_state}")
  anim_save(paste0(filepath, filename, "_a.gif"), animation=anim)
  # W
  anim <- ggplot(subset(mbm, time<t_trunc[4]), aes(x=time, color=FuncName, fill=FuncName)) +
            geom_density(alpha=0.4) +
            transition_states(W,
                              transition_length = 2,
                              state_length = 1) +
            ease_aes('cubic-in-out') + # Slow start and end for a smoother look
            labs(title="Density of Microbenchmark Results",
                 x="Time (ms)", y="Density",
                 subtitle = "w: {closest_state}")
  anim_save(paste0(filepath, filename, "_w.gif"), animation=anim)
  # SV (only for variable drift rate)
  if (t_idx==6) {
    anim <- ggplot(subset(mbm, time<t_trunc[5]), aes(x=time, color=FuncName, fill=FuncName)) +
              geom_density(alpha=0.4) +
              transition_states(SV,
                                transition_length = 2,
                                state_length = 1) +
              ease_aes('cubic-in-out') + # Slow start and end for a smoother look
              labs(title="Density of Microbenchmark Results",
                   x="Time (ms)", y="Density",
                   subtitle = "sv: {closest_state}")
    anim_save(paste0(filepath, filename, "_sv.gif"), animation=anim)
  }
}
