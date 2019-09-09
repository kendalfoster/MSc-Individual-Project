// Density functions for constant drift rate (i.e. sv=0)

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#define _USE_MATH_DEFINES
#include <math.h>




////////////////////////////////////////////////////////////////////////////////
/////////////////////// Small Time /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
///////////// Number of Terms Approximations

// BGK2014
// [[Rcpp::export]]
int ks_BGK(double t, double w, double eps) {
  double u_eps, arg, k1;
  int k;
  u_eps = std::min(-1.0, log(2*M_PI*t*t*eps*eps)); // Safe bound so that (idk that's all)
  arg = -t*(u_eps - sqrt(-2*u_eps - 2)); // sqrt(x) with x > 0
  k1 = (sqrt(2*t) - w)/2;
  if (arg > 0) { // If arg > 0, set k2 and calculate k
    double k2 = (sqrt(arg) - w)/2;
    k = ceil(std::max(k1, k2));
  }
  else { // Otherwise, we don't need k2
    k = ceil(k1);
  }
  return k;
}

// Navarro2009
// [[Rcpp::export]]
int ks_Nav(double t, double eps)
{
  double arg = 2*sqrt(2*M_PI*t)*eps;
  if (arg < 1) { // if error threshold is set low enough
    double ks = 2 + sqrt(-2*t*log(arg));
    return ceil(std::max(ks, sqrt(t)+1)); // ensure boundary conditions are met
  }
  return 2; // else error threshold set too high; return minimal kappa
}


////////////////////////////////////////////////////////////////////////////////
///////////// Infinite Sum Approximations

// BGK2014 style truncated sum
// [[Rcpp::export]]
double small_sum_2014(double t, double a, double w, int ks)
{
  double f;
  for (int j = ks; j > 0; j--) { // iterate through all ks
    f += (2*j + w) * exp(-a*a*(2*j + w)*(2*j + w)/(2*t))
        -(2*j - w) * exp(-a*a*(2*j - w)*(2*j - w)/(2*t));
  }
  f += w * exp(-a*a*w*w/(2*t)); // add case j=0
  return f;
}

// BGK2017 style truncated sum (uses ceil(k/2) as many terms as BGK2014)
// [[Rcpp::export]]
double small_sum_2017(double t, double a, double w, int ks)
{
  double sum = 0;
  int j = 0;
  double gamma = -a*a/(2*t);
  double rj;
  if (ks % 2 == 1) { // if ks is odd
    while (j < ks) { // start at j=0
      // j is even
      rj = j + w;
      sum += rj * exp(gamma*rj*rj);
      j++;
      // j is odd
      rj = j + 1 - w;
      sum -= rj * exp(gamma*rj*rj);
      j++;
    }
  }
  else { // ks is even
    sum += w * exp(gamma*w*w); // start at j=0
    j++;
    while (j < ks) { // start at j=1
      // j is odd
      rj = j + 1 - w;
      sum -= rj * exp(gamma*rj*rj);
      j++;
      // j is even
      rj = j + w;
      sum += rj * exp(gamma*rj*rj);
      j++;
    }
  }
  return sum;
}

// term < eps BGK2014 style truncated sum
// [[Rcpp::export]]
double small_sum_eps_14(double t, double a, double w, double eps)
{
  double gamma = -a*a/(2*t);
  double term = w * exp(gamma * w * w); // always include j=0 term
  double sum = term;
  int j = 0;
  while (fabs(term) > eps) {
    j++;
    term = (w + 2*j) * exp(gamma * (w + 2*j) * (w + 2*j))
          +(w - 2*j) * exp(gamma * (w - 2*j) * (w - 2*j));
    sum += term;
  }
  return sum;
}

// term < eps BGK2017 style truncated sum
// [[Rcpp::export]]
double small_sum_eps_17(double t, double a, double w, double eps)
{
  double term = w * exp(-a*a*w*w/(2*t)); // always include j=0 term
  double sum = term;
  int j = 0;
  double gamma = -a*a/(2*t);
  double rj;
  while (term > eps) { // start at j=1
    // j is odd
    j++;
    rj = j + 1 - w;
    sum -= rj * exp(gamma*rj*rj);
    // j is even
    j++;
    rj = j + w;
    term = rj * exp(gamma*rj*rj);
    sum += term;
  }
  return sum;
}


////////////////////////////////////////////////////////////////////////////////
///////////// Densities

// Use term < eps BGK2014style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fs_eps_2014(Rcpp::NumericVector rt, double v, double a, double w,
                           double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = DBL_EPSILON;
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
    out[i] = mult * small_sum_eps_14(t=t, a=a, w=w, eps=eps);
  }
  return out;
}

// Use term < eps BGK2017style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fs_eps_2017(Rcpp::NumericVector rt, double v, double a, double w,
                           double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = DBL_EPSILON;
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
    out[i] = mult * small_sum_eps_17(t=t, a=a, w=w, eps=eps);
  }
  return out;
}

// Use Navarro2009 number of terms for 2014 style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fs_Nav_2014(Rcpp::NumericVector rt, double v, double a, double w,
                                double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = DBL_EPSILON;
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  int ks;
  double t, mult;
  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    ks = ks_Nav(t=t, eps=eps); // get number of terms in sum approximation
    mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
    out[i] = mult * small_sum_2014(t=t, a=a, w=w, ks=ks);
  }
  return out;
}

// Use Navarro2009 number of terms for 2017 style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fs_Nav_2017(Rcpp::NumericVector rt, double v, double a, double w,
                                double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = DBL_EPSILON;
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  int ks;
  double t, mult;
  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    ks = ks_Nav(t=t, eps=eps); // get number of terms in sum approximation
    mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
    out[i] = mult * small_sum_2017(t=t, a=a, w=w, ks=ks);
  }
  return out;
}

// Use BGK2014 number of terms for 2014 style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fs_BGK_2014(Rcpp::NumericVector rt, double v, double a, double w,
                                double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = DBL_EPSILON;
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  int ks;
  double t, mult;
  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    ks = ks_BGK(t=t, w=w, eps=eps); // get number of terms in sum approximation
    mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t); // first two lines of Equation (1)
    out[i] = mult * small_sum_2014(t=t, a=a, w=w, ks=ks); // last line of Equation (1)
  }
  return out;
}

// Use BGK2014 number of terms for 2017 style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fs_BGK_2017(Rcpp::NumericVector rt, double v, double a, double w,
                                double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = DBL_EPSILON;
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  int ks;
  double t, mult;
  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    ks = ks_BGK(t=t, w=w, eps=eps); // get number of terms in sum approximation
    mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t); // first two lines of Equation (1)
    out[i] = mult * small_sum_2017(t=t, a=a, w=w, ks=ks); // last line of Equation (1)
  }
  return out;
}





////////////////////////////////////////////////////////////////////////////////
/////////////////////// Large Time /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
///////////// Number of Terms Approximations

// Navarro2009
// [[Rcpp::export]]
int kl_Nav(double t, double eps)
{
  double arg = M_PI*t*eps;
  if (arg < 1) { // if error threshold is set low enough
    double kl = sqrt(-2*log(arg)/(M_PI*M_PI*t));
    return ceil(std::max(kl, 1/(M_PI*sqrt(t))));
  }
  return ceil(1/(M_PI*sqrt(t)));
}


////////////////////////////////////////////////////////////////////////////////
///////////// Infinite Sum Approximations

// Navarro2009 style truncated sum
// [[Rcpp::export]]
double large_sum_Nav(double t, double a, double w, int kl)
{
  double gamma = -M_PI*M_PI*t/(2*a*a);
  double sum = 0.0;
  for (int j = 1; j <= kl; j++) {
    sum += j * sin(j*w*M_PI) * exp(gamma*j*j);
  }
  return sum;
  // double eps = DBL_EPSILON;
  // double gamma = -M_PI*M_PI*t/(2*a*a);
  // int j = 1;
  // double term = j * sin(j*w*M_PI) * exp(gamma*j*j); // always include j=1 term
  // double sum = term;
  // while (fabs(term) > eps) {
  //   j++;
  //   term = j * sin(j*w*M_PI) * exp(gamma*j*j);
  //   sum += term;
  // }
  // return sum;
}

// term < eps style truncated sum
// [[Rcpp::export]]
double large_sum_eps(double t, double a, double w, double eps)
{
  double gamma = -M_PI*M_PI*t/(2*a*a);
  int j = 1;
  double term = j * sin(j*w*M_PI) * exp(gamma*j*j); // always include j=1 term
  double sum = term;
  while (fabs(term) > eps) {
    j++;
    term = j * sin(j*w*M_PI) * exp(gamma*j*j);
    sum += term;
  }
  return sum;
}


////////////////////////////////////////////////////////////////////////////////
///////////// Densities

// Use term < eps style sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fl_eps(Rcpp::NumericVector rt, double v, double a, double w,
                           double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = DBL_EPSILON;
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    mult = M_PI * exp(-v*a*w - v*v*t/2) / (a*a);
    out[i] = mult * large_sum_eps(t=t, a=a, w=w, eps=eps);
  }
  return out;
}

// Use Navarro2009 number of terms for sum approximation
// [[Rcpp::export]]
Rcpp::NumericVector fl_Nav(Rcpp::NumericVector rt, double v, double a, double w,
                           double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = DBL_EPSILON;
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  int kl;
  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    kl = kl_Nav(t=t/(a*a), eps=eps);
    mult = M_PI * exp(-v*a*w - v*v*t/2) / (a*a);
    out[i] = mult * large_sum_Nav(t=t, a=a, w=w, kl=kl);
  }
  return out;
}





////////////////////////////////////////////////////////////////////////////////
/////////////////////// Combined Small and Large Time //////////////////////////
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
///////////// Densities

// ks = Navarro2009, kl = Navarro2009
// [[Rcpp::export]]
Rcpp::NumericVector f_Nav_Nav(Rcpp::NumericVector rt, double v, double a,
                              double w, double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = DBL_EPSILON;
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  int ks, kl;
  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    ks = ks_Nav(t=t, eps=eps);
    kl = kl_Nav(t=t, eps=eps);
    if(ks < kl) {
      mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_2014(t=t, a=a, w=w, ks=ks);
    }
    else {
      mult = M_PI * exp(-v*a*w - v*v*t/2) / (a*a);
      out[i] = mult * large_sum_Nav(t=t, a=a, w=w, kl=kl);
    }
  }
  return out;
}

// ks = BGK2014, kl = Navarro2009
// [[Rcpp::export]]
Rcpp::NumericVector f_BGK_Nav(Rcpp::NumericVector rt, double v, double a,
                              double w, double t0=0.0, double eps=0.0)
{
  if (eps <= 0.0) { // set eps to default value (doesn't work in definition)
    eps = DBL_EPSILON;
  }

  int n = rt.size(); // get number of values in rt
  Rcpp::NumericVector out(n);
  double t, mult;
  int ks, kl;
  for (int i = 0; i < n; i++) { // iterate through all values in rt and out
    t = rt[i] - t0; // subtract non-decisison time from the response time
    if (t <= 0) {
      out[i] = 0; // if t=0, return 0 instead of NaN
      continue;
    }

    ks = ks_BGK(t=t, w=w, eps=eps);
    kl = kl_Nav(t=t, eps=eps);
    if(ks < kl) {
      mult = a * exp(-v*a*w - v*v*t/2) / sqrt(2*M_PI*t*t*t);
      out[i] = mult * small_sum_2014(t=t, a=a, w=w, ks=ks);
    }
    else {
      mult = M_PI * exp(-v*a*w - v*v*t/2) / (a*a);
      out[i] = mult * large_sum_Nav(t=t, a=a, w=w, kl=kl);
    }
  }
  return out;
}
