#include <Rcpp.h>
using namespace Rcpp;

/* example change
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <boost/math/tools/roots.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <limits>

struct wstar_functor
{ // Functor also returning 1st derivative.
  wstar_functor(double const& pi, double const& p, double const& r) : pi(pi), p(p), r(r){};

  std::pair<double, double> operator()(double const& x)
  {
    // Return both f(x) and f'(x).
    double temp = pow((1-p)/(1-p*x),r);
    double fx = x - (1-pi)*temp;
    double dx = 1 - (1-pi)*p*r/(1-p*x)*temp;
    return std::make_pair(fx, dx);
  }
private:
  double pi;
  double p;
  double r;
};


double wstar_rootFinder(double pi, double p, double r)
{
  using namespace boost::math::tools;

  const int digits = std::numeric_limits<double>::digits;  // Maximum possible binary digits accuracy for type T.
  int get_digits = static_cast<int>(digits * 0.6);    // Accuracy doubles with each step, so stop when we have
  // just over half the digits correct.
  const boost::uintmax_t maxit = 20;
  boost::uintmax_t it = maxit;
  double result = newton_raphson_iterate(wstar_functor(pi, p, r), 0.5, 0.0, 1.0, get_digits, it);
  return result;
}*/

double wstar_rootFinder(double pi, double p, double r)
{
  Environment pkg = Environment::namespace_env("TransPhylo");
  Function wstar_rootFinderR = pkg["wstar_rootFinderR"];
  RObject res = wstar_rootFinderR(pi,p,r);
  return(as<double>(res));
}

//returns log(exp(u)+exp(v))
// [[Rcpp::export]]
double log_sum_exp(double u, double v)
{
  if (u == R_NegInf) return v;
  if (v == R_NegInf) return u;
  return(std::max(u, v) + log(exp(u - std::max(u, v)) + exp(v - std::max(u, v))));
}

//returns log(exp(u)-exp(v))
// [[Rcpp::export]]
double log_subtract_exp(double u, double v) {
  // if(u <= v) throw(Rcpp::exception("error!! computing the log of a negative number"));
  if(v == R_NegInf)
    return u;
  return u + log1p(-exp(v-u));
}

// [[Rcpp::export]]
double log_sum_exp_vec(NumericVector w)
{
  double total=w[0];
  for (int i = 1; i < w.size(); ++i){
    total=log_sum_exp(total, w[i]);
  }
  return(total);
}

/*
double alphastar(int d, double p, double r, double wstar)
{
  if(std::abs(r-1.0)<1e-6) // Exact solution available
    return log(1-p)+d*log(p)-(d+1)*log(1-p*wstar);
    //return log((1-p))-log(1-p*wstar)+d*log(p/(1-p*wstar));

  int k = d;
  std::vector<double> ltoSum;

  boost::math::negative_binomial_distribution<double> nbinom(r,1-p);
  while(true){

    double dnb = pdf(nbinom,k);
    double lterm = log(dnb) + k*log(wstar);


    ltoSum.push_back(lterm);
    if(lterm <log(1e-8)) break; // Achieved desired accuracy

    k++;
  }

  NumericVector ltoSumR = wrap(ltoSum); // Convert toSum to Rcpp NumericVector
  NumericVector v(k-d+1);
  for(int i=0; i<v.size(); i++) v[i] = i+d;

  return log_sum_exp_vec(lchoose(v, d) + ltoSumR) - d*log(wstar);
}
*/

/*
double alpha(int d, double p, double r,double wbar_tinf)
{
  if(std::abs(r-1.0)<1e-6) // Exact solution available
    return log((1-p))-log_subtract_exp(0.0,log(p)+wbar_tinf)+d*(log(p)-log_subtract_exp(0.0, log(p)+wbar_tinf));
  //    return log(1-p)+d*log(p)-(d+1)*log(1-p*exp(wbar_tinf));

  int k = d;
  std::vector<double> ltoSum;

  while(true){

    double dnb = R::dnbinom(k,r,1-p,1);
    double lterm = dnb + k*wbar_tinf;

    ltoSum.push_back(lterm);
    if(log(Rf_choose(k,d))+lterm < log(1e-8)) break; // Achieved desired accuracy

    k++;

    if(k>1e6) throw(Rcpp::exception("too many iterations, giving up!"));

  }
  NumericVector ltoSumR = wrap(ltoSum); // Convert toSum to Rcpp NumericVector
  NumericVector v(k-d+1);
  for(int i=0; i<v.size(); i++) v[i] = i+d;

  return (log_sum_exp_vec( lchoose(v, d)+ltoSumR ) - d*wbar_tinf);
}
*/

// alpha() computes equation (10) in TransPhylo paper
double alpha(int d, double p, double r,double wbar_tinf)
{
  if(std::abs(r-1.0)<1e-6) // Exact solution available
    return log((1-p))-log_subtract_exp(0.0,log(p)+wbar_tinf)+d*(log(p)-log_subtract_exp(0.0, log(p)+wbar_tinf));
  double toret=0;
  double term=R::dnbinom(d,r,1-p,1);
  double eterm=exp(term);
  for (int k=d;k<1e4;k++) {
    toret+=eterm;
    term+=log(p*(k+r)/(k-d+1))+wbar_tinf;
    eterm=exp(term);
    if (eterm<0.001*toret) break;
  }
  return(log(toret));
}



// here is my code to Tp ---------------------------------------------------------------------------------
// Define the step function
// [[Rcpp::export]]
NumericVector cutPi(NumericVector pi, double initialDate, double infectionDate, double delta, int n) {
  int initialIndex = std::round((infectionDate - initialDate) / delta);
  
  // If initialIndex is negative or exceeds the size of pi, return an empty vector
  if (initialIndex < 0 || initialIndex >= pi.size())
    return NumericVector();
  
  int endIndex = initialIndex + n - 1;
  
  // If endIndex exceeds the size of pi, set it to the last index of pi
  if (endIndex >= pi.size())
    endIndex = pi.size() - 1;
  
  NumericVector cut(pi.begin() + initialIndex, pi.begin() + endIndex + 1);
  return cut;
}


// [[Rcpp::export]]
double roundToPrecision(double value, double precision) {
  double roundedValue = round(value / precision) * precision;
  return roundedValue;
}


// [[Rcpp::export]]
double step_function_3(double x, double threshold1, double threshold2, double value1, double value2, double value3) {
  double y;
  if (x < threshold1) {
    y = value1;
  } else if (x < threshold2) {
    y = value2;
  } else {
    y = value3;
  }
  return y;
}

// Define the gamma function
// [[Rcpp::export]]
double Dgamma(double calendar_time, double shape, double scale, double tinf) {
  return R::dgamma( calendar_time-tinf, shape, scale, false);
}

// Define the integral function using Rcpp
// [[Rcpp::export]]
double integral_function(NumericVector calendar_time, double tinf, double shape, double scale, double value1, double value2, double value3, double threshold1, double threshold2) {
  int n = calendar_time.size();
  NumericVector y(n);

  for (int i = 1; i < n; i++) {
    y[i] = Dgamma(calendar_time[i], shape, scale, tinf) * step_function_3(calendar_time[i], threshold1, threshold2, value1, value2, value3);
  }

  double sum = 0.0;
  for (int i = 1; i < n; i++) {
    double dx = calendar_time[i] - calendar_time[i-1];
    double area = (y[i] + y[i-1]) * dx / 2.0;
    sum += area;
  }

  return sum;
}

// Define the integral function using Rcpp
// [[Rcpp::export]]
NumericVector integral_values_cpp(NumericVector calendar_time, NumericVector tinf, double shape, double scale, double value1, double value2, double value3, double threshold1, double threshold2) {
  int n_infections = tinf.size();
  NumericVector integral_values(n_infections);

  for (int i = 0; i < n_infections; i++) {
    double t = tinf[i];
    integral_values[i] = integral_function(calendar_time, t, shape, scale, value1, value2, value3, threshold1, threshold2);
  }

  return integral_values;
}

// end of my code to Tp ----------------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector wbar(double tinf, double dateT, double rOff, double pOff, double pi, double shGen, double scGen, double shSam, double scSam, double delta_t, int isTp, NumericVector Pi, double dateInitial)
{
  int n = std::round((dateT-tinf)/delta_t);
  NumericVector grid(n);
  NumericVector pi2(n);
  for(int i=0; i<n; ++i) // use the left point of each subinterval
    grid[i] = dateT-n*delta_t+i*delta_t;
  // NumericVector pi2;
  if (isTp == 1){
    // pi = 0;
    pi2 = pi*pgamma(dateT-grid, shSam, scSam);
    // Rprintf("pi2[10] for Tp: %f\n", pi2[10]);
    // Rprintf("pi2[20] for Tp: %f\n", pi2[20]);
  }
// Here is my Code to Tp ----------------------------------------------------------------------------------
  // if(isTp==0){
  //   // NumericVector calendar(n);
  //   // for(int i=0; i<n; ++i)
  //   //   calendar [i] = tinf + i*delta_t;
  //   // double value1 = 0.2;
  //   // double value2 = 0.8;
  //   // double value3 = 0.3;
  //   // double threshold1 = 2006 ;
  //   // double threshold2 = 2007;
  //   // pi2 =  integral_values_cpp(calendar, calendar, shSam, scSam, value1, value2, value3, threshold1, threshold2);
  //   // pi2 = cutPi(Pi, dateInitial, tinf, delta_t, n);
  //   // pi2 = pi;
  //   for (int i=0; i<n; ++i)
  //     pi2[i] = grid[i];
  //   if (n >= 0 && n <= 100) {
  //     pi2 = rep(0.9, n);
  //   } else if (n > 100 && n <= 200) {
  //     pi2 = rep(0.1, n);
  //   } else if (n > 200) {
  //     pi2 = rep(0.9, n);
  //   } else {
  //     // Handle other cases if needed
  //     Rcpp::stop("Invalid value of n.");
  //   }
    
    if (isTp == 0) {
      double value1 = Pi[0];
      double value2 = Pi[1];
      double value3 = Pi[2];
      double threshold1 =  2005.2;
      double threshold2 =  2006;
      for (int i = 0; i < n; ++i) {
        double g_value = grid[i]; // Replace with the actual calculation of g[i]
        // Rprintf("pi2[120]: %f\n", grid[i]);
        
        if (g_value <= threshold1) {
          pi2[i] = value1;
          Rprintf("pi2: %f\n", value1);
        } else if (g_value > threshold1 && g_value <= threshold2) {
          pi2[i] = value2;
          Rprintf("pi2: %f\n", value2);
        } else {
          pi2[i] = value3;
          Rprintf("pi2: %f\n", value3);
        }
      }
    }
    
    // for(int i=0; i<n; ++i)
    //   pi2[i] = grid[i];
    
    // Rprintf("dateIntial: %f\n", dateInitial);
    // Rprintf("tinf: %f\n", tinf);
    // Rprintf("pi2[10]: %f\n", pi2[10]);
    // Rprintf("pi2[120]: %f\n", pi2[120]);
    // Rprintf("shape: %f\n", shSam);
    // Rprintf("scale: %f\n", scSam);
    // Rprintf("calendar1: %f\n", calendar[1]);
    // Rprintf("calendar2: %f\n", calendar[2]);
    // Rprintf("calendar100: %f\n", calendar[100]);
    // Rprintf("delta: %f\n", delta_t);
    // Rprintf("Value: %d\n", n);
    // Rprintf("isTp: %d\n", isTp);
    // Rprintf("dateT-tinf: %f\n",  std::round((dateT-tinf)/delta_t));
  
  if(isTp==2){
    double value1 = Pi[0];
    double value2 = Pi[1];
    double value3 = Pi[2];
    double threshold1 =  Pi[3];
    double threshold2 =  Pi[4];
    for (int i = 0; i < n; ++i) {
      double g_value = grid[i]; // Replace with the actual calculation of g[i]
      // Rprintf("pi2[120]: %f\n", grid[i]);
      
      if (g_value <= threshold1) {
        pi2[i] = value1;
        // Rprintf("pi2: %f\n", value1);
      } else if (g_value > threshold1 && g_value <= threshold2) {
        pi2[i] = value2;
        // Rprintf("pi2: %f\n", value2);
      } else {
        pi2[i] = value3;
        // Rprintf("pi2: %f\n", value3);
      }
    }
    
  }
  
  if(isTp==3){
    double value1 = 0.01;
    double value2 = 0.01;
    double value3 = 0.9;
    double threshold1 =  2000.00;
    double threshold2 =  2000.50;
    for (int i = 0; i < n; ++i) {
      double g_value = grid[i]; // Replace with the actual calculation of g[i]
      // Rprintf("pi2[120]: %f\n", grid[i]);
      
      if (g_value <= threshold1) {
        pi2[i] = value1;
        // Rprintf("pi2: %f\n", value1);
      } else if (g_value > threshold1 && g_value <= threshold2) {
        pi2[i] = value2;
        // Rprintf("pi2: %f\n", value2);
      } else {
        pi2[i] = value3;
        // Rprintf("pi2: %f\n", value3);
      }
    }
    
  }
// end of my code to Tp -----------------------------------------------------------------------------------

  NumericVector F = 1-pgamma(dateT-grid, shGen, scGen);

  NumericVector w(n), out(n);

  IntegerVector seq = seq_len(n);
  NumericVector gam = dgamma(as<NumericVector>(seq)*delta_t,shGen,scGen);
  double sumPrev = 0.5 * gam[0];
  out[n-1]=std::min(1.0,F[n-1]+sumPrev*delta_t);
  for(int i=n-1; i>0; --i){
    w[i] = (1-pi2[i]) * pow((1-pOff)/(1-pOff*out[i]), rOff);

    sumPrev = 0.0;
    for(int j=0; j<n-i; ++j)
      sumPrev += gam[j]*w[i+j];
    sumPrev += 0.5 * gam[n-i];
    out[i-1] = std::min(1.0,F[i-1] + sumPrev*delta_t);
  }
  return log(out);
}

//' Calculates the log-probability of a transmission tree
//' @param ttree Transmission tree
//' @param rOff First parameter of the negative binomial distribution for offspring number
//' @param pOff Second parameter of the negative binomial distribution for offspring number
//' @param pi probability of sampling an infected individual
//' @param shGen Shape parameter of the Gamma probability density function representing the generation time
//' @param scGen Scale parameter of the Gamma probability density function representing the generation time
//' @param shSam Shape parameter of the Gamma probability density function representing the sampling time
//' @param scSam Scale parameter of the Gamma probability density function representing the sampling time
//' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
//' @param delta_t Grid precision
//' @return Probability of the transmission tree
//' @export
// [[Rcpp::export]]
double probTTree(NumericMatrix ttree, double rOff, double pOff, double pi,
                 double shGen, double scGen, double shSam, double scSam,
                 double dateT, double delta_t=0.01, int isTp = 0, NumericVector Pi = NumericVector::create(0.5), double dateInitial = 0){
  // Rprintf("Pi for Tp: %f\n", Pi[0]);
  int numCases = ttree.nrow();

  if(shGen*scGen<0.001) throw(Rcpp::exception("error!! mean of gamma is too small."));

  if(dateT == INFINITY){ // finished outbreak
    double wstar = wstar_rootFinder(pi, pOff, rOff);

    NumericVector lsstatus = ifelse(is_na(ttree(_,1)), log(1-pi), log(pi)+dgamma(ttree(_,1)-ttree(_,0),shSam,scSam,1));

    std::map<int, std::vector<int> > infMap; // Map from infector to infected
    std::vector<std::vector<int> > progeny(numCases);
    for(int i=0; i<numCases; ++i){
      if(ttree(i,2) == 0) continue; // Found root node i

      progeny[ttree(i,2)-1].push_back(i); // C++ index starts from 0
      infMap[ttree(i,2)-1] = progeny[ttree(i,2)-1];
    }
    double accum = 0.0;
    for(int i=0; i<numCases; ++i){
      accum += alpha(progeny[i].size(), pOff, rOff, log(wstar));

      for(unsigned int j=0; j<progeny[i].size(); ++j){
        accum += R::dgamma(ttree(progeny[i][j],0) - ttree(i,0), shGen, scGen, 1);
      }
    }
    return sum(lsstatus) + accum;
  }
  else{
    // Ongoing outbreak -- observation ends at finite dateT
    NumericVector lprobSam = log(pi)+pgamma(dateT-ttree(_,0),shSam,scSam,1,1);
    for(int i=0; i<lprobSam.size(); ++i){
      lprobSam[i] = log_subtract_exp(0.0,lprobSam[i]);
    }
    NumericVector lsstatus = ifelse(is_na(ttree(_,1)), lprobSam, log(pi)+dgamma(ttree(_,1)-ttree(_,0),shSam,scSam,1));
    std::map<int, std::vector<int> > infMap; // Map from infector to infected
    std::vector<std::vector<int> > progeny(numCases);
    for(int i=0; i<numCases; ++i){
      if(ttree(i,2) == 0) continue; // Found root node i

      progeny[ttree(i,2)-1].push_back(i); // C++ index starts from 0
      infMap[ttree(i,2)-1] = progeny[ttree(i,2)-1];
    }
    double accum = 0.0;
    double tinfmin = min(ttree(_,0));
    tinfmin = roundToPrecision(tinfmin, delta_t);
    NumericVector wbar0 = wbar(tinfmin, dateT, rOff, pOff, pi, shGen, scGen, shSam, scSam, delta_t, isTp, Pi, dateInitial);

    double gridStart = dateT-wbar0.size()*delta_t;

    for(int i=0; i<numCases; ++i){

      accum += alpha(progeny[i].size(), pOff, rOff,wbar0[std::min(wbar0.size()-1.0,std::round((ttree(i,0) - gridStart)/delta_t))]);
      for(unsigned int j=0; j<progeny[i].size(); ++j){
        accum += (R::dgamma(ttree(progeny[i][j],0)-ttree(i,0), shGen, scGen, 1) - R::pgamma(dateT-ttree(i,0), shGen, scGen, 1, 1));
      }
    }
    return sum(lsstatus) + accum;
  }
}
