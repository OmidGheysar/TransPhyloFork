#include <Rcpp.h>
using namespace Rcpp;

/*
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

// [[Rcpp::export]]
double roundToPrecision(double value, double precision) {
  double roundedValue = round(value / precision) * precision;
  return roundedValue;
}

// Calculate the integral based on the conditions
// [[Rcpp::export]]
NumericVector calculate_integral(NumericVector tinput, double change_point, double alpha, double beta, double p1, double p2, double Tend) {
  int n = tinput.size();
  NumericVector integral_values(n);
  
  for (int i = 0; i < n; i++) {
    double timeInf = tinput[i];
    if (timeInf <= change_point) {
      double cdf_change_point_minus_timeInf = R::pgamma(change_point - timeInf, alpha, beta, 1, 0);
      double integral_1 = p1 * cdf_change_point_minus_timeInf;
      double integral_2 = p2 * (R::pgamma(Tend - timeInf, alpha, beta, 1, 0) - cdf_change_point_minus_timeInf);
      integral_values[i] = integral_1 + integral_2;
    } 
    else {
      double integral_1 = p2 * (R::pgamma(Tend - timeInf, alpha, beta, 1, 0) - R::pgamma(0, alpha, beta, 1, 0));
      integral_values[i] = integral_1;
    }
  }
  return integral_values;
}
// [[Rcpp::export]]
NumericVector wbar(double tinf, double dateT, double rOff, double pOff, double pi, double shGen, double scGen, double shSam, double scSam, double delta_t, int isTp, NumericVector time_data, NumericVector prob_data)
{
  
  int n = std::round((dateT-tinf)/delta_t);
  NumericVector grid(n);
  NumericVector pi2(n);
  // this is for tuncation intergral equation 8
  NumericVector piInt(n); 
  for(int i=0; i<n; ++i) // use the left point of each subinterval
    grid[i] = dateT-n*delta_t+i*delta_t;
  // this is the integral of equation 8
  piInt = pgamma(dateT-grid, shSam, scSam);
  // NumericVector pi2;
  if (isTp == 1){
    // pi = 0;
    pi2 = pi*pgamma(dateT-grid, shSam, scSam);
    // Rprintf("pi2[10] for Tp: %f\n", pi2[10]);
    // Rprintf("pi2[20] for Tp: %f\n", pi2[20]);
  }
  // Here is my Code to Tp ----------------------------------------------------------------------------------
  
  if(isTp==2){
    // Rprintf("tinf: %f\n", tinf);
    NumericVector tinput= grid;
    double change_point= time_data[0];
    double alph= shSam;
    double beta= scSam;
    double p1= prob_data[0];
    double p2= prob_data[1];
    double Tend= dateT;
    pi2 = calculate_integral(grid,change_point, alph, beta, p1, p2, Tend);
    // Rprintf("tinf0: %f\n", pi2[0]);
    // Rprintf("tinf100: %f\n", pi2[100]);
    // Rprintf("tinf600: %f\n", pi2[600]);
    // for (int i = 0; i < pi2.size(); i++) {
    //   Rprintf("pi2[%d]: %f\n", i, pi2[i]);
    // }
  }
  
  if(isTp==3){
    // Find the index where time_data matches the first element of grid
    int matchingIndex = -1;
    for(int i=0; i<time_data.size(); ++i) {
      if(time_data[i] == grid[0]) {
        matchingIndex = i;
        break;
      }
    }
    
    // Check if a matching index was found
    if(matchingIndex != -1) {
      // Assign values from prob_data to pi2 starting from the found index
      for(int i = matchingIndex; i < time_data.size(); ++i) {
        int pi2Index = i - matchingIndex;
        if(pi2Index < n) { // Ensure we don't go out of bounds for pi2
          pi2[pi2Index] = prob_data[i];
          // Rprintf("pi2[pi2Index]: %f\n", pi2[pi2Index]);
        } else {
          // If there's no corresponding index in pi2 for the prob_data, break the loop
          break;
        }
      }
      // Here I will add truncation the equation in the equation 8 in the Transphylo paper 
      for (int i = 0; i < pi2.size(); ++i) {
        pi2[i] = pi2[i] * piInt[i];
      }
    } else {
      Rprintf("there is no match between these two sets");
      // Handle the case where no matching index is found
      // For example, you might want to set pi2 to some default values or throw an error
    }
    
    // ... rest of your code ...
    
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


// [[Rcpp::export]]
NumericVector Ptime(NumericVector ttree_col0, NumericVector time_data, NumericVector prob_data) {
  int n = ttree_col0.size();
  NumericVector pTime(n);
  
  // Ensure time_data is sorted (we assume it is always sorted as per the problem statement)
  // For each element in ttree_col0
  for (int i = 0; i < n; ++i) {
    double t = ttree_col0[i];
    
    // Use binary search to find the closest value
    auto it = std::lower_bound(time_data.begin(), time_data.end(), t);
    
    int closest_index;
    if (it == time_data.end()) {
      closest_index = time_data.size() - 1; // If t is greater than any element in time_data
    } else if (it == time_data.begin()) {
      closest_index = 0; // If t is less than any element in time_data
    } else {
      int idx = it - time_data.begin();
      // Check the closest one
      if (std::abs(time_data[idx] - t) < std::abs(time_data[idx - 1] - t)) {
        closest_index = idx;
      } else {
        closest_index = idx - 1;
      }
    }
    
    // Pick the value from prob_data at the closest index
    pTime[i] = prob_data[closest_index];
  }
  
  return pTime;
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
                  double dateT, double delta_t=0.01, int isTp = 1, NumericVector time_data = NumericVector::create(0.5), NumericVector prob_data = NumericVector::create(0.5))
 {
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
     // NumericVector lprobSam = log(pi)+pgamma(dateT-ttree(_,0),shSam,scSam,1,1);
     if(isTp == 2){
       // Rprintf("inside isTp =3");
       NumericVector turncWS = pgamma(dateT-ttree(_,0),shSam,scSam,1,1);
       // NumericVector PtimeOut = Ptime(ttree(_,0),time_data, prob_data);
       // NumericVector lprobSam = log(PtimeOut);
       NumericVector tinput= ttree(_,0);
       // NumericVector tinput = pgamma(dateT-ttree(_,0),shSam,scSam,1,1);
       double change_point= time_data[0];
       double alph= shSam;
       double beta= scSam;
       double p1= prob_data[0];
       double p2= prob_data[1];
       double Tend= dateT; 
       
       NumericVector Ptime = calculate_integral(tinput, change_point, alph, beta, p1, p2, Tend);
       NumericVector lprobSam = log(Ptime);
       for(int i=0; i<lprobSam.size(); ++i){
         lprobSam[i] = log_subtract_exp(0.0,lprobSam[i]);
       }
       // NumericVector lsstatus = ifelse(is_na(ttree(_,1)), lprobSam, log(pi)+dgamma(ttree(_,1)-ttree(_,0),shSam,scSam,1));
       NumericVector lsstatus = ifelse(is_na(ttree(_,1)), lprobSam, log(Ptime)+dgamma(ttree(_,1)-ttree(_,0),shSam,scSam,1) - turncWS);
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
       NumericVector wbar0 = wbar(tinfmin, dateT, rOff, pOff, pi, shGen, scGen, shSam, scSam, delta_t, isTp, time_data, prob_data);
       
       double gridStart = dateT-wbar0.size()*delta_t;
       
       for(int i=0; i<numCases; ++i){
         
         accum += alpha(progeny[i].size(), pOff, rOff,wbar0[std::min(wbar0.size()-1.0,std::round((ttree(i,0) - gridStart)/delta_t))]);
         for(unsigned int j=0; j<progeny[i].size(); ++j){
           // accum += R::dgamma(ttree(progeny[i][j],0)-ttree(i,0), shGen, scGen, 1);
           accum += (R::dgamma(ttree(progeny[i][j],0)-ttree(i,0), shGen, scGen, 1) - R::pgamma(dateT-ttree(i,0), shGen, scGen, 1, 1));
         }
       }
       return sum(lsstatus) + accum;
     }else {
       // Rprintf("inside isTp not 3");
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
       NumericVector wbar0 = wbar(tinfmin, dateT, rOff, pOff, pi, shGen, scGen, shSam, scSam, delta_t, isTp, time_data, prob_data);
       
       double gridStart = dateT-wbar0.size()*delta_t;
       
       for(int i=0; i<numCases; ++i){
         
         accum += alpha(progeny[i].size(), pOff, rOff,wbar0[std::min(wbar0.size()-1.0,std::round((ttree(i,0) - gridStart)/delta_t))]);
         for(unsigned int j=0; j<progeny[i].size(); ++j){
           // accum += R::dgamma(ttree(progeny[i][j],0)-ttree(i,0), shGen, scGen, 1);
           accum += (R::dgamma(ttree(progeny[i][j],0)-ttree(i,0), shGen, scGen, 1)- R::pgamma(dateT-ttree(i,0), shGen, scGen, 1, 1));
         }
       }
       return sum(lsstatus) + accum;
     }
   }
 }