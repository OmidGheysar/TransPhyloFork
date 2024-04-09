
library(Rcpp)

# Your R code to define Pi_values and call processPi()
Pi_values <- list(
  time = c(1, 2, 3),
  probs = c(0.5, 0.3, 0.2)
)

sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void processPi(List Pi_values) {
  NumericVector time = as<NumericVector>(Pi_values["time"]);
  NumericVector probs = as<NumericVector>(Pi_values["probs"]);

  for (int i = 0; i < time.size(); ++i) {
    Rprintf("Time: %f, Probability: %f\\n", time[i], probs[i]);
  }
}
')

processPi(Pi_values)












