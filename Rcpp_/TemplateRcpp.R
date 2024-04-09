library(Rcpp)

# Define the Rcpp function with default values inside the function body
sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// Helper function to check if a list element is null
bool isNullListElement(List list, std::string name) {
  return Rf_isNull(list[name]);
}

// [[Rcpp::export]]
void processPi(List Pi_values = List::create(Named("time"), Named("probs"))) {
  // Check if time or probs are not provided, and assign default values
  NumericVector time;
  NumericVector probs;
  
  if (isNullListElement(Pi_values, "time")) {
    // Default time values if not provided
    time = NumericVector::create(1, 2, 3);
  } else {
    time = as<NumericVector>(Pi_values["time"]);
  }
  
  if (isNullListElement(Pi_values, "probs")) {
    // Default probability values if not provided
    probs = NumericVector::create(0.5, 0.3, 0.2);
  } else {
    probs = as<NumericVector>(Pi_values["probs"]);
  }

  for (int i = 0; i < time.size(); ++i) {
    Rprintf("Time: %f, Probability: %f\\n", time[i], probs[i]);
  }
}
')


# Call processPi with custom values
custom_Pi_values <- list(
  time = c(4, 5, 6),
  probs = c(0.4, 0.2, 0.1)
)
processPi(custom_Pi_values)
