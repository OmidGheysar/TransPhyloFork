# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

probSubtree <- function(tab, rate) {
    .Call('_TransPhylo_probSubtree', PACKAGE = 'TransPhylo', tab, rate)
}

#' Calculate the probability of a phylogenetic tree given a transmission tree
#' @param ctree Combined phylogenetic/transmission tree
#' @param neg Within-host coalescent rate
#' @param w Vector of hosts for which to calculate the probability, or nothing for all
#' @return Probability of phylogeny given transmission tree
#' @export
probPTreeGivenTTree <- function(ctree, neg, w = integer(0)) {
    .Call('_TransPhylo_probPTreeGivenTTree', PACKAGE = 'TransPhylo', ctree, neg, w)
}

coalescent <- function(leaves, nodes, alpha) {
    .Call('_TransPhylo_coalescent', PACKAGE = 'TransPhylo', leaves, nodes, alpha)
}

log_sum_exp <- function(u, v) {
    .Call('_TransPhylo_log_sum_exp', PACKAGE = 'TransPhylo', u, v)
}

log_subtract_exp <- function(u, v) {
    .Call('_TransPhylo_log_subtract_exp', PACKAGE = 'TransPhylo', u, v)
}

log_sum_exp_vec <- function(w) {
    .Call('_TransPhylo_log_sum_exp_vec', PACKAGE = 'TransPhylo', w)
}

roundToPrecision <- function(value, precision) {
    .Call('_TransPhylo_roundToPrecision', PACKAGE = 'TransPhylo', value, precision)
}

wbar <- function(tinf, dateT, rOff, pOff, pi, shGen, scGen, shSam, scSam, delta_t, isTp, time_data, prob_data, dateInitial) {
    .Call('_TransPhylo_wbar', PACKAGE = 'TransPhylo', tinf, dateT, rOff, pOff, pi, shGen, scGen, shSam, scSam, delta_t, isTp, time_data, prob_data, dateInitial)
}

#' Calculates the log-probability of a transmission tree
#' @param ttree Transmission tree
#' @param rOff First parameter of the negative binomial distribution for offspring number
#' @param pOff Second parameter of the negative binomial distribution for offspring number
#' @param pi probability of sampling an infected individual
#' @param shGen Shape parameter of the Gamma probability density function representing the generation time
#' @param scGen Scale parameter of the Gamma probability density function representing the generation time
#' @param shSam Shape parameter of the Gamma probability density function representing the sampling time
#' @param scSam Scale parameter of the Gamma probability density function representing the sampling time
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @param delta_t Grid precision
#' @return Probability of the transmission tree
#' @export
probTTree <- function(ttree, rOff, pOff, pi, shGen, scGen, shSam, scSam, dateT, delta_t = 0.01, isTp = 0L, time_data = as.numeric( c(0.5)), prob_data = as.numeric( c(0.5)), dateInitial = 0) {
    .Call('_TransPhylo_probTTree', PACKAGE = 'TransPhylo', ttree, rOff, pOff, pi, shGen, scGen, shSam, scSam, dateT, delta_t, isTp, time_data, prob_data, dateInitial)
}

