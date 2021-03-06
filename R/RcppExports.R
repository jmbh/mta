# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

distmat <- function(id, x, y, n) {
    .Call('mta_distmat', PACKAGE = 'mta', id, x, y, n)
}

rootChoose <- function(n, k, root) {
    .Call('mta_rootChoose', PACKAGE = 'mta', n, k, root)
}

rootChooseLookup <- function(n, k, lookup) {
    .Call('mta_rootChooseLookup', PACKAGE = 'mta', n, k, lookup)
}

lookup <- function(n, root) {
    .Call('mta_lookup', PACKAGE = 'mta', n, root)
}

rootCombLookup <- function(ns, lookup) {
    .Call('mta_rootCombLookup', PACKAGE = 'mta', ns, lookup)
}

stabExp <- function(ns, lookup) {
    .Call('mta_stabExp', PACKAGE = 'mta', ns, lookup)
}

f_rescale_c <- function(x, y, npts) {
    .Call('mta_f_rescale_c', PACKAGE = 'mta', x, y, npts)
}

spatialRescale_c <- function(trs, npts) {
    .Call('mta_spatialRescale_c', PACKAGE = 'mta', trs, npts)
}

