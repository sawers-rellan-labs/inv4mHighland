# This file is part of the standard setup for testthat. It ensures that
# testthat is run when R CMD check is run on the package.
library(testthat)
library(inv4mHighland)

test_check("inv4mHighland")