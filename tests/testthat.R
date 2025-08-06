# This file is part of the standard setup for testthat. It ensures that
# testthat is run when R CMD check is run on the package.
library(testthat)
library(inv4mhighland)

test_check("inv4mhighland")