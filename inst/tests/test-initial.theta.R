context("initial.theta.R")


set.seed(1)
y = rnorm(10); vy=var(y)
S = matrix(rnorm(20),10,2)
p.nugget = 0.2
p.range = 0.1

init = as.numeric(initial.theta(y,S,p.nugget,p.range))


test_that("nugget computed correctly",{
  expect_equal(init[1], p.nugget*vy)
})

test_that("partial.sill computed correctly",{
  expect_equal(init[2], (1-p.nugget)*vy)
})

test_that("range computed correctly",{
  expect_equal(init[3], 0.613532)
})

