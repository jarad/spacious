context("initial.theta.R")


set.seed(1)
y = rnorm(10); vy=var(y)
S = matrix(rnorm(20),10,2)
p.nugget = 0.2
p.range = 0.1

init = initial.theta(y,S,p.nugget,p.range)

test_that("sill computed correctly",{
  expect_equals(init[1], vy)
})

test_that("nugget computed correctly",{
  expect_equals(init[2], p.nugget*vy)
})

test_that("partial.sill computed correctly",{
  expect_equals(init[3], (1-p.nugget)*vy)
  expect_equals(init[3], init[1]-init[2])
})

test_that("range computed correctly",{
  expect_equals(init[4], 0.613532)
})



