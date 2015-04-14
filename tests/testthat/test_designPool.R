context("designPool")

test_that("designPool:string arguments", {
	expect_error(designPool(type = "complete"))
	expect_error(designPool(type = com))
})

test_that("designPool:time limitation", {
	expect_error(designPool(type = "com", time = 11, traj = "linear", budget = 100000, unitcost = 20))
	expect_error(designPool(type = "both", time = 11, traj = "linear", budget = 100000, unitcost = 20))
	expect_error(designPool(type = "com", time = 1, traj = "linear", budget = 100000, unitcost = 20))
	expect_error(designPool(type = "both", time = 3, traj = "quadratic", budget = 100000, unitcost = 20))
})

test_that("designPool:sample size", {
	expect_error(designPool(type = "com", time = 3, traj = "linear", budget = 50, unitcost = 20)) # N >= 1
	expect_error(designPool(type = "miss", time = 3, traj = "linear", budget = 30, unitcost = 20)) 
	
	expect_error(designPool(type = "miss", time = 3, traj = "linear", budget = 30, unitcost = 20, nfix = 100)) 
	expect_error(designPool(type = "miss", time = 3, traj = "linear", unitcost = 20, nfix = 100)) 
	expect_error(designPool(type = "miss", time = 3, traj = "linear", budget = 30, nfix = 100)) 
	
	expect_error(designPool(type = "com", time = 3, traj = "linear", budget = 30, unitcost = 20, nfix = 100)) 
	expect_error(designPool(type = "com", time = 3, traj = "linear", unitcost = 20, nfix = 100)) 
	expect_error(designPool(type = "com", time = 3, traj = "linear", budget = 30, nfix = 100)) 
	
	expect_error(designPool(type = "both", time = 3, traj = "linear", budget = 30, unitcost = 20, nfix = 100)) 
	expect_error(designPool(type = "both", time = 3, traj = "linear", unitcost = 20, nfix = 100)) 
	expect_error(designPool(type = "both", time = 3, traj = "linear", budget = 30, nfix = 100)) 
	
	expect_error(designPool(type = "com", time = 3, traj = "linear", nfix = 0.5))  # nfix >= 1
	expect_error(designPool(type = "miss", time = 3, traj = "linear", nfix = 0.5))  # nfix >= 1
	expect_error(designPool(type = "both", time = 3, traj = "linear", nfix = 0.5))  # nfix >= 1
})

test_that("designPool:attrition", {
	expect_error(designPool(type = "both", time = 5, traj = "quadratic", budget = 100000, unitcost = 20, attrition = 0.5))
	expect_error(designPool(type = "both", time = 5, traj = "quadratic", budget = 100000, unitcost = 20, attrition = 0))
	expect_error(designPool(type = "both", time = 5, traj = "quadratic", budget = 100000, unitcost = 20, attrition = 1))
	expect_error(designPool(type = "both", time = 5, traj = "quadratic", budget = 100000, unitcost = 20, attrition = rep(0.5, 6)))
	
	
	expect_error(designPool(type = "both", time = 5, traj = "quadratic", budget = 100000, unitcost = 20, attrition = rep(1, 5)))
	expect_error(designPool(type = "both", time = 5, traj = "quadratic", budget = 100000, unitcost = 20, attrition = rep(NA, 5)))
	expect_error(designPool(type = "both", time = 5, traj = "quadratic", budget = 100000, unitcost = 20, attrition = c(rep(0, 3),0.5,1.5)))
})
