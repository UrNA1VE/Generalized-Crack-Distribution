
IG_random = function(n, lambda, theta) {
  u = runif(n, 0, 1)
  z = rnorm(n, 0, 1)
  v = lambda * theta + (theta / 2) * (z^2 - sqrt(z^4 + 4 * lambda * z^2)) 
  compare = (lambda * theta) / (lambda * theta + v)
  index = u < compare
  another_choice = (lambda^2) * (theta^2) / v
  IG = c(v[index], another_choice[!index])
  return(sample(IG))
}


LB_random = function(n, lambda, theta) {
  z2 = rnorm(n, 0, 1)
  IG = IG_random(n, lambda, theta)
  LB = IG + theta * z2^2
  return(sample(LB))
}


Crack_random = function(n, lambda, theta, p) {
  IG = IG_random(n, lambda, theta)
  LB = LB_random(n, lambda, theta)
  u = runif(n, 0, 1)
  index = u < p
  Crack = c(IG[index], LB[!index])
  return(sample(Crack))
}

# The density function of LB2
f_LB2 = function(x, lambda, theta) {
  first_part = sqrt(x / theta) / (theta * sqrt(2*pi) * (lambda + 1))
  second_part = (-1/2) * (sqrt(x/theta) - lambda * sqrt(theta/x))^2
  return(first_part*(exp(second_part)))
}

# Step 1: C = exp(lambda), g(y) = Gamma(shape = 3/2, beta = 2*theta)
# Step 2: generate u ~ uniform(0, 1)
# Step 3: generate y ~ g(y)
# Step 4: if u <= fLB2(y)/(C * g(y)), x = y else go to step 2
# I avoid for loop just to improve the efficiency
LB2_random = function(n, lambda, theta) {
  LB2 = c()
  while(length(LB2) < n) {
    u = runif(n, 0, 1)
    gamma = rgamma(n, shape = 3/2, scale = 2*theta)
    compare = f_LB2(gamma, lambda, theta) / (exp(lambda) * dgamma(gamma, shape = 3/2, scale = 2*theta))
    index = u <= compare
    LB2 = c(LB2, gamma[index])
  }
  return(sample(LB2)[1:n])
}


GCR_random = function(n, lambda, theta, p, q) {
  r = 1 - p - q
  u = runif(n, 0, 1)
  Crack = Crack_random(n, lambda, theta, p / (1 - r))
  LB2 = LB2_random(n, lambda, theta)
  index = u <= (1 - r)
  GCR = c(Crack[index], LB2[!index])
  return(sample(GCR))
}



# IG_random(100, 1, 1)
# LB_random(100, 1, 1)
# Crack_random(100, 1, 1, 0.5)
# LB2_random(100, 1, 1)
# GCR_random(100, 1, 1, 0.3, 0.3)

# density function for GCR distribution
f_GCR = function(x, lambda, theta, p, q) {
  r = 1 - p - q
  first_part = sqrt(theta/(2*pi))
  second_part = (p*lambda*(x^(-3/2)) + (q/theta)*x^(-1/2) + (r/((lambda + 1)*theta^(2)))*x^(1/2))
  third_part = (-1/2) * (sqrt(x/theta) - lambda*(sqrt(theta/x)))^(2)
  return(first_part * second_part * exp(third_part))
}


# origin_bins = function(x, k) {
#   y1 = min(x)
#   ylast = max(x)
#   bins = c(0, y1)
#   delta = (ylast - y1) / (k - 2)
#   if (k > 3) {
#     for (i in 1: (k - 3)) {
#       bins = c(bins, c(y1 + i * delta))
#     }
#   }
#   bins = c(bins, c(ylast, Inf))
#   return(bins)
#   
# }
# 
# BS_bins = function(k, lambda, theta) {
#   BS_percentile = function(p, lambda, theta) {
#     alpha = 1 / sqrt(lambda)
#     beta = theta * lambda
#     return((beta / 4) * ((alpha * qnorm(p) + sqrt((alpha * qnorm(p))** 2 +4)) ** 2))
#   }
#   bins = (1 : (k - 1)) / k
#   bins = BS_percentile(bins, lambda, theta)
#   bins = c(0, bins, Inf)
#   return(bins)
# }


# go through chi2 goodness of fit by calling GCR_test()
# step 0 fix parameters
# GCR_test = function(n = 127, lambda = 2, theta = 0.5, p = 1/3, q = 1/3, bins_method = "BS") {
# 
#   # step 1 generate random numbers and get k by Sturges' Rule
#   k = round(log2(n) + 1)
#   if (k <= 2) {
#     stop("Not enough number")
#   }
#   GCR_numbers = GCR_random(n, lambda, theta, p, q)
# 
#   # step 2 get class boundaries, bins are the boundaries like c(0, y1...., yk-1, Inf)
#   #origin means using ymax, ymin and delta to generate bins.
#   if (bins_method == "origin") {
#     bins = origin_bins(GCR_numbers, k)
#   }
#   #BS means using BS distribution to generate bins.
#   else if (bins_method == "BS") {
#     bins = BS_bins(k, lambda, theta)
#   }
# 
#   # step 3 get expected values E like c(n*P1, n*p2 ... n*pk)
#   f = function(x) {
#     return(f_GCR(x, lambda, theta, p, q))
#   }
#   E = c()
#   for (i in 1: k) {
#     Pi = integrate(f, bins[i], bins[i + 1])$value
#     E = c(E, n * Pi)
#   }
# 
#   # step 3.5 check and merge bins whose expected value is less than 5
# 
#   old_bins = bins
#   old_E = E
# 
#   i = 1
#   while (i <= length(E)){
#     if ((i < length(E)) && (E[i] < 5)) {
#       bins = bins[-c(i + 1)]
#       E[i + 1] = E[i + 1] + E[i]
#       E = E[-c(i)]
#       k = k - 1
# 
#     }
#     else if ((i == length(E)) && (E[i] < 5)) {
#       bins = bins[-c(i)]
#       E[i - 1] = E[i - 1] + E[i]
#       E = E[-c(i)]
#       k = k - 1
#       i = i - 1
# 
#     }
#     else if (E[i] >= 5) {
#       i = i + 1
# 
#     }
#   }
# 
# 
#   # step 4 count observations O like c(O1, O2... Ok)
#   O = rep(0, k)
#   for (i in 1: n) {
#     x = GCR_numbers[i]
#     for (j in 1 : k) {
#       if (bins[j + 1] > x) {
#         break
#       }
#     }
#     O[j] = O[j] +1
#   }
# 
#   # step 5 compute chi2 statistics
#   chi2 = sum((O - E)^(2) / E)
#   p_value = pchisq(chi2, k - 1, lower.tail = FALSE)
# 
#   res = list(GCR_numbers, bins, E, O, p_value, old_bins, old_E)
#   names(res) = c("GCR_numbers", "bins", "expected", "observed", "p_value", "old_bins", "old_E")
# 
#   return(res)
# }


# set.seed(2021)
# res = GCR_test()
# GCR_test() outputs a list which contains GCR random numbers, bins, expected values, observed values and p-values
# for the chi-square test. By calling res$GCR_numbers, res$bins, res$expected, res$observed, res$p_values to check
# its corresponding value. Calling res$old_bins, res$old_E to check the original bins and expected values(before merging).
# res$old_E
# res$GCR_numbers
# res$bins
# res$expected
# res$observed
# res$p_value
# res$old_E
# res$old_bins

# res = GCR_test(n = 1000, lambda = 2.5, theta = 1.5, p = 1/2, q = 1/5)
# res = GCR_test(n = 200, lambda = 2.5, theta = 0.5, p = 1/3, q = 1/5)

# I set using BS distribution(alpha = 1 / sqrt(lambda) and beta = theta * lambda) to generate bins 
# as default but still preserve the original bin method. If you want to use the ymax, ymin and delta
# way to generate bins, simply call bins_method = "origin".

# res = GCR_test(n = 1000, lambda = 2.5, theta = 1.5, p = 1/2, q = 1/5, bins_method="origin")
# 
MLE_l = function(lambda, theta, p, q, x) {
  n = length(x)
  r = 1 - q - p
  first = (n / 2) * (log(theta) - log(2*pi))
  second = p * lambda * (x**(-3/2)) + (q/theta) * (x**(-1/2)) + (r/((lambda +1)*(theta**2)))*(x**(1/2))
  second = sum(log(second))
  third = n*lambda - (1/2)*sum(x/theta + (lambda**2)*theta/x)
  return(first + second + third)
}

# MLE_f = function(lambda, theta, p, q, x) {
#   n = length(x)
#   second = (p*(theta**2) - r*(x**2)/(lambda + 1)**2) / (p*lambda*(theta**2) + q*theta*x + r*(x**2)/(lambda+1))
#   return(n - sum(second) - lambda*theta*sum(1/x))
# }
# 
# MLE_g = function(lambda, theta, p, q, x) {
#   r = 1 - q - p
#   n = length(x)
#   first = n/(2*theta) + (2*n)/theta
#   second = (2*p*lambda*theta + q*x) / (p*lambda*(theta**2) + q*theta*x + r*(x**2)/(lambda +1))
#   third = 1/(2*(theta**2)) * sum(x) - (lambda**2)/2 * sum(1/x)
#   return(first - sum(second) +third)
# }
# 
# MLE_h = function(lambda, theta, p, q, x) {
#   return(MLE_f(lambda, theta, p, q, x)**2 + MLE_g(lambda, theta, p, q, x)**2)
# }




# Round one and round two
# The MLE_test_onetime function would print the time it needed to perform the computation of lambda and theta, and output a result from optim function.
MLE_test_onetime = function(p, q, x) {

  l = function(paras) {
    lambda = paras[1]
    theta = paras[2]
    return(-1*MLE_l(lambda, theta, p = p, q = q, x = x))
  }

  # time1 = Sys.time()
  res = optim(par = c(1, 1), fn = l, method = "L-BFGS-B", lower = 0, upper = Inf)
  # time2 = Sys.time()
  # print("MLE method costs: ")
  # print(time2 - time1)
  return(res)

}
x = GCR_random(n = 100, lambda = 3, theta = 2, p = 0.3, q = 0.2)
MLE_test_onetime(0.3, 0.2, x)


#Round three
# MLE_helper is  function to compute lambda, theta and mle according to all combinations of p and q in one grid. 
# MLE_helper = function(grid, x) {
#   m = nrow(grid)
#   res = data.frame(matrix(0, m, 3))
#   colnames(res) = c("lambda", "theta", "mle")
#   for (i in 1:m) {
#     p = grid[i, "p"]
#     q = grid[i, "q"]
#
#     curr_res = MLE_test_onetime(p, q, x)
#     res[i, "lambda"] = curr_res$par[1]
#     res[i, "theta"] = curr_res$par[2]
#     res[i, "mle"] = -curr_res$value
#
#   }
#   grid = cbind(grid, res)
#
#   # if you want to check all possible values, print(grid) would help
#   # print(grid)
#
#   target = which.max(grid[, "mle"])
#   ret = list(lambda=grid[target, "lambda"], theta=grid[target, "theta"], p=grid[target, "p"], q=grid[target, "q"], mle=grid[target, "mle"])
#
#   # if you want to check every choice made by choosing highest mle, print(ret) would help
#   # print(ret)
#   return(ret)
# }


# grid_maker is a function that list out all possible combinations of q and p, and put them in a grid.
# grid_maker = function(p_center, q_center, delta, bins) {
#   grid = data.frame(matrix(0, 0, 2))
#   colnames(grid) = c("p", "q")
#   for (i in -bins:bins) {
#     p = p_center + i * delta
#     for (j in -bins:bins) {
#       q = q_center + j * delta
#       if (p >= 0 && p <= 1 && q >= 0 && q <= 1 && p + q <= 1) {
#         grid = rbind(grid, data.frame(p = p, q = q))
#       }
#     }
#   }
#   return(grid)
# }

# MLE_test is a function to implement round 3, with the help of MLE_helper and grid_maker functions.
# MLE_test = function(n = 1000, true_lambda = 1, true_theta = 3, true_p = 0.3, true_q = 0.2, random_seed = 100, bins = 5) {
#   set.seed(random_seed)
#   x = GCR_random(n, true_lambda, true_theta, true_p, true_q)
#   p_center = 0.5
#   q_center = 0.5
#   delta = 0.5 / bins
#   while (delta >= 1e-5) {
#     grid = grid_maker(p_center, q_center, delta, bins)
#     if (nrow(grid) == 0) {
#       break
#     }
#     choice = MLE_helper(grid, x)
#
#     p_center = choice$p
#     q_center = choice$q
#     delta = delta/bins - 1e-5
#   }
#   return(choice)
# }



# An alternative and faster method for round 3, the same idea as round 1.  

MLE_test_another = function(n = 1000, true_lambda = 1, true_theta = 3, true_p = 0.3, true_q = 0.2, random_seed = 100) {
  set.seed(random_seed)
  x = GCR_random(n, true_lambda, true_theta, true_p, true_q)
  l = function(paras) {
    lambda = paras[1]
    theta = paras[2]
    p = paras[3]
    q = paras[4]
    return(-1*MLE_l(lambda, theta, p , q, x = x))
  }
  res = optim(par = c(2, 1, 0.2, 0.4), fn = l, method = "L-BFGS-B", lower = 0, upper = Inf)
  ret = list(lambda=res$par[1], theta=res$par[2], p=res$par[3], q=res$par[4], mle=-res$value)
  return(ret)
  
}

# For MLE_test, if you want to check all values and highest values, you can print(grid) and print(ret) respectively in the MLE_helper function.
# Bins = n,  q and p divide into 2n parts.
# MLE_test(n = 1000, true_lambda = 1, true_theta = 3, true_p = 0.3, true_q = 0.2, random_seed = 100, bins = 20)

# MLE_test_another only has final results
MLE_test_another(n = 1000, true_lambda = 3, true_theta = 2, true_p = 0.3, true_q = 0.2, random_seed = 100)

p = c()
q = c()
l = c()
t = c()
for (i in 1:100) {
  
  res = MLE_test_another(n = 1000, true_lambda = 3, true_theta = 2, true_p = 0.3, true_q = 0.2, random_seed = sample(1:2147483647, 1))
  p = c(p, res$p)
  q = c(q, res$q)
  l = c(l, res$lambda)
  t = c(t, res$theta)


}

