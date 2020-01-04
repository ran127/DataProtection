library(mitools)
library(tidyverse)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
#*****************************************Simulation Function***********************************************************#
# r as target correlation between U1 and U2
# n as the total number of simulation
# k as the number of observations for each X(as data augmentation part)
# var as the variance for error added to models(chosen from {1, 3, 5})
setwd("~Desktop/FALL2019/MATH410")
simulation <- function(iteration){
  r <- 0.5                            # Target (Spearman) correlation
  n <- 500                            # Number of samples

  # Generate U
  gen.gauss.cop <- function(r, n){
    rho <- 2 * sin(r * pi/6)          # Pearson correlation
    P <- toeplitz(c(1, rho))          # Correlation matrix
    d <- nrow(P)                      # Dimension
    U <- pnorm(matrix(rnorm(n*d), ncol = d) %*% chol(P))
    return(U)
  }
  
  # U1 and U2 are uniformly correlated
  U <- gen.gauss.cop(r, n)
  U_1<-U[,1]
  U_2<-U[,2]

  # Z_cov has uniform distribution with (-2, 2)
  Z_cov <- runif(n, -2, 2)

  # Generate X1, X2
  miu <- Z_cov^2 + U_1
  X1 <- rnorm(n, miu, 1)
  X2 <- rnorm(n, miu, 1)

  # Z_resp is from bernoulli
  l <- 1 + X1 + X2 + Z_cov
  p_z <- pnorm(l)                     # calculate the cdf of l
  Z_resp <- rbinom(n, 1, p_z)

  # Mask Model
  k = 15                              # number of oberservation for each X fixed in paper
  var = 1                             # variance for error added to model {1, 3, 5}
  X1p = t(rep(1, k) %*% t.default(X1))
  W1 = X1p + matrix(rnorm(500*k, 0, var), ncol=k)

  X2p = t(rep(1, k) %*% t.default(X2))
  W2 = X2p + matrix(rnorm(500*k, 0, var), ncol=k)


  # Stan model
  library(rstan)
  rstan_options(auto_write = TRUE)
  model_Mx <- stan_model('Model.stan')
  fit_1 <- sampling(model_Mx, list(N=n, K=k, Z_cov=Z_cov,Z_resp = Z_resp,
                                  W1 = W1, W2=W2, X1=X1, X2=X2,alpha = 1))

  # For loop to create 5 synthetic datasets from draws 100,200,300,400,500 from sampler
  list_of_draws <- extract(fit_1)
  X1data_list <- list(c(list_of_draws$X1_new[100,]), c(list_of_draws$X1_new[200,]), c(list_of_draws$X1_new[300,]),
                      c(list_of_draws$X1_new[400,]), c(list_of_draws$X1_new[500,]))
  X2data_list <- list(c(list_of_draws$X2_star_new[100,]), c(list_of_draws$X2_star_new[200,]), 
                      c(list_of_draws$X2_star_new[300,]),c(list_of_draws$X2_star_new[400,]), 
                      c(list_of_draws$X2_star_new[500,]))
  mylist<-list()
  for(i in 1:5){
    X_1_result <- X1data_list[[i]]
    X_2_result <- X2data_list[[i]]
    mylist[[i]]<- glm(Z_resp~Z_cov+X_1_result+X_2_result+1, family=binomial(link="probit"))
  }
  
  # True values' result
  trueresult <- as_tibble(coef(summary(glm(Z_resp~Z_cov+X1+X2, family=binomial(link="probit")))), rownames="Parameter")%>% 
    rename(TrueEstimate = Estimate, TrueSD = `Std. Error` ) %>% select(TrueEstimate, TrueSD)
  trueresult <- transform(trueresult, TrueCILow = TrueEstimate - 1.96*TrueSD)
  trueresult <- transform(trueresult, TrueCIHigh = TrueEstimate + 1.96*TrueSD)
  
  # Synthetic datasets' result
  sampleresult <- as_tibble(summary(MIcombine(mylist)),rownames="Parameter") %>% rename(SynEstimate = results, SynSD= se,
                                                                                        SynLower = `(lower`, SynUpper = `upper)`)
  sampleresult <- transform(sampleresult, SynCILow = SynEstimate - 1.96*SynSD)
  sampleresult <- transform(sampleresult, SynCIHigh = SynEstimate + 1.96*SynSD)
  
  finalresult <- bind_cols(sampleresult, trueresult)
  
  # Calculate the standard difference
  finalresult <- transform(finalresult, standardDiff = abs(TrueEstimate - SynEstimate)/TrueSD)
  
  # Calculate data utility
  data_utility <- select(finalresult, Parameter, SynCILow, SynCIHigh, TrueCILow, TrueCIHigh)
  
  data_utility <- data_utility %>% mutate(DU = (min(TrueCIHigh, SynCIHigh) - max(TrueCILow, SynCILow))/(TrueCIHigh-TrueCILow) +
                            (min(TrueCIHigh, SynCIHigh) - max(TrueCILow, SynCILow))/(SynCIHigh-SynCILow))
  finalresult <- finalresult %>% mutate(DU= data_utility$DU)
  
  match_rate <- match_rate_fn(X1, X2, X1data_list, X2data_list)
  finalresult <- finalresult %>% mutate(MR = match_rate)
  # Generate final result by adding iterations
  finalresult <- finalresult %>% mutate(Iter=iteration)  
}

#Match Rate function
match_rate_fn <- function(X1, X2, X1data_list, X2data_list){
  match_array <- array(rep(0, 500*500*5), dim=c(500, 500, 5))
  d <- 0.1
  # initalize the true range array X[1][1] represents the lower bound for X1's 1st record, 
  # X[1][4] represents the upper bound for X2's 1st record
  true_array <- array(rep(0, 500*4), dim=c(500, 4))
  X1_array <- array(X1, dim = c(500,1))
  X2_array <- array(X2, dim = c(500,1))
  for (i in 1:500){
    true_array[i,1] = X1_array[i] - d
    true_array[i,2] = X1_array[i] + d
    true_array[i,3] = X2_array[i] - d
    true_array[i,4] = X2_array[i] + d
  }
  
  for (i in 1:5){
    X1_syn = X1data_list[[i]]
    X2_syn = X2data_list[[i]]
    for (j in 1:500){
      X1_j <- X1_syn[j]
      X2_j <- X2_syn[j]
      count <- 0
      for (k in 1:500){
        if(true_array[k,1] <= X1_j  && X1_j<= true_array[k,2] && true_array[k,3] <= X2_j && X2_j<= true_array[k, 4]){
          match_array[j, k, i] <- 1
        }
      }
    }
  }
  
  for (i in 1:5){
    for (j in 1:500){
      sum <- sum(match_array[j,,i])
      for (k in 1:500){
        if(match_array[j,k,i] ==1){
          match_array[j,k,i] <-1/sum
        }
      }
    }
  }
  twod_matching <- array(rep(0,500*500), dim=c(500, 500))
  for (j in 1:500){
    for(k in 1:500){
      twod_matching[j,k] = mean(match_array[j,k,])
    }
  }
  result <- rep(0, 500)
  match_rate <-0
  for (i in 1:500){
    max_prob <- max(twod_matching[i,])
    index_len <- length(which(twod_matching[i,]== max(twod_matching[i,])))
    result[i] <- max_prob*index_len
    if(result[i] == 1){
      match_rate <- match_rate +1
    }
  }
  match_rate <- match_rate/500
}




results <- map_dfr(1:25, simulation)
# Analyze the result
# Part 1 Data Utility- Overlap
results %>% group_by(Parameter)%>%
  summarise(Avg_DU = mean(DU), Max_DU = max(DU), Min_DU = min(DU)) 
ggplot(results, aes(x=Parameter, y=DU)) + stat_boxplot(geom = 'errorbar', width = 0.25) + 
  geom_boxplot(fill ='lightpink')

# Part 2 Data Utility- Standard Difference
results %>% group_by(Parameter) %>%
  summarise(Avg_SD = mean(standardDiff), Max_SD=max(standardDiff), Min_SD = min(standardDiff))
ggplot(results, aes(x=Parameter, y=standardDiff)) + stat_boxplot(geom = 'errorbar', width = 0.25) + 
  geom_boxplot(fill ='lightpink')

# Part 3 Analyzing the Parameters
X1_res <- results %>% filter(Parameter == "X_1_result")
X2_res <- results %>% filter(Parameter == "X_2_result")
interc <- results %>% filter(Parameter == "(Intercept)")
zcov_res <- results %>% filter(Parameter == "Z_cov")

p1 <- ggplot(X1_res, aes(y = SynEstimate)) + stat_boxplot(geom = 'errorbar', width = 0.25) + geom_boxplot(fill ='lightpink')
p2 <- ggplot(X1_res, aes(y = TrueEstimate)) + stat_boxplot(geom = 'errorbar', width = 0.25) + geom_boxplot(fill ='lightpink')
grid.arrange(p1, p2, nrow=1)

p3 <- ggplot(X2_res, aes(y = SynEstimate)) + stat_boxplot(geom = 'errorbar', width = 0.25) + geom_boxplot(fill ='lightpink')
p4 <- ggplot(X2_res, aes(y = TrueEstimate)) + stat_boxplot(geom = 'errorbar', width = 0.25) + geom_boxplot(fill ='lightpink')
grid.arrange(p3, p4, nrow=1)

p5 <- ggplot(interc, aes(y = SynEstimate)) + stat_boxplot(geom = 'errorbar', width = 0.25) + geom_boxplot(fill ='lightpink')
p6 <- ggplot(interc, aes(y = TrueEstimate)) + stat_boxplot(geom = 'errorbar', width = 0.25) + geom_boxplot(fill ='lightpink')
grid.arrange(p5, p6, nrow=1)

p7 <- ggplot(zcov_res, aes(y = SynEstimate)) + stat_boxplot(geom = 'errorbar', width = 0.25) + geom_boxplot(fill ='lightpink')
p8 <- ggplot(zcov_res, aes(y = TrueEstimate)) + stat_boxplot(geom = 'errorbar', width = 0.25) + geom_boxplot(fill ='lightpink')
grid.arrange(p7, p8, nrow=1)

p9 <- ggplot(X1_res, aes(x="SynX1", y=SynEstimate)) + geom_point()
p10 <- ggplot(X1_res, aes(x="TrueX1", y=TrueEstimate)) + geom_point()
grid.arrange(p9, p10, nrow=1)

p11 <- ggplot(X2_res, aes(x="SynX2", y=SynEstimate)) + geom_point()
p12 <- ggplot(X2_res, aes(x="TrueX2", y=TrueEstimate)) + geom_point()
grid.arrange(p11, p12, nrow=1)

p12 <- ggplot(interc, aes(x="SynIntercept", y=SynEstimate)) + geom_point()
p13 <- ggplot(interc, aes(x="Trueintercept", y=TrueEstimate)) + geom_point()
grid.arrange(p12, p13, nrow=1)

p14 <- ggplot(zcov_res, aes(x="SynZcov", y=SynEstimate)) + geom_point()
p15 <- ggplot(zcov_res, aes(x="TrueZcov", y=TrueEstimate)) + geom_point()
grid.arrange(p14, p15, nrow=1)

# Part 4 Analyzing the True Match Rate
subresults <- results %>% filter(Parameter == 'Z_cov')
ggplot(subresults, aes(y=MR)) + stat_boxplot(geom = 'errorbar', width = 0.25) + 
  geom_boxplot(fill ='lightpink')
subresults %>%
  summarise(Avg_MR = mean(MR), Max_DU = max(MR), Min_DU = min(MR)) 







