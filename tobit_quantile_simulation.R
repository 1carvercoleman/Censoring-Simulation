# Hi Zac, press CTR + ALT + R to begin :)

# Ideas for graphing: just connect the dots on alpha for bias (found - truth) and variance (found - mean)^2

library(sn)
library(stats)
library(tibble)
library(VGAM)
library(quantreg)

#alpha represents slantness, with a positive value making the bell curve steeper (ie more L-shaped) on the right side of the mean, and a negative value steeper on the left side (0 gives standard normal distribution)

#omega represents scale, with a larger value making a more spread-out distribution (only positive numbers; 1 given standard normal distribution)

#beta refers to true marginal effect of x on y

#cutoff_l refers to a value for a left-side cutoff, cutoff_r for a right-side

#n refers to number of observations

#function will produce a dataframe with y, a column of 1s x_0, and x_1

skewed_df <- function(alpha=0, omega=1, beta=0, cutoff_l=NA, cutoff_r=NA, n=1000){
  x_0 <- rep.int(1,n)
  x_1 <- runif(n,-1,1)
  #To keep things simple, I am making x_1 uniformly distributed around 0
  y <- beta*x_1 + rsn(n=n, xi=0, omega=omega, alpha=alpha, tau=0, dp=NULL)
  if (!is.na(cutoff_l)){
    y <- ifelse(y <= quantile(y, cutoff_l), quantile(y, cutoff_l), y)
  }
  if (!is.na(cutoff_r)){
    y <- ifelse(y >= quantile(y, cutoff_r), quantile(y, cutoff_r), y)
  }
  df <- tibble(y, x_0, x_1)
  return(df)
}
# # 
# # 
# cutoff <- .8
# 
# 
# output1 <- c()
# output <- c()
# cutoff_val <- c()
# 
# #for (cutoff in c(0.7, 0.75, 0.8, 0.85, 0.9, 0.95)) {
#   for (i in 1:500) {
#     df <- skewed_df(alpha=-3, omega=1, beta=.5, cutoff_l = NA, cutoff_r = cutoff, 1000)
#     #my_model <- rq(y ~ x_1, data= df, tau = .5)
#     tobit <- vglm(y ~ x_1, tobit(Upper = max(df$y), Lower = min(df$y)), data = df)
#     
#     #output <- append(output, my_model$coefficients[2])
#     #cutoff_val <- append(cutoff_val, cutoff)
#     output1[i] <- as.numeric(coef(summary(tobit))[3], 1)
#   }
# #}
# 
# ggplot2::ggplot(as.data.frame(output1), ggplot2::aes(x = output1)) +
#   ggplot2::geom_density()
# 
# 
# 
# final_data <- as.data.frame(cbind(cutoff_val, output1))
# rownames(final_data) <- 1:nrow(final_data)
# tapply(final_data$output, as.factor(final_data$cutoff_val), mean)
# 
# ggplot2::ggplot(as.data.frame(output), ggplot2::aes(x = output)) +
#   ggplot2::geom_density()
# 
# mean(output1)
# plot(df$x_1,df$y)
# write.csv(df, "data.csv")




# 
# df_2 <- skewed_df(alpha=-2, omega=3, beta=5, cutoff_l = -10, cutoff_r = 10, 1000)
# plot(df_2$x_1,df_2$y)



# Set beta, number of bootstraps, and sample size 
beta <- .5
boot <- 500
N <- 1000

#Initialize dataframes
quantile_data <- setNames(data.frame(matrix(ncol = 4, nrow = 0))
                          , c("Alpha", "Omega", "Cutoff", "Coefficient"))
tobit_data <- quantile_data

counter <- 1
alpha_range <- c(-3, -2, -1, 0, 1, 2, 3)
omega_range <- 1:7
cutoff_range <- c(.65, .7, .75, .8, .85, .9, .95)

for (alpha in alpha_range) {
  for (omega in omega_range) {
    for (cutoff in cutoff_range) {
      for (i in 1:boot) {
        # Generate dataset
        df <- skewed_df(alpha=alpha, omega=omega, beta=beta, cutoff_l = NA, cutoff_r = cutoff, N)
        
        # Quantile Regression
        my_model <- rq(y ~ x_1, data= df, tau = .5)
        
        # Tobit Regression
        tobit <- vglm(y ~ x_1, tobit(Upper = max(df$y), Lower = min(df$y)), data = df)
        
        # Append coefficient to dataset
        insert_data <- c(alpha, omega, cutoff, my_model$coefficients[2])
        insert_data <- as.data.frame(rbind(insert_data))
        colnames(insert_data) <- colnames(quantile_data)
        rownames(insert_data) <- c()
        quantile_data <- as.data.frame(rbind(quantile_data, insert_data))
        insert_data$Coefficient[1] <- as.numeric(coef(summary(tobit))[3], 1)
        tobit_data <- as.data.frame(rbind(tobit_data, insert_data))
        
        counter <- counter + 1
        fileConn <- file("progress.txt")
        writeLines(paste0("Simulation is ", round((counter / total) * 100, 2), "% complete"), fileConn)
        close(fileConn)
      }
    }
  }
}

# Save datasets as csvs
write.csv(quantile_data, "quantile.csv", row.names = FALSE)
write.csv(tobit_data, "tobit.csv", row.names = FALSE)


print("Simulation Complete. CSV's are located in the following directory:")
print(getwd())
