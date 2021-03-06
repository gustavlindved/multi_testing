# Simulation of different distribution parameters

# Specify different deltas
deltas <- c(0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.35,0.5, 0.7)
sample_size = length(deltas)

# Initialize list to hold results
error_rates <- list()

methods <- c('No_correction', 'Bonferroni', 'Holm', 'BH', 'Q_value')
col_names <- deltas

# Dataframe for total proportion of wrong conclusions
error_rates$total_false <- setNames(data.frame(matrix(NA, nrow=5, ncol=sample_size), row.names = methods), col_names)

# Dataframe for False Discovery rate
error_rates$FDR <- setNames(data.frame(matrix(NA, nrow=5, ncol=sample_size), row.names = methods), col_names)

# Dataframe for Ratio between type 1 and type 2 errors
error_rates$ratio <- setNames(data.frame(matrix(NA, nrow=5, ncol=sample_size), row.names = methods), col_names)

# Simulation parameters
num_simulations <- 15
num_tests <- 4000
alpha <- 0.05

# Run simulations
for (delta_id in 1:sample_size){

  #Distribution parameters
  delta <- deltas[delta_id]
  cat('num_tests is', delta)
  mu1 <- 1
  mu2 <- 1 + (1 * delta)
  sd <- 0.2
  k <- 60
  prop_discoveries <- 0.05
  
  # Create dataframe to store p-values
  p_table <- data.frame(matrix(NA, nrow=num_tests, ncol=num_simulations))
  
  for (j in 1:num_simulations){
    for (i in 1:num_tests){
      if (i < (prop_discoveries*num_tests)){
        p_table[i,j] <- t.test(rnorm(k, mean = mu2, sd = sd), mu = mu1)$p.value
      } else {
        p_table[i,j] <- t.test(rnorm(k, mean = mu1, sd = sd), mu = mu1)$p.value
      }
    }
  }
  
  # Save unadjusted p-values (split into true null and alternative hypotheses)
  split <- prop_discoveries*num_tests
  alt_true <- p_table[1:split,]
  null_true <- p_table[(split+1):num_tests,]
  
  result_table <- data.frame(matrix(NA, nrow=5, ncol=4))
  colnames(result_table) <- c('TP', 'FP', 'TN', 'FN')
  rownames(result_table) <- c('unadjusted', 'bonferroni', 'holm', 'BH', 'qvalue')
  
  result_table$TP[1] <- sum(alt_true < alpha)
  result_table$FN[1] <- sum(alt_true >= alpha)
  result_table$FP[1] <- sum(null_true < alpha)
  result_table$TN[1] <- sum(null_true >= alpha)
  
  # Bonferroni
  result_table$TP[2] <- sum(alt_true < (alpha / num_tests))
  result_table$FN[2] <- sum(alt_true >= (alpha / num_tests))
  result_table$FP[2] <- sum(null_true < (alpha / num_tests))
  result_table$TN[2] <- sum(null_true >= (alpha / num_tests))
  
  # Holm
  p_table_adjust <- data.frame(matrix(NA, nrow=num_tests, ncol=num_simulations))
  for (sim in 1:num_simulations){
    p_table_adjust[,sim] <- p.adjust(p_table[,sim], method = 'holm')
  }
  
  alt_true <- p_table_adjust[1:split,]
  null_true <- p_table_adjust[(split+1):num_tests,]
  
  result_table$TP[3] <- sum(alt_true < (alpha))
  result_table$FN[3] <- sum(alt_true >= (alpha))
  result_table$FP[3] <- sum(null_true < (alpha))
  result_table$TN[3] <- sum(null_true >= (alpha))
  
  # BH
  p_table_adjust <- data.frame(matrix(NA, nrow=num_tests, ncol=num_simulations))
  for (sim in 1:num_simulations){
    p_table_adjust[,sim] <- p.adjust(p_table[,sim], method = 'BH')
  }
  
  alt_true <- p_table_adjust[1:split,]
  null_true <- p_table_adjust[(split+1):num_tests,]
  
  result_table$TP[4] <- sum(alt_true < (alpha))
  result_table$FN[4] <- sum(alt_true >= (alpha))
  result_table$FP[4] <- sum(null_true < (alpha))
  result_table$TN[4] <- sum(null_true >= (alpha))
  
  # Q
  library('qvalue')
  p_table_adjust <- data.frame(matrix(NA, nrow=num_tests, ncol=num_simulations))
  for (sim in 1:num_simulations){
    if (num_tests > 99){
      qobj <- qvalue(p = p_table[,sim])
    } else {
      cat('Using BH approach instead of Q')
      qobj <- qvalue(p = p_table[,sim], pi0 = 1)
    }
    p_table_adjust[,sim] <- qobj$qvalues
  }
  
  alt_true <- p_table_adjust[1:split,]
  null_true <- p_table_adjust[(split+1):num_tests,]
  
  result_table$TP[5] <- sum(alt_true < (alpha))
  result_table$FN[5] <- sum(alt_true >= (alpha))
  result_table$FP[5] <- sum(null_true < (alpha))
  result_table$TN[5] <- sum(null_true >= (alpha))
  
  
  # Average per simulation
  avg_results <- result_table / num_simulations
  
  # Calculate proportion of false conclusions
  error_rates$total_false[,delta_id] <- (avg_results$FP + avg_results$FN) / num_tests
  
  # Calculate FDR
  error_rates$FDR[,delta_id] <- (avg_results$FP / (avg_results$TP + avg_results$FP))
  
  # Calculate Type 1/2 ratio
  error_rates$ratio[,delta_id] <- avg_results$FP / avg_results$FN
}  


# Plotting 
library(ggplot2)
library("reshape2")

# Proportion of wrong conclusions
plot_data <- data.frame(t(error_rates$total_false))
plot_data$delta <- as.numeric(row.names(plot_data))
df_wrong_melt <- melt(plot_data,id.vars = "delta", measure.vars = c("No_correction", 
                                                                             "Bonferroni", 'Holm', 'BH', 'Q_value'))
colnames(df_wrong_melt) <- c("delta", "Method", "value")

ggplot(df_wrong_melt, aes(x = delta, y = value)) + 
  geom_line(aes(color = Method),size = 1.5) + ggtitle("Proportion of wrong conclusions") +
  xlab("Delta") + ylab("Error rate") + scale_color_brewer(palette="Set2") + 
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"), title = element_text(size=18), legend.text = element_text(size = 14))
ggsave('Figures/wrong_delta.pdf')

# False discovery rate
plot_data <- data.frame(t(error_rates$FDR))
plot_data$delta <- as.numeric(row.names(plot_data))
df_FDR_melt <- melt(plot_data,id.vars = "delta", measure.vars = c("No_correction", 
                                                                           "Bonferroni", 'Holm', 'BH', 'Q_value'))
colnames(df_FDR_melt) <- c("delta", "Method", "value")
ggplot(df_FDR_melt, aes(x = delta, y = value)) + 
  geom_line(aes(color = Method),size = 1.5) + ggtitle("False Discovery rate") +
  xlab("delta") + ylab("FDR") + scale_color_brewer(palette="Set2")+
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"), title = element_text(size=18), legend.text = element_text(size = 14))
ggsave('Figures/FDR_delta.pdf')

# Ratio of errors
plot_data <- data.frame(t(error_rates$ratio))
plot_data$delta <- as.numeric(row.names(plot_data))
df_ratio_melt <- melt(plot_data,id.vars = "delta", measure.vars = c("No_correction", 
                                                                             "Bonferroni", 'Holm', 'BH', 'Q_value'))
colnames(df_ratio_melt) <- c("delta", "Method", "value")
ggplot(df_ratio_melt, aes(x = delta, y = value)) + 
  geom_line(aes(color = Method),size = 1.5) + ggtitle("Ratio between type I and type II errors") +
  xlab("delta") + ylab("False positive / False Negative")+ scale_y_log10() + scale_color_brewer(palette="Set2") + 
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"), title = element_text(size=18), legend.text = element_text(size = 14))
ggsave('Figures/ratio_delta.pdf')


