# load package
pacman::p_load(readxl, stringr, haven, dplyr, tidyr, zoo, ggplot2, data.table, broom, survival, splines, tidyverse, parallel, tidyr, ggspatial, cowplot, wesanderson, ggsci, scales,
               purrr, patchwork, tableone, lubridate, rSPARCS, dlnm, gtsummary, officer, flextable, RColorBrewer, knitr, forecast, mvmeta, sf, tsModel, ggforce, ggplot2)

# import demo data
pathogen2011_2019<- read.csv("/data/demo_data.csv")

# Create time-based hierarchical data
out.data<- CXover.data(data=pathogen2011_2019, date="fb_date", ID="id")

# match climate public data
out.data.merge<- Reduce(function(x, y) left_join(x, y), list(out.data, temperature.data.lag))

# Conditional Logistic Regression - Distributed Lag Nonlinear Model
# Create data frame with temperature
temp.data<- out.data.merge[,c("TEMP", paste0("TEMP_lag", 1:21))]
# Generate crossbasis matrix 
cbtemp<-crossbasis(temp.data, lag=c(0,21), argvar=list(fun="bs",degree=2,df=3), arglag=list(knots=logknots(c(0,21),3)))
# Fit the conditional logistic regression model 
dlnm_fit<- clogit(status ~  cbtemp+ ns(QLML_ma7, df = 5) + ns(PM25_ma7, df = 5) + legal_holidays + strata(id), data = out.data.merge, method = 'approximate')
summary(dlnm_fit)

# AIC、BIC
AIC(dlnm_fit)
BIC(dlnm_fit)

# Plot the overall cumulative association
set_min_temp<- quantile(temp.data$TEMP, 0.05, na.rm = T)
set_centered<- quantile(temp.data$TEMP, 0.5, na.rm = T)
set_max_temp<- quantile(temp.data$TEMP, 0.95, na.rm = T)
pred<- crosspred(cbtemp, dlnm_fit, from = set_min_temp, to = set_max_temp, cen = set_centered, cumul=TRUE, bylag = 0.1)
redcum<- crossreduce(cbtemp, dlnm_fit, type="overall", from = set_min_temp, to = set_max_temp, cen = set_centered, lag=21)

plot(redcum, 
     col="black", 
     lty=1, 
     ci.arg=list(col="#7f93a860"),
     lwd=2,
     xlab="Atmospheric temperature",
     ylab="OR (95%CI)",
     # ylim=c(0.0,1.6),
     cex.lab=2,
     cex.axis=2)

# MMT
mmt <- pred$predvar[which.min(pred$allRRfit)]; mmt

# Plot the overall cumulative association
redcum<- crossreduce(cbtemp, dlnm_fit, type="overall", from = set_min_temp, to = set_max_temp, cen = mmt, lag=21)
# cumulative association plot
pdf("./Cumulative exposure–response associations between temperature and flu.pdf", width = 6, height = 5)
par(oma=c(1,1,1,1), mar=c(5,5,5,5))
plot<-plot(redcum, 
           col="#933c30", 
           lty=1, 
           ci.arg=list(col="#e4e7ea"),
           lwd=2,
           xlab="Atmospheric temperature (°C)",
           ylab="OR (95%CI)",
           xaxt = "n", # 
           cex.lab=2,
           cex.axis=1.5)
x_tick_positions <- seq(-15, 25, by = 5)
x_tick_labels <- c("-15","-10", "-5", "0", "5", "10", "15", "20","25")
axis(1, at = x_tick_positions, labels = x_tick_labels, cex.axis = 1.5)
abline(v = c(mmt), lty = 3, col = "gray60", cex.axis = 1)
text(x = 6, y = 3.5, labels = "Temperature=16.5°C", col = "black", cex = 1.5)
dev.off()

# cumulative association and hist plot
out.data.merge2<- subset(out.data.merge, status==1)
pdf("./Hist between temperature and flu.pdf", width = 7, height = 7)
par(oma=c(1,1,1,1), mar=c(5,5,3,5))
plot<-plot(redcum, 
           col="#933c30", 
           lty=1, 
           ci.arg=list(col="#e4e7ea"),
           lwd=2,
           xlab="Atmospheric temperature (°C)",
           ylab="OR (95%CI)",
           xlim=c(-10,30),
           ylim=c(-0.8,3.6),
           axes=F,
           cex.lab=2,
           cex.axis=1.5)
x_tick_positions <- seq(-10, 30, by = 5)
x_tick_labels <- c("-10", "-5", "0", "5", "10", "15", "20", "25", "30")
axis(1, at = x_tick_positions, labels = x_tick_labels, cex.axis = 1.5)
axis(2, at=c(1:7*0.5), cex.axis = 1.5)
text(x = 15, y = 4.5, labels = "Temperature=27°C", col = "black", cex = 1.5)
par(new=T)
hist(subset(out.data.merge2, TEMP>=set_min_temp&TEMP<=set_max_temp)$TEMP,xlim=c(-10,30),ylim=c(0,80000),axes=F,ann=F,col="#e4e7ea",breaks= 40,border="white")
axis(4,at=0:5*5000,cex=1)
abline(v = c(mmt), lty = 3, col = "gray60", cex.axis = 1)
mtext("Freq",4,line=2.5,at=12000,cex=1.5)
dev.off()

# Plot the contour
pred<- crosspred(cbtemp, dlnm_fit, from = set_min_temp, to = set_max_temp, cen = mmt, cumul=TRUE, bylag = 0.1)
pdf("./Lag association pattern from contour plots.pdf", width = 9, height = 7.3)
par(oma=c(1,1,1,1),
    mar=c(5,7,5,5),
    mgp=c(3.5,1,0)) 

nlag <- 21
y <- pred$predvar
x <- seq(0, nlag, 0.1)
z <- t(pred$matRRfit)

# pal <- rev(brewer.pal(11, "BrBG"))
pal <- rev(c("#933c30", "#ad6b62", "#b98179", "#c99e98", "#d7b7b3", "#e5d0ce", "#ffffff", "#e4e7ea", "#c2c8cf", "#a6afb9", "#84909e", "#5f6e81", "#31435c"))
levels <- pretty(z, 20)
col1 <- colorRampPalette(pal[1:7])
col2 <- colorRampPalette(pal[7:13])
cols <- c(col1(sum(levels < 1)), col2(sum(levels > 1)))

filled.contour(x,y,z,
               plot.title = title(xlab="Lag (days)",
                                  ylab="Atmospheric temperature (°C)",
                                  main = "",
                                  col.lab = "black", 
                                  col.main = "black", 
                                  cex.lab=2, 
                                  cex.main=2, 
                                  # font.lab=2, 
                                  font.main=2), 
               col = cols,
               levels = levels, 
               key.title = title("OR", cex.main = 2), 
               key.axes = axis(4, cex.axis=1.5), 
               plot.axes = { axis(1, cex.axis=2, at = 0:nlag, labels = 0:nlag) 
                 axis(2, cex.axis=2)}
)
# mtext(side = 2, at = max(y)*1.1, text = "a", las = 2, cex = 1.2, line = 2)
dev.off()

# plot the slice
quantile(temp.data$TEMP, 0.1, na.rm = T)

pdf("./Lag association pattern from sclice plots.pdf", width = 8, height = 7)
par(oma=c(1,1,1,1), mar=c(5,5,5,5))
pred1 <- crosspred(cbtemp, dlnm_fit, from = set_min_temp, to = set_max_temp, cen=mmt, bylag=0.1, cumul=TRUE)
plot(pred1, 
     "slices", 
     var=-5, 
     col="grey20", 
     xlab="Lag days", 
     ylab="OR (95%CI)", 
     ci.arg=list(col="#7f93a860"), 
     cex.lab=2,
     cex.axis=2,
     lwd=2,
     main="Lag pattern of temperature = -5°C",
     cex.main=1.5)
dev.off()


# xgboost
library(xgboost) 
library(caret)   
library(Metrics)  
library(tibble)  
library(shapviz) 
library(SHAPforxgboost)

features<- c("TEMP_ma7", "QLML_ma7", "SPEEDLML_ma7", "PM1_ma7", "PM25_ma7", "PM10_ma7", "O3_ma7", "NO2_ma7", "SO2_ma7")

X <- as.matrix(out.data.merge[, features, drop = FALSE])
y <- as.numeric(as.character(out.data.merge$status))

dtrain <- xgb.DMatrix(data = X, label = y)

xgb_model <- xgb.train(
  params = list(
    max_depth = 10,
    eta = 1,
    nthread = 4,
    objective = "binary:logistic"
  ),
  data  = dtrain,
  nrounds = 100
)

set.seed(12)
idx <- sample(seq_len(nrow(X)), 5000)
shap_obj <- shapviz(
  xgb_model,
  X_pred = X[idx, , drop = FALSE]
)

shap_mat <- as.matrix(shap_obj$S)
mean_abs_shap <- colMeans(abs(shap_mat))
importance_df <- tibble(
  feature    = names(mean_abs_shap),       
  importance = as.numeric(mean_abs_shap)    
) %>%
  arrange(desc(importance))              

shap_mat_long <- shap_mat %>%
  data.frame() %>% 
  mutate(sample = row_number()) %>% 
  pivot_longer(
    cols = -sample, 
    names_to = "feature",
    values_to = "shap_value"
  )


# randomforest
library(randomForestSRC)

set.seed(123456)

out.data.merge2$status <- as.factor(out.data.merge2$status) 
rf_model <- rfsrc(
  formula = as.factor(status) ~ TEMP_ma7 + QLML_ma7 + SPEEDLML_ma7 + PM1_ma7 + PM25_ma7 + PM10_ma7 + O3_ma7 + NO2_ma7 + SO2_ma7,
  data = out.data.merge2,
  family = "class",
  cluster = id,  
  ntree = 200,
  mtry = 3,
  importance = TRUE, 
  na.action = "na.omit"
)

var_imp <- data.frame(
  variable = names(rf_model[["importance"]][, "all"]),
  importance = rf_model[["importance"]][, "all"],
  row.names = NULL
) %>% arrange(desc(importance)) 


# lstm model
library(tidyverse)
library(lubridate)
library(keras)
library(Metrics)
library(corrplot)
library(ggplot2)
library(caret)


time_steps_list <- c(3)
tuning_results <- tibble()


for (time_steps in time_steps_list) {
  
  y <- model_data$ILI_cases
  
  X_A <- model_data %>% select(starts_with("ILI_cases_lag"), legal_holidays, weekend) %>% as.matrix()
  X_B <- model_data %>% select(starts_with("ILI_cases_lag"), legal_holidays, weekend, starts_with("TEMP_lag"), starts_with("QLML_lag"), starts_with("PM25_lag"), starts_with("O3_lag")) %>% as.matrix()
  X_C <- model_data %>% select(starts_with("ILI_cases_lag"), legal_holidays, weekend, starts_with("TEMP_lag"), starts_with("QLML_lag"), starts_with("PM25_lag"), starts_with("O3_lag"), starts_with("p10_2d_lag"), starts_with("hp_lag")) %>% as.matrix()
  
  scale_vars_A <- setdiff(colnames(X_A), c("legal_holidays", "weekend"))
  scale_vars_B <- setdiff(colnames(X_B), c("legal_holidays", "weekend"))
  scale_vars_C <- setdiff(colnames(X_C), c("legal_holidays", "weekend"))
  
  scale_data <- function(X, scale_vars) {
    scaler <- preProcess(X[, scale_vars, drop = FALSE], method = "range")
    X_scaled <- X
    X_scaled[, scale_vars] <- predict(scaler, X[, scale_vars, drop = FALSE])
    list(X_scaled = X_scaled, scaler = scaler)
  }
  
  scale_A <- scale_data(X_A, scale_vars_A)
  scale_B <- scale_data(X_B, scale_vars_B)
  scale_C <- scale_data(X_C, scale_vars_C)
  
  X_A_scaled <- scale_A$X_scaled
  X_B_scaled <- scale_B$X_scaled
  X_C_scaled <- scale_C$X_scaled
  
  create_lstm_input <- function(X_scaled, y, time_steps) {
    n <- nrow(X_scaled)
    if (n < time_steps) stop(paste0("Insufficient data volume, unable to support", time_steps, "day stride length need", time_steps, "rowss，only", n, "rows"))
    
    X_lstm <- array(0, dim = c(n - time_steps + 1, time_steps, ncol(X_scaled)))
    y_lstm <- numeric(n - time_steps + 1)
    
    for (i in 1:(n - time_steps + 1)) {
      X_lstm[i, , ] <- X_scaled[i:(i + time_steps - 1), ]
      y_lstm[i] <- y[i + time_steps - 1]  
    }
    list(X_lstm = X_lstm, y_lstm = y_lstm)
  }
  
  lstm_A <- create_lstm_input(X_A_scaled, y, time_steps)
  lstm_B <- create_lstm_input(X_B_scaled, y, time_steps)
  lstm_C <- create_lstm_input(X_C_scaled, y, time_steps)
  
  train_ratio <- 0.8
  train_size <- floor(train_ratio * length(lstm_A$y_lstm))
  
  X_A_train <- lstm_A$X_lstm[1:train_size, , ]
  X_A_test <- lstm_A$X_lstm[(train_size + 1):length(lstm_A$y_lstm), , ]
  y_A_train <- lstm_A$y_lstm[1:train_size]
  y_A_test <- lstm_A$y_lstm[(train_size + 1):length(lstm_A$y_lstm)]
  
  X_B_train <- lstm_B$X_lstm[1:train_size, , ]
  X_B_test <- lstm_B$X_lstm[(train_size + 1):length(lstm_B$y_lstm), , ]
  y_B_train <- lstm_B$y_lstm[1:train_size]
  y_B_test <- lstm_B$y_lstm[(train_size + 1):length(lstm_B$y_lstm)]
  
  X_C_train <- lstm_C$X_lstm[1:train_size, , ]
  X_C_test <- lstm_C$X_lstm[(train_size + 1):length(lstm_C$y_lstm), , ]
  y_C_train <- lstm_C$y_lstm[1:train_size]
  y_C_test <- lstm_C$y_lstm[(train_size + 1):length(lstm_C$y_lstm)]
  
  build_lstm_model <- function(input_shape) {
    model <- keras_model_sequential() %>%
      layer_lstm(units = 64, activation = "tanh", return_sequences = FALSE, input_shape = input_shape) %>%
      layer_dense(units = 32, activation = "relu") %>%
      layer_dropout(rate = 0.2) %>%
      layer_dense(units = 1)
    
    model %>% compile(
      optimizer = optimizer_adam(learning_rate = 0.001),
      loss = "mse",
      metrics = c("mae")
    )
    return(model)
  }
  
  early_stop <- callback_early_stopping(monitor = "val_loss", patience = 10, restore_best_weights = TRUE)
  
  model_A <- build_lstm_model(input_shape = c(time_steps, ncol(X_A_scaled)))
  start_time_A <- Sys.time()  
  history_A <- model_A %>% fit(
    x = X_A_train, y = y_A_train,
    validation_data = list(X_A_test, y_A_test),
    epochs = 50, batch_size = 32, callbacks = list(early_stop), verbose = 0
  )
  train_time_A <- difftime(Sys.time(), start_time_A, units = "min") 
  
  model_B <- build_lstm_model(input_shape = c(time_steps, ncol(X_B_scaled)))
  start_time_B <- Sys.time()
  history_B <- model_B %>% fit(
    x = X_B_train, y = y_B_train,
    validation_data = list(X_B_test, y_B_test),
    epochs = 50, batch_size = 32, callbacks = list(early_stop), verbose = 0
  )
  train_time_B <- difftime(Sys.time(), start_time_B, units = "min")
  
  model_C <- build_lstm_model(input_shape = c(time_steps, ncol(X_C_scaled)))
  start_time_C <- Sys.time()
  history_C <- model_C %>% fit(
    x = X_C_train, y = y_C_train,
    validation_data = list(X_C_test, y_C_test),
    epochs = 50, batch_size = 32, callbacks = list(early_stop), verbose = 0
  )
  train_time_C <- difftime(Sys.time(), start_time_C, units = "min")
  
  evaluate_model <- function(model, X_test, y_test, train_time, time_steps, model_type) {
    y_pred <- predict(model, X_test) %>% as.vector()
    tibble(
      time_steps = time_steps,
      model_type = model_type,
      RMSE = rmse(y_test, y_pred),
      MAE = mae(y_test, y_pred),
      R2 = cor(y_test, y_pred)^2,
      MSE = mse(y_test, y_pred),
      train_time_min = as.numeric(train_time),
      Convergence_epoch_number = which.min(history_A$metrics$val_loss)  
    )
  }
  
  current_results <- bind_rows(
    evaluate_model(model_A, X_A_test, y_A_test, train_time_A, time_steps, "A"),
    evaluate_model(model_B, X_B_test, y_B_test, train_time_B, time_steps, "B"),
    evaluate_model(model_C, X_C_test, y_C_test, train_time_C, time_steps, "C")
  )
  
  tuning_results <- bind_rows(tuning_results, current_results)
}

write.csv(tuning_results, "./early_warning_results.csv", row.names = FALSE, fileEncoding = "GBK")

