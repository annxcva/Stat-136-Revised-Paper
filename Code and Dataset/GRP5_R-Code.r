library(car)
library(dplyr)
library(ggplot2)
library(lmtest)
library(MASS)
library(nortest)
library(olsrr)
library(purrr)
library(tidyr)
library(tseries)
library(whitestrap)

data <- read.csv(
  "[...]/GRP5_Dataset.csv",
  colClasses = rep(c("numeric"), times = 22), header = TRUE
)

data %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
    xlab("Value") + ylab("Frequency") +
    facet_wrap(~ key, scales = "free") +
    geom_histogram() +
    theme_bw()

data %>%
  gather(-PALMA, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = PALMA)) +
  xlab("Value") + ylab("PALMA") +
  geom_point() +
  stat_smooth() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

# Full Model
model <- lm(PALMA ~ ., data)
summary(model)

##### NONLINEARITY
resettest(model)
plot(model, 1:2)
crPlots(model)

### Box-Cox transformation
boxcox <- boxcox(model, lambda = seq(-2,2,0.1))
lambda <- boxcox$x[which.max(boxcox$y)]
data1 <- data %>% mutate(
  PALMA_BC = (PALMA^lambda-1)/lambda
)

# Transformed Model (1)
model1 <- lm(PALMA_BC ~ POV + I(POV^2) + log(UNEMP) + UNDEREMP + log(TRV_FOR) +
  TRV_DOM + I(TRV_DOM^2) + COMPET + GDP_MFG + GDP_CONSTR + GDP_WRT +
  log(GDP_TRANSP) + GDP_AFS + GDP_ICT + log(GDP_FIN) + GDP_REOD + GDP_BUS +
  I(GDP_BUS^2) + sqrt(LDF) + sqrt(LDF_EXP) + log(REV_TAX_RP) +
  sqrt(REV_TAX_BUS) + VIS_GRP + MIN_GRP, data1)
summary(model1)

resettest(model1)
plot(model1, 1:2)
crPlots(model1)

##### AUTOCORRELATION
dwtest(model1)

##### MULTICOLLINEARITY
cor(data1)
ols_vif_tol(model1)
ols_eigen_cindex(model1)

### Centering variables with added polynomial terms
data2 <- data1 %>% mutate(
  POV_CEN = scale(POV, center = TRUE, scale = FALSE),
  TRV_DOM_CEN = scale(TRV_DOM, center = TRUE, scale = FALSE),
  GDP_BUS_CEN = scale(GDP_BUS, center = TRUE, scale = FALSE)
)

# Transformed Model (2)
model2 <- lm(PALMA_BC ~ POV_CEN + I(POV_CEN^2) + log(UNEMP) + UNDEREMP +
  log(TRV_FOR) + TRV_DOM_CEN + I(TRV_DOM_CEN^2) + COMPET + GDP_MFG +
  GDP_CONSTR + GDP_WRT + log(GDP_TRANSP) + GDP_AFS + GDP_ICT + log(GDP_FIN) +
  GDP_REOD + GDP_BUS_CEN + I(GDP_BUS_CEN^2) + sqrt(LDF) + sqrt(LDF_EXP) +
  log(REV_TAX_RP) + sqrt(REV_TAX_BUS) + VIS_GRP + MIN_GRP, data2)
summary(model2)

ols_vif_tol(model2)
ols_eigen_cindex(model2)

### Backward elimination
stepAIC(model2, data2, direction = "backward")

# Reduced Model (3)
model3 <- lm(PALMA_BC ~ POV_CEN + I(POV_CEN^2) + log(UNEMP) + UNDEREMP +
  TRV_DOM_CEN + I(TRV_DOM_CEN^2) + GDP_MFG + GDP_CONSTR + log(GDP_TRANSP) +
  GDP_BUS_CEN + I(GDP_BUS_CEN^2) + sqrt(LDF) + sqrt(REV_TAX_BUS) + VIS_GRP +
  MIN_GRP, data2)
summary(model3)

ols_vif_tol(model3)
ols_eigen_cindex(model3)

##### NONNORMALITY
resid <- resid(model3)
plot(model3, 2)
hist(resid, breaks = 30, main = NULL, xlab = "Residuals", col = "lightblue")

ks.test(resid, "pnorm", mean = mean(resid), sd = sd(resid))
shapiro.test(resid)
ad.test(resid)
cvm.test(resid)
jarque.bera.test(resid)

##### HETEROSKEDASTICITY
gqtest(model3)
white_test(model3)
bptest(model3)

##### OUTLIERS & INFLUENTIAL OBSERVATIONS
influence.measures(model3)
plot(model3, 3:4)
ols_plot_dffits(model3)
ols_plot_dfbetas(model3)

# Final Model
# PALMA_BC ~ POV_CEN + I(POV_CEN^2) + log(UNEMP) + UNDEREMP + TRV_DOM_CEN +
# I(TRV_DOM_CEN^2) + GDP_MFG + GDP_CONSTR + log(GDP_TRANSP) + GDP_BUS_CEN +
# I(GDP_BUS_CEN^2) + sqrt(LDF) + sqrt(REV_TAX_BUS) + VIS_GRP + MIN_GRP

# All diagnostic measures
summary(model3)
resettest(model3)
plot(model3)
crPlots(model3)
dwtest(model3)
ols_vif_tol(model3)
ols_eigen_cindex(model3)
hist(resid, breaks = 30, main = NULL, xlab = "Residuals", col = "lightblue")
ks.test(resid, "pnorm", mean = mean(resid), sd = sd(resid))
shapiro.test(resid)
ad.test(resid)
cvm.test(resid)
jarque.bera.test(resid)
gqtest(model3)
white_test(model3)
bptest(model3)