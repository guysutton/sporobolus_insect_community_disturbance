library(gllvm)
library(mvabund)

# Get community matrix again
com_matrix <- com_matrix[,-7]
com_matrix

# Define species abundance matrix
y <- com_matrix

# Define predictor variables 
x <- org.data %>%
  dplyr::select(disturb = disturb_allow_overwinter,
                season = season)

x <- as.matrix(x)
x

# Model without predictors, poisson errors
fit_null_p <- gllvm::gllvm(y, family = "poisson")
fit_null_nb <- gllvm::gllvm(y, family = "negative.binomial")

# Plot residuals for poisson model
par(mfrow=c(3,2), mar = c(4, 4, 2, 1))
plot(fit_null_p, var.colours = 1)

# Plot residuals for NB model
plot(fit_null_nb, var.colours = 1)

# Clearly the NB is inappropriate (see upper tail of QQ plot).
# For all subsequent models, use Poisson distribution. 

###
# - Fit global model
###

# Specify model with predictor variables 

criterias <- NULL
for(i in 0:5){
  fiti <- gllvm::gllvm(y, x, family = "poisson",
                       num.lv = i, 
                       starting.val = "random",
                       sd.errors = FALSE,
                       formula = ~ disturb + season,
                       seed = 1234)
  criterias[i + 1] <- 
    summary(fiti)$AICc
            names(criterias)[i + 1] = i}

# Fit global model
fit_glb <- gllvm::gllvm(y, x, family = "poisson",
                        num.lv = 2, 
                        formula = ~ disturb * season, 
                        seed = 1234)

# Fit global model
fit_glb_no <- gllvm::gllvm(y, x, family = "poisson",
                        num.lv = 2, 
                        formula = ~ disturb + season, 
                        seed = 1234)

fit_dis <- gllvm::gllvm(y, x, family = "poisson",
                         num.lv = 2, 
                         formula = ~ disturb, 
                         seed = 1234)

fit_sea <- gllvm::gllvm(y, x, family = "poisson",
                        num.lv = 2, 
                        formula = ~ season, 
                        seed = 1234)

# Fit global model
fit_null <- gllvm::gllvm(y, family = "poisson",
                           num.lv = 2, 
                           seed = 1234)

# Is interaction required?
anova(fit_glb, fit_glb_no)

# Is disturb required?
anova(fit_null, fit_dis)

# Extract coefficients 
gllvm::coefplot.gllvm(fit_env)
                      #xlim.list = list(NULL, NULL, c(-4, 4),
                      #                 mfrow = c(1,1)))
# 
coef(fit_env)
confint(fit_env)
