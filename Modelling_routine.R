##################################################################################################################
#                                                                                                                #
#                           Modelling of Canopy Greenness Temporal Variability                                   #
#                                                                                                                #
##################################################################################################################


###################### Load dataset ########################

# The elemental data derive from the IEFC database of the Catalan Forest and Ecological Inventory. Please,
# check the original database to access full details: http://www.creaf.uab.es/iefc. The data used in this 
# script should be cited and accessed as follows:

# Diniz, E. S., Fernández-Martínez, M., Vayreda, J. y Rodríguez-Penedo, E. (2025). Biomass, elemental 
# composition and environmental characteristics of forests in Catalonia (V1). CORA. Repositorio de Dades 
# de Recerca. https://doi.org/10.34810/data2607



data=read.csv("IEFC_data.csv", sep=";") 
head(data)


############### Select variables to be analyzed ############


# Select the variables to be analyzed and put them into a new data frame. Below, for instance, we will select 
# element concentration data from leaves only and total (all aboveground organs)


library(tidyverse)

df = data%>%
  dplyr::select(Species, CoordenadaUTMX, CoordinateUTMY, Leaves_Biomass, C_total_conc,
	        N_total_conc, P_total_conc, S_total_conc, Ca_total_conc, Mg_total_conc
                K_total_conc, N_Leaves_conc, C_Leaves_conc, P_Leaves_conc, S_Leaves_conc,
                Ca_Leaves_conc, Mg_Leaves_conc, K_Leaves_conc, MAT, MAP, precip_iav, 
                temp_iav, Temp.Season, Temp.Wet.Qt, Prec.Season, Prec.Wet.Month,
                C, pH, N_soil, P_soil, CEC, EVI_mean, EVI_iav)
head(df)



################## Generate the variables of element stocks for the leaves ######################


df$C_st_leaf = df$C_Leaves_conc * df$Leaves_Biomass
df$N_st_leaf = df$N_Leaves_conc * df$Leaves_Biomass
df$P_st_leaf = df$P_Leaves_conc * df$Leaves_Biomass
df$K_st_leaf = df$K_Leaves_conc * df$Leaves_Biomass
df$S_st_leaf = df$S_Leaves_conc * df$Leaves_Biomass
df$Mg_st_leaf = df$Mg_Leaves_conc * df$Leaves_Biomass
df$Ca_st_leaf = df$Ca_Leaves_conc * df$Leaves_Biomass


################################### Construct global models #######################################   

# GAMM models for predicting interanual canopy greenness (EVI_iav) using different leaf stocks, 
# climate, and soil variables. These models are further submitted to Akaike (AIC) selection 
# and model averaging.

### Example: modelling interanual canopy greenness (EVI_iav) from the leaf element stocks ####

library(mgcv)
set.seed(100) 

model<-uGamm(EVI_iav ~ C_st_leaf + N_st_leaf + P_st_leaf + K_st_leaf + S_st_leaf + 
             Ca_st_leaf + Mg_st_leaf + MAT + MAP + precip_iav + temp_iav + 
             Temp.Season + Temp.Wet.Qt + Prec.Season, Prec.Wet.Month + N_soil + 
             P_soil + CEC + s(CoordinateUTM_X, bs ="ds", k=9) + s(CoordinateUTM_Y, bs ="ds",k=10), 
             random=list(Species = ~ 1), data=df, family="gaussian", method = "ML", 
             optimizer=c("outer","newton"), control = gam.control(maxit = 200)) 



### Check the model parameters

summary(model$gam)

###### Check the standardized coefficients ######      

library(MuMIn)
std.coef(model$lme, partial = FALSE)

#### Check the model residuals:
par(mfrow = c(2,2))
gam.check(model$gam)


######## Check the residuals spatial independence #######   

# 1-extract residuals of the model
Y=residuals(model$gam)  
Y

# 2-Visualize a spatial variogram of the residuals

library(gstat)

res=lm(Y ~ 1, data=df)

var.dat_resid <- variogram(Y~1, loc= ~CoordinateUTMX+CoordinateUTMY, data=df)
plot(var.dat_resid)

summary(var.dat_resid)

v.fit = fit.variogram(var.dat_resid, vgm(NA,"Gau",NA,NA))
summary(v.fit)

plot(var.dat_resid, v.fit)


############################ Akaike selection for the global model ################################


#### Selection from the previous for interanual canopy greenness (EVI_iav) from leaf stocks.


##### Make parallel processing #####

library(parallel)
library(doParallel)

ncores <- detectCores() - 10
cl <- makeCluster(ncores)
clusterExport(clust, c("uGamm")) 
clusterEvalQ(cl, library(MuMIn))
clusterExport(cl, varlist = c("df"))

#### Run AIC selection ####

library(MuMIn)


Sys.time() # quantify processing time
a <- Sys.time() # start time

na.action = "na.fail"
aic <- dredge(model, extra = "R^2") 

Sys.time() - a # end time


# Naming and saving the dredge output.

save(aic, file="~/Results/AIC_leaf_stocks_EVI_iav.Rdata")

# Reading the saved output

df2=load("~/Results/AIC_leaf_stocks_EVI_iav.Rdata", ex <- new.env())
ls.str(ex) 
df2=ex$aic
df2 


##################################### Model averaging #############################################

# Extracting further valuable information from models equally robust (i.e., delta < 4) to predict 
# interanual canopy greenness (EVI_iav) from leaf element stocks

ma=summary(model.avg(df2, subset = delta < 4))
ma

# Saving the output
save(ma, file="~/Results/MA_leaf_stocks_EVI_iav.Rdata")

# Reading the saved output
df3=load("~/Results/MA_leaf_stocks_EVI_iav.Rdata", ex2 <- new.env())
ls.str(ex2) 
df4=ex2$ma
df4

df4$coefmat.full # extracting only the full average coefficients
df4$sw # extracting only the variable importance 


########################### Visualization of the results of model averaging ########################

#### Example: visualization of the average models

df3=load("~/Results/MA_leaf_stocks_EVI_iav.Rdata", ex2 <- new.env())
ls.str(ex2) 
df4=ex2$ma
df4

df5=as.data.frame(df4$coefmat.full)
df5

# Attributing the name "Variable" the column "Estimate" 
library(tibble)
df6 <- tibble::rownames_to_column(df5, "Variable")
df6

# converting it to a df object
df7=as.data.frame(df6) 
head(df7)

# selecting the desired variables through indicating the rows. E.g.,
df8=df7[c(2,3,4,5),] 
df8

# Removing the end part of the variables names

df8$Variable<-gsub("_Leaves","",as.character(df8$Variable))
df8

# Reordering the variables' values (increasing or decreasing)
library(dplyr)

df9=df8 %>% arrange(desc(Estimate))
df9  


##### Generating a bar graph ######

library(ggplot2)

g=df9 %>%
  mutate(name = fct_reorder(Variable, desc(Estimate))) %>%
  ggplot( aes(x=Estimate, y=Variable)) +
  geom_bar(stat="identity", fill="#023e8a", alpha=.6, width=.4) +
  ylab("Leaf Element Stocks") +
  xlab("Coefficient")+
  theme_bw()
g

# saving the barplot
ggsave("Fig_MA_leaf_stocks_EVI_iav.jpg", width = 4, height = 4)


################################### SEM (Structural Equation Model) ###############################

df=read.csv("IEFC_data.csv", sep=";")
na.omit(df)
head(df)

########## Create sets of predictors from PCA axes to be used in the SEM equations ##########

# --- PCA replacements via prcomp() ---
pca_elem <- prcomp(df[, c(38, 40, 42, 44)], scale. = TRUE) # select the nutrient stocks chosen 
# in the best model (lowest AIC) in the GAMM model above. 
elem     <- as.data.frame(pca_elem$x) # scores

pca_clim <- prcomp(df[, c(45:52)], scale. = TRUE)
clim     <- as.data.frame(pca_clim$x)

pca_soil <- prcomp(df[, c(55, 56, 57)], scale. = TRUE)
soil     <- as.data.frame(pca_soil$x)

#---- Check it out the factor loading and variables' contributions to each axis


# Example for the set of predictors derived from leaf element stocks from the best GAMM model 
# (lowest AICc) for explaining interannual EVI (EVI_iav)

# Assume the PCA object from prcomp is called pca_elem and use the first axis for representing
# the elemental predictor set

# 1. Eigenvalues (variance per PC)
eigvals <- (pca_elem$sdev)^2 # eigenvalues
sqrt_eig <- pca_elem$sdev # square roots of eigenvalues (sdev)

# 2. Raw loadings (rotation matrix)
raw_loadings <- pca_elem$rotation # columns = PCs, rows = variables

# 3. Correlation-style loadings
# In prcomp, factor loadings = eigenvectors * sqrt(eigenvalue)
corr_loadings <- sweep(raw_loadings, 2, sqrt_eig, "*")  

# 4. Compute variable contributions (% per PC)
var_contrib <- (corr_loadings^2) / rowSums(corr_loadings^2) * 100

# 5. Print results
print("Correlation-style loadings:")
print(corr_loadings)

print("Variable contributions (% per PC):")
print(var_contrib)

############################### Run the SEM model ##############################

##### --- build the SEM data.frame ---
dfsem <- tibble(
  elem    = elem$PC1,
  clim    = clim$PC1,
  soil    = soil$PC1,
  EVI_av  = df$EVI_av,
  EVI_iav = df$EVI_iav,
  sps     = factor(df$Species) 
)

######### ---- Fit non-saturated model ---- ########
fit_psem_reduced1 <- psem(
  lme(EVI_iav ~  soil + elem + EVI_av, random = ~1 | sps, data = dfsem, method="ML"),
  lme(EVI_av  ~  soil + elem,          random = ~1 | sps, data = dfsem, method="ML"),
  lme(soil    ~ clim,                        random = ~1 | sps, data = dfsem, method="ML")
)

summary(fit_psem_reduced1)# Configure the model parameters and coefficients 

plot(fit_psem_reduced1) # ordinary visualization with path diagram


####### Compute total effects ########


## Totals + exact  MC(Monte-Carlo)-based indirects
set.seed(123)
R <- 1000   # Monte-Carlo draws for uncertainty quantification

# 1) pull coefs and detect columns robustly
# coefs() from piecewiseSEM gives regression table; we detect Estimate, SE, t-stat.
cf <- coefs(fit_psem_reduced1)
nm <- tolower(names(cf))

# Estimate: regression coefficients (fixed effects); basis for all effects
est_idx <- which(grepl("^estimate$|^est$|estimate", nm))[1]
if(is.na(est_idx)) stop("Couldn't find an Estimate column in coefs(). Names: ", paste(names(cf), collapse=", "))
est_col <- names(cf)[est_idx]

# SE: standard errors (if available); if missing, reconstructed via Estimate/statistic
se_idx <- which(grepl("std.err|std.error|se\\b|se$|std\\.", nm))
se_col <- if(length(se_idx)) names(cf)[se_idx[1]] else NA

# t/z statistic: allows recovery of SE if not directly provided
stat_idx <- which(grepl("t.value|tvalue|tstat|t.stat|z.value|zstat", nm))
stat_col <- if(length(stat_idx)) names(cf)[stat_idx[1]] else NA

# 2) build edges table
# each edge = one regression coefficient (Predictor → Response)
edges <- data.frame(
  Predictor = as.character(cf$Predictor),
  Response  = as.character(cf$Response),
  Estimate  = suppressWarnings(as.numeric(cf[[est_col]])),
  stringsAsFactors = FALSE
)

# 3) SE extraction
if(!is.na(se_col)) {
  edges$SE <- suppressWarnings(as.numeric(cf[[se_col]]))
} else {
  edges$SE <- NA_real_
}

# recover missing SEs via Estimate / t-stat (since t = beta / SE)
if(any(is.na(edges$SE)) && !is.na(stat_col)) {
  stat_vals <- suppressWarnings(as.numeric(cf[[stat_col]]))
  miss <- which(is.na(edges$SE))
  use_idx <- miss[!is.na(stat_vals[miss]) & abs(stat_vals[miss]) > 1e-12]
  if(length(use_idx) > 0) {
    edges$SE[use_idx] <- abs(edges$Estimate[use_idx] / stat_vals[use_idx])
  }
}

# replace any residual NA/zero SEs with tiny constant (avoid simulation crash)
edges$SE[is.na(edges$SE)] <- 1e-6
edges$SE[edges$SE == 0]   <- 1e-6

# 4) adjacency and point totals
# adjacency A = direct coefficients; total effects = (I - A)^(-1) - I
nodes <- unique(c(edges$Predictor, edges$Response))
n <- length(nodes)
A <- matrix(0, n, n, dimnames = list(nodes, nodes))
for(i in seq_len(nrow(edges))) A[ edges$Predictor[i], edges$Response[i] ] <- edges$Estimate[i]
I <- diag(n)

Tot_point <- tryCatch(solve(I - A) - I, error = function(e) NULL)
if(is.null(Tot_point)) stop("(I - A) is singular for point estimates.")

# extract pairwise total/indirect from point estimates
row_idx <- match(edges$Predictor, nodes)
col_idx <- match(edges$Response, nodes)
Total_point    <- Tot_point[cbind(row_idx, col_idx)]
Indirect_point <- Total_point - edges$Estimate

# 5) Monte-Carlo sims
# simulate edge coefficients ~ Normal(Estimate, SE), propagate into totals
m <- nrow(edges)
tot_sims <- matrix(NA_real_, nrow = R, ncol = m)
direct_sims <- matrix(NA_real_, nrow = R, ncol = m)

for(r in 1:R) {
  ddraws <- rnorm(m, mean = edges$Estimate, sd = edges$SE)
  direct_sims[r, ] <- ddraws
  A_sim <- matrix(0, n, n, dimnames = list(nodes, nodes))
  A_sim[cbind(row_idx, col_idx)] <- ddraws
  S <- tryCatch(solve(I - A_sim) - I, error = function(e) NULL)
  if(!is.null(S)) tot_sims[r, ] <- S[cbind(row_idx, col_idx)]
}

indirect_sims <- tot_sims - direct_sims

# 6) summaries
# For each effect: mean, SD, 95% CI, and empirical p-value (sign test)
total_mean <- colMeans(tot_sims, na.rm = TRUE)
total_sd   <- apply(tot_sims, 2, sd, na.rm = TRUE)
ci_low_tot <- apply(tot_sims, 2, quantile, probs = 0.025, na.rm = TRUE)
ci_high_tot<- apply(tot_sims, 2, quantile, probs = 0.975, na.rm = TRUE)
total_emp_p <- sapply(seq_len(m), function(i){
  v <- tot_sims[, i]; v <- v[!is.na(v)]
  if(length(v)==0) return(NA)
  2 * min(mean(v > 0), mean(v < 0)) # prob(sign different from 0)
})

ind_mean <- colMeans(indirect_sims, na.rm = TRUE)
ind_sd   <- apply(indirect_sims, 2, sd, na.rm = TRUE)
ci_low_ind <- apply(indirect_sims, 2, quantile, probs = 0.025, na.rm = TRUE)
ci_high_ind<- apply(indirect_sims, 2, quantile, probs = 0.975, na.rm = TRUE)
ind_emp_p <- sapply(seq_len(m), function(i){
  v <- indirect_sims[, i]; v <- v[!is.na(v)]
  if(length(v)==0) return(NA)
  2 * min(mean(v > 0), mean(v < 0))
})

# 7) final table
# Combines direct, indirect, total (point + MC), CIs and p-values
out <- data.frame(
  Predictor    = edges$Predictor,
  Response     = edges$Response,
  Direct       = round(edges$Estimate, 6),
  Direct.SE    = round(edges$SE, 6),
  Indirect_pt  = round(Indirect_point, 6),
  Total_pt     = round(Total_point, 6),
  Indirect_MC_mean = round(ind_mean, 6),
  Indirect_MC_CI_low = round(ci_low_ind, 6),
  Indirect_MC_CI_high= round(ci_high_ind, 6),
  Indirect_MC_p = round(ind_emp_p, 4),
  Total_MC_mean   = round(total_mean, 6),
  Total_MC_CI_low = round(ci_low_tot, 6),
  Total_MC_CI_high= round(ci_high_tot, 6),
  Total_MC_p      = round(total_emp_p, 4),
  stringsAsFactors = FALSE
)

out <- out[order(out$Predictor, -abs(out$Total_MC_mean)), ]
rownames(out) <- NULL
print(out, row.names = FALSE)

# save structured results
res_mc <- list(edges = edges,
               Tot_point = Tot_point,
               sims = list(total = tot_sims,
                           direct = direct_sims,
                           indirect = indirect_sims),
               table = out)
res_mc

######### Template option for showing total effects #########


library(ggplot2)

# Plot using Total effects and MC bootstrapped CIs
df_plot <- out

# Filter only effects where Response == "EVI_iav"
df_plot2 <- subset(out, Response == "EVI_iav")

g <- ggplot(out, aes(x = Total_MC_mean, y = Predictor)) +
  geom_point(color = "blue", size = 3) +
  geom_errorbarh(aes(xmin = Total_MC_CI_low, xmax = Total_MC_CI_high), height = 0.2) +
  facet_wrap(~Response, scales = "fixed") +   # <– fixed scale, not free
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Total Effect (bootstrapped CI)",
    y = "Predictor",
    title = "Total Effects (Asymmetric Bootstrapped CIs, Fixed Scale)"
  ) +
  theme_minimal(base_size = 14)

print(g)
