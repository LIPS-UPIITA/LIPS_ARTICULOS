# ============================================================
# VALIDATION OF MATLAB META-ANALYSIS
# sEMG handgrip - sarcopenia
# validation_metafor.R
# R version: 4.6.1
# Package: metafor 5.0-1
#
# Primary implementation:
#   MATLAB R2024b
#
# Independent verification:
#   R / metafor
# ============================================================

rm(list = ls())
cat("\014")

library(metafor)

# ------------------------------------------------------------
# 1. STUDY-LEVEL DATA
# ------------------------------------------------------------

dat <- data.frame(
  
  Study = c(
    "He(2024)_F",
    "Li(2024)_F",
    "He(2024)_M",
    "Li(2024)_M",
    "Sepulveda(2025)",
    "He(2024)",
    "Li(2024)",
    "Hu(2021)"
  ),
  
  Group = c(
    "Elderly_Female",
    "Elderly_Female",
    "Elderly_Male",
    "Elderly_Male",
    "Elderly",
    "Elderly",
    "Elderly",
    "Elderly"
  ),
  
  nControl = c(68,29,30,16,22,98,45,5),
  
  meanControl = c(
    20.4,22.4,32.1,34.1,
    24.5,23.9,26.6,29.6
  ),
  
  sdControl = c(
    3.6,3.0,5.3,5.4,
    7.4,7.0,6.9,8.88
  ),
  
  nCase = c(36,31,19,17,13,55,48,5),
  
  meanCase = c(
    14.8,14.4,24.1,23.8,
    19.8,18.0,17.7,20.0
  ),
  
  sdCase = c(
    2.6,2.8,2.9,2.2,
    5.9,5.2,5.2,3.74
  )
)

# ------------------------------------------------------------
# 2. STUDY-LEVEL EFFECT SIZE
# ------------------------------------------------------------

# Mean Difference:
# Control - Sarcopenia

dat$MD <-
  dat$meanControl - dat$meanCase

# Sampling variance:
#
# Var(MD) =
# SD_control^2 / n_control +
# SD_case^2 / n_case

dat$Variance <-
  (dat$sdControl^2 / dat$nControl) +
  (dat$sdCase^2 / dat$nCase)

dat$SE <- sqrt(dat$Variance)

cat("\n========================================\n")
cat("STUDY-LEVEL EFFECTS\n")
cat("========================================\n")

print(
  dat[,c(
    "Study",
    "Group",
    "MD",
    "SE",
    "Variance"
  )],
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. FUNCTION TO RUN BOTH MODELS
# ------------------------------------------------------------

run_meta <- function(data, group_name) {
  
  d <- subset(data, Group == group_name)
  
  cat("\n\n========================================\n")
  cat("GROUP:", group_name, "\n")
  cat("k =", nrow(d), "\n")
  cat("========================================\n")
  
  # -----------------------------
  # DerSimonian-Laird
  # -----------------------------
  
  fit_DL <- rma(
    yi = MD,
    vi = Variance,
    data = d,
    method = "DL",
    test = "z"
  )
  
  # -----------------------------
  # REML
  # -----------------------------
  
  fit_REML <- rma(
    yi = MD,
    vi = Variance,
    data = d,
    method = "REML",
    test = "z"
  )
  
  cat("\n--- DER SIMONIAN-LAIRD ---\n")
  print(fit_DL)
  
  cat("\n--- REML ---\n")
  print(fit_REML)
  
  # -----------------------------
  # Prediction interval
  # -----------------------------
  
  if (nrow(d) >= 4) {
    
    pred_DL <- predict(
      fit_DL,
      digits = 6
    )
    
    pred_REML <- predict(
      fit_REML,
      digits = 6
    )
    
    cat("\nPrediction interval DL:\n")
    print(pred_DL)
    
    cat("\nPrediction interval REML:\n")
    print(pred_REML)
    
  } else {
    
    cat(
      "\nPrediction intervals NOT reported for k < 4.\n"
    )
  }
  
  # -----------------------------
  # Return results
  # -----------------------------
  
  list(
    data = d,
    DL = fit_DL,
    REML = fit_REML
  )
}

# ------------------------------------------------------------
# 4. THREE PRESPECIFIED INDEPENDENT SUBGROUPS
# ------------------------------------------------------------

res_Female <- run_meta(
  dat,
  "Elderly_Female"
)

res_Male <- run_meta(
  dat,
  "Elderly_Male"
)

res_Elderly <- run_meta(
  dat,
  "Elderly"
)

# ------------------------------------------------------------
# 5. SENSITIVITY ANALYSIS
#    Excluding Hu(2021)
# ------------------------------------------------------------

dat_no_Hu <- subset(
  dat,
  !(Group == "Elderly" & Study == "Hu(2021)")
)

cat("\n\n========================================\n")
cat("SENSITIVITY ANALYSIS\n")
cat("Elderly without Hu(2021)\n")
cat("k =", sum(dat_no_Hu$Group == "Elderly"), "\n")
cat("========================================\n")

fit_DL_no_Hu <- rma(
  yi = MD,
  vi = Variance,
  data = dat_no_Hu,
  subset = Group == "Elderly",
  method = "DL",
  test = "z"
)

fit_REML_no_Hu <- rma(
  yi = MD,
  vi = Variance,
  data = dat_no_Hu,
  subset = Group == "Elderly",
  method = "REML",
  test = "z"
)

cat("\n--- DL without Hu ---\n")
print(fit_DL_no_Hu)

cat("\n--- REML without Hu ---\n")
print(fit_REML_no_Hu)

# ------------------------------------------------------------
# 6. SUMMARY TABLE
# ------------------------------------------------------------

extract_results <- function(fit, group, method) {
  
  data.frame(
    Group = group,
    Method = method,
    k = fit$k,
    Theta = as.numeric(fit$b),
    CI_low = fit$ci.lb,
    CI_high = fit$ci.ub,
    tau2 = fit$tau2,
    Q = fit$QE,
    Q_df = fit$k - 1,
    Q_p = fit$QEp,
    I2 = fit$I2
  )
}

summary_table <- rbind(
  
  extract_results(
    res_Female$DL,
    "Elderly_Female",
    "DL"
  ),
  
  extract_results(
    res_Female$REML,
    "Elderly_Female",
    "REML"
  ),
  
  extract_results(
    res_Male$DL,
    "Elderly_Male",
    "DL"
  ),
  
  extract_results(
    res_Male$REML,
    "Elderly_Male",
    "REML"
  ),
  
  extract_results(
    res_Elderly$DL,
    "Elderly",
    "DL"
  ),
  
  extract_results(
    res_Elderly$REML,
    "Elderly",
    "REML"
  ),
  
  extract_results(
    fit_DL_no_Hu,
    "Elderly_no_Hu",
    "DL"
  ),
  
  extract_results(
    fit_REML_no_Hu,
    "Elderly_no_Hu",
    "REML"
  )
)

cat("\n\n========================================\n")
cat("FINAL VALIDATION TABLE\n")
cat("========================================\n")

print(
  summary_table,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. SAVE VALIDATION TABLE
# ------------------------------------------------------------

write.csv(
  summary_table,
  "metafor_validation_results.csv",
  row.names = FALSE
)

cat("\nValidation table saved as:\n")
cat("metafor_validation_results.csv\n")

# ------------------------------------------------------------
# 8. SAVE COMPLETE R OBJECT
# ------------------------------------------------------------

saveRDS(
  list(
    data = dat,
    Female = res_Female,
    Male = res_Male,
    Elderly = res_Elderly,
    Elderly_no_Hu_DL = fit_DL_no_Hu,
    Elderly_no_Hu_REML = fit_REML_no_Hu
  ),
  "metafor_validation_objects.rds"
)

cat("\nComplete R objects saved as:\n")
cat("metafor_validation_objects.rds\n")

cat("\n========================================\n")
cat("VALIDATION COMPLETED\n")
cat("========================================\n")