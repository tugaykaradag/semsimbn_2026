# model 3 200 icin 
{
  bas <- Sys.time()
#gozlem sayisini degistir sadece
library(lavaan)
model_pop <- '
# Measurement model
B =~ 0.8*B1 + 0.8*B2 + 0.8*B3
D =~ 0.8*D1 + 0.8*D2 + 0.8*D3
E =~ 0.8*E1 + 0.8*E2 + 0.8*E3

# Structural model (Model3: E->B, B->D)
B ~ 0.6*E
D ~ 0.6*B

# (Optional) latent variances
B ~~ 1*B
D ~~ 1*D
E ~~ 1*E

# (Optional) indicator residual variances (since loading=0.8 => residual var=0.36)
B1 ~~ 0.36*B1
B2 ~~ 0.36*B2
B3 ~~ 0.36*B3
D1 ~~ 0.36*D1
D2 ~~ 0.36*D2
D3 ~~ 0.36*D3
E1 ~~ 0.36*E1
E2 ~~ 0.36*E2
E3 ~~ 0.36*E3
'

model_fit <- '
# Measurement model (free loadings)
B =~ B1 + B2 + B3
D =~ D1 + D2 + D3
E =~ E1 + E2 + E3

# Structural model (free regressions) - Model3
B ~ E
D ~ B
'
#fit control function

# 1) istersen "tum yapisal yollar anlamli mi" kontrol??:
#    - sadece regressions: op == "~"
#    - ister measurement ( "=~") da dahil edebilirsin (parametreden a??t??m)
max_pvalue_check <- function(fit, include_measurement = FALSE) {
  pe <- parameterEstimates(fit)
  ops_keep <- if (include_measurement) c("~", "=~") else c("~")
  pe <- pe[pe$op %in% ops_keep, ]
  # Sabit/variance gibi seyleri karsilastirmamak icin:
  pe <- pe[is.finite(pe$pvalue), ]
  if (nrow(pe) == 0) return(NA_real_)
  max(pe$pvalue, na.rm = TRUE)
}

# 2) Fit filtresi (senin kosullarin; 0.95 -> 0.90 guncellendi)
pass_fit_full <- function(fit,
                          p_max = 0.05,
                          srmr_max = 0.05,
                          rmsea_max = 0.05,
                          min_fit = 0.90,
                          chisq_df_max = 3,
                          include_measurement_p = FALSE) {
  
  # Eger fit test="none" ile uretildiyse fitMeasures alinamaz -> direkt reject
  if (isTRUE(lavInspect(fit, "options")$test == "none")) {
    return(list(ok = FALSE, f_meas = NA, pmax = NA_real_))
  }
  
  f_meas <- fitMeasures(fit, c("srmr","gfi","agfi","cfi","nfi","nnfi","ifi","rmsea","chisq","df"))
  pmax <- max_pvalue_check(fit, include_measurement = include_measurement_p)
  
  ok <- is.finite(pmax) && (pmax <= p_max) &&
    is.finite(f_meas["srmr"]) && (f_meas["srmr"] <= srmr_max) &&
    is.finite(f_meas["rmsea"]) && (f_meas["rmsea"] <= rmsea_max) &&
    all(is.finite(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")])) &&
    all(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")] >= min_fit) &&
    all(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")] <= 1) &&
    is.finite(f_meas["chisq"]) && is.finite(f_meas["df"]) &&
    (f_meas["chisq"] / f_meas["df"] < chisq_df_max)
  
  list(ok = ok, f_meas = f_meas, pmax = pmax)
}

#produce - fit SEM - if yes save

generate_sem_datasets_filtered <- function(model_pop, model_fit, n, K = 1000,
                                           max_attempts = 200000,
                                           seed = NULL,
                                           verbose_every = 50,
                                           include_measurement_p = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  accepted <- vector("list", K)
  log_df <- data.frame(
    accepted_id = integer(0),
    attempt = integer(0),
    pmax = numeric(0),
    srmr = numeric(0),
    gfi = numeric(0),
    agfi = numeric(0),
    cfi = numeric(0),
    nfi = numeric(0),
    nnfi = numeric(0),
    ifi = numeric(0),
    rmsea = numeric(0),
    chisq = numeric(0),
    df = numeric(0),
    chisq_df = numeric(0)
  )
  
  acc <- 0
  attempt <- 0
  
  while (acc < K) {
    attempt <- attempt + 1
    if (attempt > max_attempts) stop("max_attempts asildi. Esikleri gevset veya max_attempts artir.")
    
    dat <- tryCatch(simulateData(model_pop, sample.nobs = n), error = function(e) NULL)
    if (is.null(dat)) next
    
    fit <- tryCatch(
      lavaan::sem(model_fit, data = dat, std.lv = TRUE, test = "standard"),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    # <-- kritik kontrol
    topt <- lavInspect(fit, "options")$test
    if (isTRUE(topt == "none")) {
      # bu run'da fitMeasures mumkun degil, reject
      next
    }
    
    chk <- pass_fit_full(fit, include_measurement_p = include_measurement_p)
    
    if (isTRUE(chk$ok)) {
      acc <- acc + 1
      accepted[[acc]] <- dat
      
      fm <- chk$f_meas
      log_df <- rbind(log_df, data.frame(
        accepted_id = acc,
        attempt = attempt,
        pmax = chk$pmax,
        srmr = fm["srmr"],
        gfi  = fm["gfi"],
        agfi = fm["agfi"],
        cfi  = fm["cfi"],
        nfi  = fm["nfi"],
        nnfi = fm["nnfi"],
        ifi  = fm["ifi"],
        rmsea = fm["rmsea"],
        chisq = fm["chisq"],
        df = fm["df"],
        chisq_df = fm["chisq"]/fm["df"]
      ))
      
      if (acc %% verbose_every == 0) {
        cat(sprintf("Accepted %d/%d (attempt %d) | pmax=%.4f CFI=%.3f RMSEA=%.3f SRMR=%.3f\n",
                    acc, K, attempt, chk$pmax, fm["cfi"], fm["rmsea"], fm["srmr"]))
      }
    }
  }
  
  list(datasets = accepted, fitlog = log_df, n = n, K = K, total_attempts = attempt)
}


#run example
out_200 <- generate_sem_datasets_filtered(model_pop, model_fit, n = 200, K = 1000)

# Kac denemede 999+1 kabul topladi
out_200$total_attempts

# Kabul edilenlerin fit ozetleri
summary(out_200$fitlog[, c("pmax","cfi","rmsea","srmr","chisq_df")])

#aggregation
# out_200$datasets : list of data.frames (A1..E3)

aggregate_to_latent_means <- function(df) {
  data.frame(
    B = rowMeans(df[, c("B1","B2","B3")], na.rm = TRUE),
    D = rowMeans(df[, c("D1","D2","D3")], na.rm = TRUE),
    E = rowMeans(df[, c("E1","E2","E3")], na.rm = TRUE)
  )
}



# BN pipeline'a girecek liste (her biri A..E olan 200 gozlemlik data.frame)
data_list_bn_200 <- lapply(out_200$datasets, aggregate_to_latent_means)

#veri seti olusturuldu. artik algoritmalari deneyelim

# ground truth (Model3)
true_edges <- c(
  "E B",
  "B D"
)

nodes <- c("B","D","E")
all_edges <- as.vector(outer(nodes, nodes, FUN = paste))
all_edges <- all_edges[!grepl("(B B|D D|E E)", all_edges)]


#algorithm wrappers
library(bnlearn)

algo_fns <- list(
  "H2PC"   = function(d) h2pc(d),
  "MMHC"   = function(d) mmhc(d),
  "RSMAX2" = function(d) rsmax2(d),
  "HC"     = function(d) hc(d),
  "PCS"    = function(d) pc.stable(d),
  "Tabu"   = function(d) tabu(d),
  "HPC"    = function(d) hpc(d),
  "SPC"    = function(d) si.hiton.pc(d),
  "IAMB"   = function(d) iamb(d),
  "I-IAMB" = function(d) inter.iamb(d),
  "IAMB-F" = function(d) iamb.fdr(d),
  "F-IAMB" = function(d) fast.iamb(d),
  "MMPC"   = function(d) mmpc(d),
  "GS"     = function(d) gs(d)
)

#extract estimated edges

extract_edges <- function(bn_fit) {
  am <- amat(bn_fit)
  idx <- which(am == 1, arr.ind = TRUE)
  if (nrow(idx) == 0) return(character(0))
  paste(rownames(am)[idx[,1]], colnames(am)[idx[,2]])
}

#metrics for single dataset
compute_shd_directed <- function(pred, truth) {
  # pred, truth: character vector "X Y"
  
  all_nodes <- unique(unlist(strsplit(c(pred, truth), " ")))
  all_edges <- as.vector(outer(all_nodes, all_nodes, paste))
  all_edges <- all_edges[!grepl("^(.+) \\1$", all_edges)]
  
  FP <- sum(pred %in% setdiff(all_edges, truth))
  FN <- sum(truth %in% setdiff(all_edges, pred))
  
  FP + FN
}

compute_metrics <- function(pred, truth) {
  
  TP <- sum(pred %in% truth)
  FP <- sum(pred %in% setdiff(all_edges, truth))
  FN <- sum(!(truth %in% pred))
  TN <- length(all_edges) - TP - FP - FN
  
  ACC <- (TP + TN) / (TP + FP + FN + TN)
  SEN <- ifelse(TP+FN==0, NA, TP/(TP+FN))
  PRE <- ifelse(TP+FP==0, NA, TP/(TP+FP))
  SPE <- ifelse(TN+FP==0, NA, TN/(TN+FP))
  F1  <- ifelse(is.na(PRE+SEN), NA, 2*PRE*SEN/(PRE+SEN))
  
  Pc <- (((TP+FN)*(TP+FP))+((TN+FN)*(TN+FP)))/(TP+FP+FN+TN)^2
  CK <- (ACC - Pc)/(1 - Pc)
  
  MCC <- ifelse(
    (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)==0, NA,
    (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  )
  
  SHD <- compute_shd_directed(pred, truth)
  
  c(ACC,SEN,PRE,SPE,F1,CK,MCC,SHD)
}


#999+1x14  run
run_all_algorithms <- function(data_list) {
  
  metrics <- c("ACC","SEN","PRE","SPE","F1","CK","MCC","SHD")
  results <- list()
  
  for (alg in names(algo_fns)) {
    cat("Running:", alg, "\n")
    mat <- matrix(NA, nrow = length(data_list), ncol = length(metrics))
    
    for (i in seq_along(data_list)) {
      fit <- tryCatch(algo_fns[[alg]](data_list[[i]]), error = function(e) NULL)
      if (is.null(fit)) next
      
      pred <- extract_edges(fit)
      mat[i,] <- compute_metrics(pred, true_edges)
    }
    
    colnames(mat) <- metrics
    results[[alg]] <- as.data.frame(mat)
  }
  results
}

res_200 <- run_all_algorithms(data_list_bn_200)

#create final table
fmt_num <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_character_, length(x))
  
  ok <- !is.na(x)
  
  # |x| < 1  ??? 2 decimals, no leading zero
  small <- ok & abs(x) < 1
  out[small] <- sub(
    "^0",
    "",
    formatC(round(x[small], 2), format = "f", digits = 2)
  )
  
  # |x| ??? 1 ??? 1 decimal
  big <- ok & abs(x) >= 1
  out[big] <- formatC(round(x[big], 1), format = "f", digits = 1)
  
  out
}

make_summary_table <- function(res) {
  
  metrics <- c("ACC","SEN","PRE","SPE","F1","CK","MCC","SHD")
  out <- data.frame(Algorithm = names(res))
  
  for (m in metrics) {
    vals <- sapply(res, function(x) mean(x[[m]], na.rm=TRUE))
    sds  <- sapply(res, function(x) sd(x[[m]], na.rm=TRUE))
    if (m == "SHD") {
      ranks <- rank(vals, ties.method = "average")
    } else {
      ranks <- rank(-vals, ties.method = "average")
    }
    
    out[[m]] <- paste0(fmt_num(vals), " (", fmt_num(sds), ")")
    out[[paste0("R_",m)]] <- ranks
  }
  
  rank_cols <- grep("^R_", names(out))
  out$Overall_Rank <- rowMeans(out[, rank_cols])
  
  out[order(out$Overall_Rank), ]
}

final_table_200 <- make_summary_table(res_200)
final_table_200
library(writexl)
write_xlsx(final_table_200,"D:/Makaleler/SEM Simulasyon/Revs/Rev13 Array/Added Analyses/model_3_final_table_200.xlsx")
}

#500 icin

#gozlem sayisini degistir sadece
library(lavaan)
model_pop <- '
# Measurement model
B =~ 0.8*B1 + 0.8*B2 + 0.8*B3
D =~ 0.8*D1 + 0.8*D2 + 0.8*D3
E =~ 0.8*E1 + 0.8*E2 + 0.8*E3

# Structural model (Model3: E->B, B->D)
B ~ 0.6*E
D ~ 0.6*B

# (Optional) latent variances
B ~~ 1*B
D ~~ 1*D
E ~~ 1*E

# (Optional) indicator residual variances (since loading=0.8 => residual var=0.36)
B1 ~~ 0.36*B1
B2 ~~ 0.36*B2
B3 ~~ 0.36*B3
D1 ~~ 0.36*D1
D2 ~~ 0.36*D2
D3 ~~ 0.36*D3
E1 ~~ 0.36*E1
E2 ~~ 0.36*E2
E3 ~~ 0.36*E3
'

model_fit <- '
# Measurement model (free loadings)
B =~ B1 + B2 + B3
D =~ D1 + D2 + D3
E =~ E1 + E2 + E3

# Structural model (free regressions) - Model3
B ~ E
D ~ B
'
#fit control function

# 1) istersen "tum yapisal yollar anlamli mi" kontrol??:
#    - sadece regressions: op == "~"
#    - ister measurement ( "=~") da dahil edebilirsin (parametreden a??t??m)
max_pvalue_check <- function(fit, include_measurement = FALSE) {
  pe <- parameterEstimates(fit)
  ops_keep <- if (include_measurement) c("~", "=~") else c("~")
  pe <- pe[pe$op %in% ops_keep, ]
  # Sabit/variance gibi seyleri karsilastirmamak icin:
  pe <- pe[is.finite(pe$pvalue), ]
  if (nrow(pe) == 0) return(NA_real_)
  max(pe$pvalue, na.rm = TRUE)
}

# 2) Fit filtresi (senin kosullarin; 0.95 -> 0.90 guncellendi)
pass_fit_full <- function(fit,
                          p_max = 0.05,
                          srmr_max = 0.05,
                          rmsea_max = 0.05,
                          min_fit = 0.90,
                          chisq_df_max = 3,
                          include_measurement_p = FALSE) {
  
  # Eger fit test="none" ile uretildiyse fitMeasures alinamaz -> direkt reject
  if (isTRUE(lavInspect(fit, "options")$test == "none")) {
    return(list(ok = FALSE, f_meas = NA, pmax = NA_real_))
  }
  
  f_meas <- fitMeasures(fit, c("srmr","gfi","agfi","cfi","nfi","nnfi","ifi","rmsea","chisq","df"))
  pmax <- max_pvalue_check(fit, include_measurement = include_measurement_p)
  
  ok <- is.finite(pmax) && (pmax <= p_max) &&
    is.finite(f_meas["srmr"]) && (f_meas["srmr"] <= srmr_max) &&
    is.finite(f_meas["rmsea"]) && (f_meas["rmsea"] <= rmsea_max) &&
    all(is.finite(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")])) &&
    all(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")] >= min_fit) &&
    all(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")] <= 1) &&
    is.finite(f_meas["chisq"]) && is.finite(f_meas["df"]) &&
    (f_meas["chisq"] / f_meas["df"] < chisq_df_max)
  
  list(ok = ok, f_meas = f_meas, pmax = pmax)
}

#produce - fit SEM - if yes save

generate_sem_datasets_filtered <- function(model_pop, model_fit, n, K = 1000,
                                           max_attempts = 200000,
                                           seed = NULL,
                                           verbose_every = 50,
                                           include_measurement_p = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  accepted <- vector("list", K)
  log_df <- data.frame(
    accepted_id = integer(0),
    attempt = integer(0),
    pmax = numeric(0),
    srmr = numeric(0),
    gfi = numeric(0),
    agfi = numeric(0),
    cfi = numeric(0),
    nfi = numeric(0),
    nnfi = numeric(0),
    ifi = numeric(0),
    rmsea = numeric(0),
    chisq = numeric(0),
    df = numeric(0),
    chisq_df = numeric(0)
  )
  
  acc <- 0
  attempt <- 0
  
  while (acc < K) {
    attempt <- attempt + 1
    if (attempt > max_attempts) stop("max_attempts asildi. Esikleri gevset veya max_attempts artir.")
    
    dat <- tryCatch(simulateData(model_pop, sample.nobs = n), error = function(e) NULL)
    if (is.null(dat)) next
    
    fit <- tryCatch(
      lavaan::sem(model_fit, data = dat, std.lv = TRUE, test = "standard"),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    # <-- kritik kontrol
    topt <- lavInspect(fit, "options")$test
    if (isTRUE(topt == "none")) {
      # bu run'da fitMeasures mumkun degil, reject
      next
    }
    
    chk <- pass_fit_full(fit, include_measurement_p = include_measurement_p)
    
    if (isTRUE(chk$ok)) {
      acc <- acc + 1
      accepted[[acc]] <- dat
      
      fm <- chk$f_meas
      log_df <- rbind(log_df, data.frame(
        accepted_id = acc,
        attempt = attempt,
        pmax = chk$pmax,
        srmr = fm["srmr"],
        gfi  = fm["gfi"],
        agfi = fm["agfi"],
        cfi  = fm["cfi"],
        nfi  = fm["nfi"],
        nnfi = fm["nnfi"],
        ifi  = fm["ifi"],
        rmsea = fm["rmsea"],
        chisq = fm["chisq"],
        df = fm["df"],
        chisq_df = fm["chisq"]/fm["df"]
      ))
      
      if (acc %% verbose_every == 0) {
        cat(sprintf("Accepted %d/%d (attempt %d) | pmax=%.4f CFI=%.3f RMSEA=%.3f SRMR=%.3f\n",
                    acc, K, attempt, chk$pmax, fm["cfi"], fm["rmsea"], fm["srmr"]))
      }
    }
  }
  
  list(datasets = accepted, fitlog = log_df, n = n, K = K, total_attempts = attempt)
}


#run example
out_500 <- generate_sem_datasets_filtered(model_pop, model_fit, n = 500, K = 1000)

# Kac denemede 999+1 kabul topladi
out_500$total_attempts

# Kabul edilenlerin fit ozetleri
summary(out_500$fitlog[, c("pmax","cfi","rmsea","srmr","chisq_df")])

#aggregation
# out_500$datasets : list of data.frames (A1..E3)

aggregate_to_latent_means <- function(df) {
  data.frame(
    B = rowMeans(df[, c("B1","B2","B3")], na.rm = TRUE),
    D = rowMeans(df[, c("D1","D2","D3")], na.rm = TRUE),
    E = rowMeans(df[, c("E1","E2","E3")], na.rm = TRUE)
  )
}



# BN pipeline'a girecek liste (her biri A..E olan 500 gozlemlik data.frame)
data_list_bn_500 <- lapply(out_500$datasets, aggregate_to_latent_means)

#veri seti olusturuldu. artik algoritmalari deneyelim

# ground truth (Model3)
true_edges <- c(
  "E B",
  "B D"
)

nodes <- c("B","D","E")
all_edges <- as.vector(outer(nodes, nodes, FUN = paste))
all_edges <- all_edges[!grepl("(B B|D D|E E)", all_edges)]


#algorithm wrappers
library(bnlearn)

algo_fns <- list(
  "H2PC"   = function(d) h2pc(d),
  "MMHC"   = function(d) mmhc(d),
  "RSMAX2" = function(d) rsmax2(d),
  "HC"     = function(d) hc(d),
  "PCS"    = function(d) pc.stable(d),
  "Tabu"   = function(d) tabu(d),
  "HPC"    = function(d) hpc(d),
  "SPC"    = function(d) si.hiton.pc(d),
  "IAMB"   = function(d) iamb(d),
  "I-IAMB" = function(d) inter.iamb(d),
  "IAMB-F" = function(d) iamb.fdr(d),
  "F-IAMB" = function(d) fast.iamb(d),
  "MMPC"   = function(d) mmpc(d),
  "GS"     = function(d) gs(d)
)

#extract estimated edges

extract_edges <- function(bn_fit) {
  am <- amat(bn_fit)
  idx <- which(am == 1, arr.ind = TRUE)
  if (nrow(idx) == 0) return(character(0))
  paste(rownames(am)[idx[,1]], colnames(am)[idx[,2]])
}

#metrics for single dataset
compute_shd_directed <- function(pred, truth) {
  # pred, truth: character vector "X Y"
  
  all_nodes <- unique(unlist(strsplit(c(pred, truth), " ")))
  all_edges <- as.vector(outer(all_nodes, all_nodes, paste))
  all_edges <- all_edges[!grepl("^(.+) \\1$", all_edges)]
  
  FP <- sum(pred %in% setdiff(all_edges, truth))
  FN <- sum(truth %in% setdiff(all_edges, pred))
  
  FP + FN
}

compute_metrics <- function(pred, truth) {
  
  TP <- sum(pred %in% truth)
  FP <- sum(pred %in% setdiff(all_edges, truth))
  FN <- sum(!(truth %in% pred))
  TN <- length(all_edges) - TP - FP - FN
  
  ACC <- (TP + TN) / (TP + FP + FN + TN)
  SEN <- ifelse(TP+FN==0, NA, TP/(TP+FN))
  PRE <- ifelse(TP+FP==0, NA, TP/(TP+FP))
  SPE <- ifelse(TN+FP==0, NA, TN/(TN+FP))
  F1  <- ifelse(is.na(PRE+SEN), NA, 2*PRE*SEN/(PRE+SEN))
  
  Pc <- (((TP+FN)*(TP+FP))+((TN+FN)*(TN+FP)))/(TP+FP+FN+TN)^2
  CK <- (ACC - Pc)/(1 - Pc)
  
  MCC <- ifelse(
    (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)==0, NA,
    (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  )
  
  SHD <- compute_shd_directed(pred, truth)
  
  c(ACC,SEN,PRE,SPE,F1,CK,MCC,SHD)
}


#999+1x14  run
run_all_algorithms <- function(data_list) {
  
  metrics <- c("ACC","SEN","PRE","SPE","F1","CK","MCC","SHD")
  results <- list()
  
  for (alg in names(algo_fns)) {
    cat("Running:", alg, "\n")
    mat <- matrix(NA, nrow = length(data_list), ncol = length(metrics))
    
    for (i in seq_along(data_list)) {
      fit <- tryCatch(algo_fns[[alg]](data_list[[i]]), error = function(e) NULL)
      if (is.null(fit)) next
      
      pred <- extract_edges(fit)
      mat[i,] <- compute_metrics(pred, true_edges)
    }
    
    colnames(mat) <- metrics
    results[[alg]] <- as.data.frame(mat)
  }
  results
}

res_500 <- run_all_algorithms(data_list_bn_500)

#create final table
fmt_num <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_character_, length(x))
  
  ok <- !is.na(x)
  
  # |x| < 1  ??? 2 decimals, no leading zero
  small <- ok & abs(x) < 1
  out[small] <- sub(
    "^0",
    "",
    formatC(round(x[small], 2), format = "f", digits = 2)
  )
  
  # |x| ??? 1 ??? 1 decimal
  big <- ok & abs(x) >= 1
  out[big] <- formatC(round(x[big], 1), format = "f", digits = 1)
  
  out
}

make_summary_table <- function(res) {
  
  metrics <- c("ACC","SEN","PRE","SPE","F1","CK","MCC","SHD")
  out <- data.frame(Algorithm = names(res))
  
  for (m in metrics) {
    vals <- sapply(res, function(x) mean(x[[m]], na.rm=TRUE))
    sds  <- sapply(res, function(x) sd(x[[m]], na.rm=TRUE))
    if (m == "SHD") {
      ranks <- rank(vals, ties.method = "average")
    } else {
      ranks <- rank(-vals, ties.method = "average")
    }
    
    out[[m]] <- paste0(fmt_num(vals), " (", fmt_num(sds), ")")
    out[[paste0("R_",m)]] <- ranks
  }
  
  rank_cols <- grep("^R_", names(out))
  out$Overall_Rank <- rowMeans(out[, rank_cols])
  
  out[order(out$Overall_Rank), ]
}

final_table_500 <- make_summary_table(res_500)
final_table_500
library(writexl)
write_xlsx(final_table_500,"D:/Makaleler/SEM Simulasyon/Revs/Rev13 Array/Added Analyses/model_3_final_table_500.xlsx")

#1000 icin 

#gozlem sayisini degistir sadece
library(lavaan)
model_pop <- '
# Measurement model
B =~ 0.8*B1 + 0.8*B2 + 0.8*B3
D =~ 0.8*D1 + 0.8*D2 + 0.8*D3
E =~ 0.8*E1 + 0.8*E2 + 0.8*E3

# Structural model (Model3: E->B, B->D)
B ~ 0.6*E
D ~ 0.6*B

# (Optional) latent variances
B ~~ 1*B
D ~~ 1*D
E ~~ 1*E

# (Optional) indicator residual variances (since loading=0.8 => residual var=0.36)
B1 ~~ 0.36*B1
B2 ~~ 0.36*B2
B3 ~~ 0.36*B3
D1 ~~ 0.36*D1
D2 ~~ 0.36*D2
D3 ~~ 0.36*D3
E1 ~~ 0.36*E1
E2 ~~ 0.36*E2
E3 ~~ 0.36*E3
'

model_fit <- '
# Measurement model (free loadings)
B =~ B1 + B2 + B3
D =~ D1 + D2 + D3
E =~ E1 + E2 + E3

# Structural model (free regressions) - Model3
B ~ E
D ~ B
'
#fit control function

# 1) istersen "tum yapisal yollar anlamli mi" kontrol??:
#    - sadece regressions: op == "~"
#    - ister measurement ( "=~") da dahil edebilirsin (parametreden a??t??m)
max_pvalue_check <- function(fit, include_measurement = FALSE) {
  pe <- parameterEstimates(fit)
  ops_keep <- if (include_measurement) c("~", "=~") else c("~")
  pe <- pe[pe$op %in% ops_keep, ]
  # Sabit/variance gibi seyleri karsilastirmamak icin:
  pe <- pe[is.finite(pe$pvalue), ]
  if (nrow(pe) == 0) return(NA_real_)
  max(pe$pvalue, na.rm = TRUE)
}

# 2) Fit filtresi (senin kosullarin; 0.95 -> 0.90 guncellendi)
pass_fit_full <- function(fit,
                          p_max = 0.05,
                          srmr_max = 0.05,
                          rmsea_max = 0.05,
                          min_fit = 0.90,
                          chisq_df_max = 3,
                          include_measurement_p = FALSE) {
  
  # Eger fit test="none" ile uretildiyse fitMeasures alinamaz -> direkt reject
  if (isTRUE(lavInspect(fit, "options")$test == "none")) {
    return(list(ok = FALSE, f_meas = NA, pmax = NA_real_))
  }
  
  f_meas <- fitMeasures(fit, c("srmr","gfi","agfi","cfi","nfi","nnfi","ifi","rmsea","chisq","df"))
  pmax <- max_pvalue_check(fit, include_measurement = include_measurement_p)
  
  ok <- is.finite(pmax) && (pmax <= p_max) &&
    is.finite(f_meas["srmr"]) && (f_meas["srmr"] <= srmr_max) &&
    is.finite(f_meas["rmsea"]) && (f_meas["rmsea"] <= rmsea_max) &&
    all(is.finite(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")])) &&
    all(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")] >= min_fit) &&
    all(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")] <= 1) &&
    is.finite(f_meas["chisq"]) && is.finite(f_meas["df"]) &&
    (f_meas["chisq"] / f_meas["df"] < chisq_df_max)
  
  list(ok = ok, f_meas = f_meas, pmax = pmax)
}

#produce - fit SEM - if yes save

generate_sem_datasets_filtered <- function(model_pop, model_fit, n, K = 1000,
                                           max_attempts = 200000,
                                           seed = NULL,
                                           verbose_every = 50,
                                           include_measurement_p = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  accepted <- vector("list", K)
  log_df <- data.frame(
    accepted_id = integer(0),
    attempt = integer(0),
    pmax = numeric(0),
    srmr = numeric(0),
    gfi = numeric(0),
    agfi = numeric(0),
    cfi = numeric(0),
    nfi = numeric(0),
    nnfi = numeric(0),
    ifi = numeric(0),
    rmsea = numeric(0),
    chisq = numeric(0),
    df = numeric(0),
    chisq_df = numeric(0)
  )
  
  acc <- 0
  attempt <- 0
  
  while (acc < K) {
    attempt <- attempt + 1
    if (attempt > max_attempts) stop("max_attempts asildi. Esikleri gevset veya max_attempts artir.")
    
    dat <- tryCatch(simulateData(model_pop, sample.nobs = n), error = function(e) NULL)
    if (is.null(dat)) next
    
    fit <- tryCatch(
      lavaan::sem(model_fit, data = dat, std.lv = TRUE, test = "standard"),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    # <-- kritik kontrol
    topt <- lavInspect(fit, "options")$test
    if (isTRUE(topt == "none")) {
      # bu run'da fitMeasures mumkun degil, reject
      next
    }
    
    chk <- pass_fit_full(fit, include_measurement_p = include_measurement_p)
    
    if (isTRUE(chk$ok)) {
      acc <- acc + 1
      accepted[[acc]] <- dat
      
      fm <- chk$f_meas
      log_df <- rbind(log_df, data.frame(
        accepted_id = acc,
        attempt = attempt,
        pmax = chk$pmax,
        srmr = fm["srmr"],
        gfi  = fm["gfi"],
        agfi = fm["agfi"],
        cfi  = fm["cfi"],
        nfi  = fm["nfi"],
        nnfi = fm["nnfi"],
        ifi  = fm["ifi"],
        rmsea = fm["rmsea"],
        chisq = fm["chisq"],
        df = fm["df"],
        chisq_df = fm["chisq"]/fm["df"]
      ))
      
      if (acc %% verbose_every == 0) {
        cat(sprintf("Accepted %d/%d (attempt %d) | pmax=%.4f CFI=%.3f RMSEA=%.3f SRMR=%.3f\n",
                    acc, K, attempt, chk$pmax, fm["cfi"], fm["rmsea"], fm["srmr"]))
      }
    }
  }
  
  list(datasets = accepted, fitlog = log_df, n = n, K = K, total_attempts = attempt)
}


#run example
out_1000 <- generate_sem_datasets_filtered(model_pop, model_fit, n = 1000, K = 1000)

# Kac denemede 999+1 kabul topladi
out_1000$total_attempts

# Kabul edilenlerin fit ozetleri
summary(out_1000$fitlog[, c("pmax","cfi","rmsea","srmr","chisq_df")])

#aggregation
# out_1000$datasets : list of data.frames (A1..E3)

aggregate_to_latent_means <- function(df) {
  data.frame(
    B = rowMeans(df[, c("B1","B2","B3")], na.rm = TRUE),
    D = rowMeans(df[, c("D1","D2","D3")], na.rm = TRUE),
    E = rowMeans(df[, c("E1","E2","E3")], na.rm = TRUE)
  )
}



# BN pipeline'a girecek liste (her biri A..E olan 1000 gozlemlik data.frame)
data_list_bn_1000 <- lapply(out_1000$datasets, aggregate_to_latent_means)

#veri seti olusturuldu. artik algoritmalari deneyelim

# ground truth (Model3)
true_edges <- c(
  "E B",
  "B D"
)

nodes <- c("B","D","E")
all_edges <- as.vector(outer(nodes, nodes, FUN = paste))
all_edges <- all_edges[!grepl("(B B|D D|E E)", all_edges)]


#algorithm wrappers
library(bnlearn)

algo_fns <- list(
  "H2PC"   = function(d) h2pc(d),
  "MMHC"   = function(d) mmhc(d),
  "RSMAX2" = function(d) rsmax2(d),
  "HC"     = function(d) hc(d),
  "PCS"    = function(d) pc.stable(d),
  "Tabu"   = function(d) tabu(d),
  "HPC"    = function(d) hpc(d),
  "SPC"    = function(d) si.hiton.pc(d),
  "IAMB"   = function(d) iamb(d),
  "I-IAMB" = function(d) inter.iamb(d),
  "IAMB-F" = function(d) iamb.fdr(d),
  "F-IAMB" = function(d) fast.iamb(d),
  "MMPC"   = function(d) mmpc(d),
  "GS"     = function(d) gs(d)
)

#extract estimated edges

extract_edges <- function(bn_fit) {
  am <- amat(bn_fit)
  idx <- which(am == 1, arr.ind = TRUE)
  if (nrow(idx) == 0) return(character(0))
  paste(rownames(am)[idx[,1]], colnames(am)[idx[,2]])
}

#metrics for single dataset
compute_shd_directed <- function(pred, truth) {
  # pred, truth: character vector "X Y"
  
  all_nodes <- unique(unlist(strsplit(c(pred, truth), " ")))
  all_edges <- as.vector(outer(all_nodes, all_nodes, paste))
  all_edges <- all_edges[!grepl("^(.+) \\1$", all_edges)]
  
  FP <- sum(pred %in% setdiff(all_edges, truth))
  FN <- sum(truth %in% setdiff(all_edges, pred))
  
  FP + FN
}

compute_metrics <- function(pred, truth) {
  
  TP <- sum(pred %in% truth)
  FP <- sum(pred %in% setdiff(all_edges, truth))
  FN <- sum(!(truth %in% pred))
  TN <- length(all_edges) - TP - FP - FN
  
  ACC <- (TP + TN) / (TP + FP + FN + TN)
  SEN <- ifelse(TP+FN==0, NA, TP/(TP+FN))
  PRE <- ifelse(TP+FP==0, NA, TP/(TP+FP))
  SPE <- ifelse(TN+FP==0, NA, TN/(TN+FP))
  F1  <- ifelse(is.na(PRE+SEN), NA, 2*PRE*SEN/(PRE+SEN))
  
  Pc <- (((TP+FN)*(TP+FP))+((TN+FN)*(TN+FP)))/(TP+FP+FN+TN)^2
  CK <- (ACC - Pc)/(1 - Pc)
  
  MCC <- ifelse(
    (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)==0, NA,
    (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  )
  
  SHD <- compute_shd_directed(pred, truth)
  
  c(ACC,SEN,PRE,SPE,F1,CK,MCC,SHD)
}


#999+1x14  run
run_all_algorithms <- function(data_list) {
  
  metrics <- c("ACC","SEN","PRE","SPE","F1","CK","MCC","SHD")
  results <- list()
  
  for (alg in names(algo_fns)) {
    cat("Running:", alg, "\n")
    mat <- matrix(NA, nrow = length(data_list), ncol = length(metrics))
    
    for (i in seq_along(data_list)) {
      fit <- tryCatch(algo_fns[[alg]](data_list[[i]]), error = function(e) NULL)
      if (is.null(fit)) next
      
      pred <- extract_edges(fit)
      mat[i,] <- compute_metrics(pred, true_edges)
    }
    
    colnames(mat) <- metrics
    results[[alg]] <- as.data.frame(mat)
  }
  results
}

res_1000 <- run_all_algorithms(data_list_bn_1000)

#create final table
fmt_num <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_character_, length(x))
  
  ok <- !is.na(x)
  
  # |x| < 1  ??? 2 decimals, no leading zero
  small <- ok & abs(x) < 1
  out[small] <- sub(
    "^0",
    "",
    formatC(round(x[small], 2), format = "f", digits = 2)
  )
  
  # |x| ??? 1 ??? 1 decimal
  big <- ok & abs(x) >= 1
  out[big] <- formatC(round(x[big], 1), format = "f", digits = 1)
  
  out
}

make_summary_table <- function(res) {
  
  metrics <- c("ACC","SEN","PRE","SPE","F1","CK","MCC","SHD")
  out <- data.frame(Algorithm = names(res))
  
  for (m in metrics) {
    vals <- sapply(res, function(x) mean(x[[m]], na.rm=TRUE))
    sds  <- sapply(res, function(x) sd(x[[m]], na.rm=TRUE))
    if (m == "SHD") {
      ranks <- rank(vals, ties.method = "average")
    } else {
      ranks <- rank(-vals, ties.method = "average")
    }
    
    out[[m]] <- paste0(fmt_num(vals), " (", fmt_num(sds), ")")
    out[[paste0("R_",m)]] <- ranks
  }
  
  rank_cols <- grep("^R_", names(out))
  out$Overall_Rank <- rowMeans(out[, rank_cols])
  
  out[order(out$Overall_Rank), ]
}

final_table_1000 <- make_summary_table(res_1000)
final_table_1000
library(writexl)
write_xlsx(final_table_1000,"D:/Makaleler/SEM Simulasyon/Revs/Rev13 Array/Added Analyses/model_3_final_table_1000.xlsx")


#2500 icin

#gozlem sayisini degistir sadece
library(lavaan)
model_pop <- '
# Measurement model
B =~ 0.8*B1 + 0.8*B2 + 0.8*B3
D =~ 0.8*D1 + 0.8*D2 + 0.8*D3
E =~ 0.8*E1 + 0.8*E2 + 0.8*E3

# Structural model (Model3: E->B, B->D)
B ~ 0.6*E
D ~ 0.6*B

# (Optional) latent variances
B ~~ 1*B
D ~~ 1*D
E ~~ 1*E

# (Optional) indicator residual variances (since loading=0.8 => residual var=0.36)
B1 ~~ 0.36*B1
B2 ~~ 0.36*B2
B3 ~~ 0.36*B3
D1 ~~ 0.36*D1
D2 ~~ 0.36*D2
D3 ~~ 0.36*D3
E1 ~~ 0.36*E1
E2 ~~ 0.36*E2
E3 ~~ 0.36*E3
'

model_fit <- '
# Measurement model (free loadings)
B =~ B1 + B2 + B3
D =~ D1 + D2 + D3
E =~ E1 + E2 + E3

# Structural model (free regressions) - Model3
B ~ E
D ~ B
'
#fit control function

# 1) istersen "tum yapisal yollar anlamli mi" kontrol??:
#    - sadece regressions: op == "~"
#    - ister measurement ( "=~") da dahil edebilirsin (parametreden a??t??m)
max_pvalue_check <- function(fit, include_measurement = FALSE) {
  pe <- parameterEstimates(fit)
  ops_keep <- if (include_measurement) c("~", "=~") else c("~")
  pe <- pe[pe$op %in% ops_keep, ]
  # Sabit/variance gibi seyleri karsilastirmamak icin:
  pe <- pe[is.finite(pe$pvalue), ]
  if (nrow(pe) == 0) return(NA_real_)
  max(pe$pvalue, na.rm = TRUE)
}

# 2) Fit filtresi (senin kosullarin; 0.95 -> 0.90 guncellendi)
pass_fit_full <- function(fit,
                          p_max = 0.05,
                          srmr_max = 0.05,
                          rmsea_max = 0.05,
                          min_fit = 0.90,
                          chisq_df_max = 3,
                          include_measurement_p = FALSE) {
  
  # Eger fit test="none" ile uretildiyse fitMeasures alinamaz -> direkt reject
  if (isTRUE(lavInspect(fit, "options")$test == "none")) {
    return(list(ok = FALSE, f_meas = NA, pmax = NA_real_))
  }
  
  f_meas <- fitMeasures(fit, c("srmr","gfi","agfi","cfi","nfi","nnfi","ifi","rmsea","chisq","df"))
  pmax <- max_pvalue_check(fit, include_measurement = include_measurement_p)
  
  ok <- is.finite(pmax) && (pmax <= p_max) &&
    is.finite(f_meas["srmr"]) && (f_meas["srmr"] <= srmr_max) &&
    is.finite(f_meas["rmsea"]) && (f_meas["rmsea"] <= rmsea_max) &&
    all(is.finite(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")])) &&
    all(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")] >= min_fit) &&
    all(f_meas[c("gfi","agfi","cfi","nfi","nnfi","ifi")] <= 1) &&
    is.finite(f_meas["chisq"]) && is.finite(f_meas["df"]) &&
    (f_meas["chisq"] / f_meas["df"] < chisq_df_max)
  
  list(ok = ok, f_meas = f_meas, pmax = pmax)
}

#produce - fit SEM - if yes save

generate_sem_datasets_filtered <- function(model_pop, model_fit, n, K = 1000,
                                           max_attempts = 200000,
                                           seed = NULL,
                                           verbose_every = 50,
                                           include_measurement_p = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  accepted <- vector("list", K)
  log_df <- data.frame(
    accepted_id = integer(0),
    attempt = integer(0),
    pmax = numeric(0),
    srmr = numeric(0),
    gfi = numeric(0),
    agfi = numeric(0),
    cfi = numeric(0),
    nfi = numeric(0),
    nnfi = numeric(0),
    ifi = numeric(0),
    rmsea = numeric(0),
    chisq = numeric(0),
    df = numeric(0),
    chisq_df = numeric(0)
  )
  
  acc <- 0
  attempt <- 0
  
  while (acc < K) {
    attempt <- attempt + 1
    if (attempt > max_attempts) stop("max_attempts asildi. Esikleri gevset veya max_attempts artir.")
    
    dat <- tryCatch(simulateData(model_pop, sample.nobs = n), error = function(e) NULL)
    if (is.null(dat)) next
    
    fit <- tryCatch(
      lavaan::sem(model_fit, data = dat, std.lv = TRUE, test = "standard"),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    # <-- kritik kontrol
    topt <- lavInspect(fit, "options")$test
    if (isTRUE(topt == "none")) {
      # bu run'da fitMeasures mumkun degil, reject
      next
    }
    
    chk <- pass_fit_full(fit, include_measurement_p = include_measurement_p)
    
    if (isTRUE(chk$ok)) {
      acc <- acc + 1
      accepted[[acc]] <- dat
      
      fm <- chk$f_meas
      log_df <- rbind(log_df, data.frame(
        accepted_id = acc,
        attempt = attempt,
        pmax = chk$pmax,
        srmr = fm["srmr"],
        gfi  = fm["gfi"],
        agfi = fm["agfi"],
        cfi  = fm["cfi"],
        nfi  = fm["nfi"],
        nnfi = fm["nnfi"],
        ifi  = fm["ifi"],
        rmsea = fm["rmsea"],
        chisq = fm["chisq"],
        df = fm["df"],
        chisq_df = fm["chisq"]/fm["df"]
      ))
      
      if (acc %% verbose_every == 0) {
        cat(sprintf("Accepted %d/%d (attempt %d) | pmax=%.4f CFI=%.3f RMSEA=%.3f SRMR=%.3f\n",
                    acc, K, attempt, chk$pmax, fm["cfi"], fm["rmsea"], fm["srmr"]))
      }
    }
  }
  
  list(datasets = accepted, fitlog = log_df, n = n, K = K, total_attempts = attempt)
}


#run example
out_2500 <- generate_sem_datasets_filtered(model_pop, model_fit, n = 2500, K = 1000)

# Kac denemede 999+1 kabul topladi
out_2500$total_attempts

# Kabul edilenlerin fit ozetleri
summary(out_2500$fitlog[, c("pmax","cfi","rmsea","srmr","chisq_df")])

#aggregation
# out_2500$datasets : list of data.frames (A1..E3)

aggregate_to_latent_means <- function(df) {
  data.frame(
    B = rowMeans(df[, c("B1","B2","B3")], na.rm = TRUE),
    D = rowMeans(df[, c("D1","D2","D3")], na.rm = TRUE),
    E = rowMeans(df[, c("E1","E2","E3")], na.rm = TRUE)
  )
}



# BN pipeline'a girecek liste (her biri A..E olan 2500 gozlemlik data.frame)
data_list_bn_2500 <- lapply(out_2500$datasets, aggregate_to_latent_means)

#veri seti olusturuldu. artik algoritmalari deneyelim

# ground truth (Model3)
true_edges <- c(
  "E B",
  "B D"
)

nodes <- c("B","D","E")
all_edges <- as.vector(outer(nodes, nodes, FUN = paste))
all_edges <- all_edges[!grepl("(B B|D D|E E)", all_edges)]


#algorithm wrappers
library(bnlearn)

algo_fns <- list(
  "H2PC"   = function(d) h2pc(d),
  "MMHC"   = function(d) mmhc(d),
  "RSMAX2" = function(d) rsmax2(d),
  "HC"     = function(d) hc(d),
  "PCS"    = function(d) pc.stable(d),
  "Tabu"   = function(d) tabu(d),
  "HPC"    = function(d) hpc(d),
  "SPC"    = function(d) si.hiton.pc(d),
  "IAMB"   = function(d) iamb(d),
  "I-IAMB" = function(d) inter.iamb(d),
  "IAMB-F" = function(d) iamb.fdr(d),
  "F-IAMB" = function(d) fast.iamb(d),
  "MMPC"   = function(d) mmpc(d),
  "GS"     = function(d) gs(d)
)

#extract estimated edges

extract_edges <- function(bn_fit) {
  am <- amat(bn_fit)
  idx <- which(am == 1, arr.ind = TRUE)
  if (nrow(idx) == 0) return(character(0))
  paste(rownames(am)[idx[,1]], colnames(am)[idx[,2]])
}

#metrics for single dataset
compute_shd_directed <- function(pred, truth) {
  # pred, truth: character vector "X Y"
  
  all_nodes <- unique(unlist(strsplit(c(pred, truth), " ")))
  all_edges <- as.vector(outer(all_nodes, all_nodes, paste))
  all_edges <- all_edges[!grepl("^(.+) \\1$", all_edges)]
  
  FP <- sum(pred %in% setdiff(all_edges, truth))
  FN <- sum(truth %in% setdiff(all_edges, pred))
  
  FP + FN
}

compute_metrics <- function(pred, truth) {
  
  TP <- sum(pred %in% truth)
  FP <- sum(pred %in% setdiff(all_edges, truth))
  FN <- sum(!(truth %in% pred))
  TN <- length(all_edges) - TP - FP - FN
  
  ACC <- (TP + TN) / (TP + FP + FN + TN)
  SEN <- ifelse(TP+FN==0, NA, TP/(TP+FN))
  PRE <- ifelse(TP+FP==0, NA, TP/(TP+FP))
  SPE <- ifelse(TN+FP==0, NA, TN/(TN+FP))
  F1  <- ifelse(is.na(PRE+SEN), NA, 2*PRE*SEN/(PRE+SEN))
  
  Pc <- (((TP+FN)*(TP+FP))+((TN+FN)*(TN+FP)))/(TP+FP+FN+TN)^2
  CK <- (ACC - Pc)/(1 - Pc)
  
  MCC <- ifelse(
    (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)==0, NA,
    (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  )
  
  SHD <- compute_shd_directed(pred, truth)
  
  c(ACC,SEN,PRE,SPE,F1,CK,MCC,SHD)
}


#999+1x14  run
run_all_algorithms <- function(data_list) {
  
  metrics <- c("ACC","SEN","PRE","SPE","F1","CK","MCC","SHD")
  results <- list()
  
  for (alg in names(algo_fns)) {
    cat("Running:", alg, "\n")
    mat <- matrix(NA, nrow = length(data_list), ncol = length(metrics))
    
    for (i in seq_along(data_list)) {
      fit <- tryCatch(algo_fns[[alg]](data_list[[i]]), error = function(e) NULL)
      if (is.null(fit)) next
      
      pred <- extract_edges(fit)
      mat[i,] <- compute_metrics(pred, true_edges)
    }
    
    colnames(mat) <- metrics
    results[[alg]] <- as.data.frame(mat)
  }
  results
}

res_2500 <- run_all_algorithms(data_list_bn_2500)

#create final table
fmt_num <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_character_, length(x))
  
  ok <- !is.na(x)
  
  # |x| < 1  ??? 2 decimals, no leading zero
  small <- ok & abs(x) < 1
  out[small] <- sub(
    "^0",
    "",
    formatC(round(x[small], 2), format = "f", digits = 2)
  )
  
  # |x| ??? 1 ??? 1 decimal
  big <- ok & abs(x) >= 1
  out[big] <- formatC(round(x[big], 1), format = "f", digits = 1)
  
  out
}

make_summary_table <- function(res) {
  
  metrics <- c("ACC","SEN","PRE","SPE","F1","CK","MCC","SHD")
  out <- data.frame(Algorithm = names(res))
  
  for (m in metrics) {
    vals <- sapply(res, function(x) mean(x[[m]], na.rm=TRUE))
    sds  <- sapply(res, function(x) sd(x[[m]], na.rm=TRUE))
    if (m == "SHD") {
      ranks <- rank(vals, ties.method = "average")
    } else {
      ranks <- rank(-vals, ties.method = "average")
    }
    
    out[[m]] <- paste0(fmt_num(vals), " (", fmt_num(sds), ")")
    out[[paste0("R_",m)]] <- ranks
  }
  
  rank_cols <- grep("^R_", names(out))
  out$Overall_Rank <- rowMeans(out[, rank_cols])
  
  out[order(out$Overall_Rank), ]
}

final_table_2500 <- make_summary_table(res_2500)
final_table_2500
library(writexl)
write_xlsx(final_table_2500,"D:/Makaleler/SEM Simulasyon/Revs/Rev13 Array/Added Analyses/model_3_final_table_2500.xlsx")
bit <- Sys.time()
bit - bas
library(beepr)
beep(sound = 8)
