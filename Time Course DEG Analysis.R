#Load required packages
library(maSigPro)
library(readxl)
library(tidyverse)
library(glm2)
library(mclust)

stepfor2 <- function (y = y, d = d, alfa = 0.05, family = gaussian(), epsilon = 1e-05) 
{
  pval <- NULL
  design <- NULL
  j = 1
  resul0 <- summary(glm2(y ~ ., data = d, family = family, epsilon = epsilon))$coefficients[, 
                                                                                           4]
  d <- as.data.frame(d[, names(resul0)[-1]])
  for (i in 1:ncol(d)) {
    sub <- cbind(design, d[, i])
    sub <- as.data.frame(sub)
    lm2 <- glm2(y ~ ., data = sub, family = family, epsilon = epsilon)
    result <- summary(lm2)
    pval[i] <- result$coefficients[, 4][j + 1]
  }
  min <- min(pval)
  while (min < alfa) {
    b <- pval == min
    c <- c(1:length(pval))
    pos <- c[b]
    pos <- pos[!is.na(pos)][1]
    design <- cbind(design, d[, pos])
    design <- as.data.frame(design)
    colnames(design)[j] <- colnames(d)[pos]
    j = j + 1
    if (ncol(d) == 2) {
      lastname <- colnames(d)[!b]
    }
    d <- as.data.frame(d[, -pos])
    if (ncol(d) == 1) {
      colnames(d) = lastname
    }
    pval <- NULL
    if (ncol(d) != 0) {
      for (i in 1:ncol(d)) {
        sub <- cbind(design, d[, i])
        sub <- as.data.frame(sub)
        lm2 <- glm2(y ~ ., data = sub, family = family, 
                   epsilon = epsilon)
        result <- summary(lm2)
        pval[i] <- result$coefficients[, 4][j + 1]
      }
      min <- min(pval, na.rm = TRUE)
    }
    else min <- 1
  }
  if (is.null(design)) {
    lm1 <- glm2(y ~ 1, family = family, epsilon = epsilon)
  }
  else {
    lm1 <- glm2(y ~ ., data = design, family = family, epsilon = epsilon)
  }
  return(lm1)
}

T.fit2 <- function (data, design = data$dis, step.method = "backward", 
                    min.obs = data$min.obs, alfa = data$Q, nvar.correction = FALSE, 
                    family = gaussian(), epsilon = 1e-05, item = "gene") 
{
  if (is.list(data)) {
    dat <- as.matrix(data$SELEC)
    dat <- rbind(c(rep(1, ncol(dat))), dat)
    groups.vector <- data$groups.vector
    groups.vector <- c(groups.vector[nchar(groups.vector) == 
                                       min(nchar(groups.vector))][1], groups.vector)
    edesign <- data$edesign
    G <- data$g
    family <- data$family
  }
  else {
    G <- nrow(data)
    data <- rbind(c(rep(1, ncol(data))), data)
    dat <- as.matrix(data)
    count.na <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[apply(dat, 1, count.na) >= min.obs, ]
    groups.vector = NULL
    edesign = NULL
  }
  dis <- as.data.frame(design)
  dat <- dat[, as.character(rownames(dis))]
  g <- (dim(dat)[1] - 1)
  n <- dim(dat)[2]
  p <- dim(dis)[2]
  vars.in <- colnames(dis)
  sol <- coefficients <- group.coeffs <- t.score <- sig.profiles <- NULL
  influ.info <- matrix(NA, nrow = nrow(dis), ncol = 1)
  rownames(influ.info) <- rownames(dis)
  if (nvar.correction) 
    alfa <- alfa/ncol(dis)
  for (i in 2:(g + 1)) {
    y <- as.numeric(dat[i, ])
    name <- rownames(dat)[i]
    if (step.method == "backward") {
      reg <- stepback(y = y, d = dis, alfa = alfa, family = family, 
                      epsilon = epsilon)
    }
    else if (step.method == "forward") {
      reg <- stepfor2(y = y, d = dis, alfa = alfa, family = family, 
                     epsilon = epsilon)
    }
    else if (step.method == "two.ways.backward") {
      reg <- two.ways.stepback(y = y, d = dis, alfa = alfa, 
                               family = family, epsilon = epsilon)
    }
    else if (step.method == "two.ways.forward") {
      reg <- two.ways.stepfor(y = y, d = dis, alfa = alfa, 
                              family = family, epsilon = epsilon)
    }
    else stop("stepwise method must be one of backward, forward, two.ways.backward, two.ways.forward")
    div <- c(1:round(g/100)) * 100
    if (is.element(i, div)) 
      print(paste(c("fitting ", item, i, "out of", g), 
                  collapse = " "))
    lmf <- glm2(y ~ ., data = as.data.frame(dis), family = family, 
                epsilon = epsilon)
    result <- summary(lmf)
    novar <- vars.in[!is.element(vars.in, names(result$coefficients[, 
                                                                    4]))]
    influ <- influence.measures(reg)$is.inf
    influ <- influ[, c(ncol(influ) - 3, ncol(influ) - 1)]
    influ1 <- which(apply(influ, 1, all))
    if (length(influ1) != 0) {
      paste.names <- function(a) {
        paste(names(a)[a], collapse = "/")
      }
      match <- match(rownames(dis), rownames(influ))
      influ <- as.data.frame(apply(influ, 1, paste.names))
      influ.info <- cbind(influ.info, influ[match, ])
      colnames(influ.info)[ncol(influ.info)] <- name
    }
    result <- summary(reg)
    if ((!(result$aic == -Inf) & !is.na(result$aic) & family$family == 
         "gaussian") | family$family != "gaussian") {
      k <- i
      model.glm.0 <- glm2(y ~ 1, family = family, epsilon = epsilon)
      if (family$family == "gaussian") {
        test <- anova(model.glm.0, reg, test = "F")
        p.value = test[6][2, 1]
      }
      else {
        test <- anova(model.glm.0, reg, test = "Chisq")
        p.value = test[5][2, 1]
      }
      bondad <- (reg$null.deviance - reg$deviance)/reg$null.deviance
      if (bondad < 0) {
        bondad = 0
      }
      beta.coeff <- result$coefficients[, 1]
      beta.p.valor <- result$coefficients[, 4]
      coeff <- rep(0, (length(vars.in) + 1))
      if (length(novar) != 0) {
        for (m in 1:length(novar)) {
          coeff[position(dis, novar[m]) + 1] <- NA
        }
      }
      p.valor <- t <- as.numeric(rep(NA, (length(vars.in) + 
                                            1)))
      if (result$coefficients[, 4][rownames(result$coefficients) == 
                                   "(Intercept)"] < alfa) {
        coeff[1] <- result$coefficients[, 1][rownames(result$coefficients) == 
                                               "(Intercept)"]
        p.valor[1] <- result$coefficients[, 4][rownames(result$coefficients) == 
                                                 "(Intercept)"]
        t[1] <- result$coefficients[, 3][rownames(result$coefficients) == 
                                           "(Intercept)"]
      }
      for (j in 2:length(coeff)) {
        if (is.element(vars.in[j - 1], rownames(result$coefficients))) {
          coeff[j] <- result$coefficients[, 1][rownames(result$coefficients) == 
                                                 vars.in[j - 1]]
          p.valor[j] <- result$coefficients[, 4][rownames(result$coefficients) == 
                                                   vars.in[j - 1]]
          t[j] <- result$coefficients[, 3][rownames(result$coefficients) == 
                                             vars.in[j - 1]]
        }
      }
      if (!all(is.na(p.valor))) {
        sol <- rbind(sol, as.numeric(c(p.value, bondad, 
                                       p.valor)))
        coefficients <- rbind(coefficients, coeff)
        t.score <- rbind(t.score, t)
        sig.profiles <- rbind(sig.profiles, y)
        h <- nrow(sol)
        rownames(sol)[h] <- name
        rownames(coefficients)[h] <- name
        rownames(t.score)[h] <- name
        rownames(sig.profiles)[h] <- name
      }
    }
  }
  if (!is.null(sol)) {
    sol <- as.data.frame(sol)
    coefficients <- as.data.frame(coefficients)
    coeffic <- coefficients
    t.score <- as.data.frame(t.score)
    sig.profiles <- as.data.frame(sig.profiles)
    colnames(sol) <- c("p-value", "R-squared", "p.valor_beta0", 
                       paste("p.valor_", vars.in, sep = ""))
    colnames(coefficients) <- c("beta0", paste("beta", vars.in, 
                                               sep = ""))
    colnames(t.score) <- c("t.score_beta0", paste("t.score_", 
                                                  vars.in, sep = ""))
    colnames(sig.profiles) <- colnames(dat)
    if (!is.null(groups.vector) & !is.null(edesign)) {
      groups <- colnames(edesign)[3:ncol(edesign)]
      degree <- (length(groups.vector)/length(groups)) - 
        1
      for (w in 1:nrow(coefficients)) {
        A <- NULL
        col.names <- NULL
        for (l in 1:length(groups)) {
          B <- reg.coeffs(coefficients = coefficients[w, 
          ], groups.vector = groups.vector, group = groups[l])
          cols <- paste(rep(groups[l], each = length(B)), 
                        paste("beta", c(0:(length(B) - 1)), sep = ""), 
                        sep = "_")
          A <- c(A, B)
          col.names <- c(col.names, cols)
        }
        group.coeffs <- (rbind(group.coeffs, A))
      }
      colnames(group.coeffs) <- col.names
      rownames(group.coeffs) <- rownames(coefficients)
    }
  }
  if (ncol(influ.info) > 2) {
    print(paste("Influence:", ncol(influ.info) - 1, "genes with influential data at slot influ.info. Model validation for these genes is recommended"))
  }
  influ.info <- influ.info[, -1]
  output <- list(sol, sig.profiles, coefficients, as.data.frame(group.coeffs), 
                 t.score, vars.in, G, g, dat, dis, step.method, groups.vector, 
                 edesign, influ.info)
  names(output) <- c("sol", "sig.profiles", "coefficients", 
                     "group.coeffs", "t.score", "variables", "G", "g", "dat", 
                     "dis", "step.method", "groups.vector", "edesign", "influ.info")
  output
}


#Initialize design matrix 
edesign = read_excel('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/time_design_all_conditions.xlsx')
edesign = column_to_rownames(edesign, var = "Client Sample ID")

#Initialize normalized data and reorder columns in normalized data to match design matrix
normalized_counts = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/normalized_counts.txt')
rownames(normalized_counts) = normalized_counts$X
normalized_counts$X = NULL

names(normalized_counts) <- gsub("\\.", "-", names(normalized_counts))
col_order <- c("A2-1", "A2-2", "A2-3", "A5-1", "A5-2", "A5-3", "A7-1", "A7-2", "A7-3", "A10-1","A10-2", "A10-3",
               "A14-1", "A14-2", "A14-3", "C2-1", "C2-2", "C2-3", "C5-1", "C5-2", "C5-3", "C7-1", "C7-2", "C7-3",
               "C10-1", "C10-2", "C10-3", "C14-1", "C14-2", "C14-3", "D2-1", "D2-2", "D2-3", "D5-1", "D5-2", "D5-3",
               "D7-1", "D7-2", "D7-3", "D10-1", "D10-2", "D10-3", "D14-1", "D14-2", "D14-3", "E2-1", "E2-2", "E2-3",
               "E5-1", "E5-2", "E5-3", "E7-1", "E7-2", "E7-3", "E10-1", "E10-2", "E10-3", "E14-1", "E14-2", "E14-3",
               "F2-1", "F2-2", "F2-3", "F5-1", "F5-2", "F5-3", "F7-1", "F7-2", "F7-3",
               "F10-1", "F10-2", "F10-3", "F14-1", "F14-2", "F14-3", "G2-1", "G2-2", "G2-3", "G5-1", "G5-2", 
               "G5-3", "G7-1", "G7-2", "G7-3", "G10-1", "G10-2", "G10-3", "G14-1","G14-2", "G14-3")
normalized_counts <- normalized_counts[,col_order]

#Time course expression analysis
design = make.design.matrix(edesign=edesign, degree=4)

fit <- p.vector(data=normalized_counts, design, counts=TRUE)
fit$i # returns the number of significant genes
SELEC_DF <- fit$SELEC # is a matrix with the significant genes and their expression values

tstep <- T.fit2(fit, step.method = "forward", alfa = 0.05)
#saveRDS("tstep", file = "/Users/ven/Documents/PhD/Rotation Botas Lab/Results/tstep_object_HD128vs200vsControl_counts=true")

sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
names(sigs$sig.genes)

sigs2 <- get.siggenes(tstep, rsq = 0.6, vars = "all")

sigs3 <- get.siggenes(tstep, rsq = 0.6, vars = "each")

##Generate and visualize results
#Visualize a Venn diagram and create a list of completely intersected genes between groups
names(sigs$summary)
VennDataFrame = suma2Venn(sigs$summary[, c(2:6)])
full_intersection_list = attr(VennDataFrame,"intersections")$`AD_B42vsNeg_Driver:AD_TauvsNeg_Driver:PDvsNeg_Driver:HD_128QvsNeg_Driver:HD_200QvsNeg_Driver`
write.table(full_intersection_list, file="/Users/ven/Documents/PhD/Rotation Botas Lab/Data/full_intersection_list.txt", sep="\t", quote=F, col.names=NA)

intersectionData = normalized_counts[FALSE,]

for (d in full_intersection_list){
  x <- subset(normalized_counts, rownames(normalized_counts) %in% d)
  intersectionData <- rbind (intersectionData, x)
}

#Visualize the consistency of the clusters and show clearly the differences between groups
results = see.genes(intersectionData, edesign = edesign, show.fit = F, dis =design$dis,
          cluster.method="Mclust" ,cluster.data = 1, k = 4)

#Export results (Venn Diagram, cluster consistency, and group differences) as PDF
for(d in dev.list()) {
  dev.set(d)
  Name = paste("/Users/ven/Documents/PhD/Rotation Botas Lab/Results/Graph", d, ".pdf", sep="")
  dev.copy(pdf, Name)
  dev.off()
}

textme()
