# The partial function of R package(propeller)
# propeller: finding statistically significant differences in cell type populations in single cell data
library(limma)
propeller = function (x = NULL, clusters = NULL, sample = NULL, group = NULL, 
                      trend = FALSE, robust = TRUE, transform = "logit") 
{
  if (is.null(x) & is.null(sample) & is.null(group) & is.null(clusters)) 
    stop("Please provide either a SingleCellExperiment object or Seurat\n        object with required annotation metadata, or explicitly provide\n        clusters, sample and group information")
  if ((is.null(clusters) | is.null(sample) | is.null(group)) & 
      !is.null(x)) {
    if (is(x, "SingleCellExperiment")) 
      y <- .extractSCE(x)
    if (is(x, "Seurat")) 
      y <- .extractSeurat(x)
    clusters <- y$clusters
    sample <- y$sample
    group <- y$group
  }
  if (is.null(transform)) 
    transform <- "logit"
  prop.list <- getTransformedProps(clusters, sample, transform)
  baseline.props <- table(clusters)/sum(table(clusters))
  group.coll <- table(sample, group)
  design <- matrix(as.integer(group.coll != 0), ncol = ncol(group.coll))
  colnames(design) <- colnames(group.coll)
  if (ncol(design) == 2) {
    message("group variable has 2 levels, t-tests will be performed")
    contrasts <- c(1, -1)
    out <- propeller.ttest(prop.list, design, contrasts = contrasts, 
                           robust = robust, trend = trend, sort = FALSE)
    out <- data.frame(BaselineProp = baseline.props, out)
    return(out[order(out$P.Value), ])
  }
  else if (ncol(design) >= 2) {
    message("group variable has > 2 levels, ANOVA will be performed")
    coef <- seq_len(ncol(design))
    out <- propeller.anova(prop.list, design, coef = coef, 
                           robust = robust, trend = trend, sort = FALSE)
    out <- data.frame(BaselineProp = as.vector(baseline.props), 
                      out)
    return(out[order(out$P.Value), ])
  }
}


getTransformedProps = function (clusters = clusters, sample = sample, transform = NULL) 
{
  if (is.null(transform)) 
    transform <- "logit"
  tab <- table(sample, clusters)
  props <- tab/rowSums(tab)
  if (transform == "asin") {
    message("Performing arcsin square root transformation of proportions")
    prop.trans <- asin(sqrt(props))
  }
  else if (transform == "logit") {
    message("Performing logit transformation of proportions")
    props.pseudo <- (tab + 0.5)/rowSums(tab + 0.5)
    prop.trans <- log(props.pseudo/(1 - props.pseudo))
  }
  return(list(Counts = t(tab), TransformedProps = t(prop.trans), 
              Proportions = t(props)))
}
propeller.ttest = function (prop.list = prop.list, design = design, contrasts = contrasts, 
                            robust = robust, trend = trend, sort = sort) 
{
  prop.trans <- prop.list$TransformedProps
  prop <- prop.list$Proportions
  if (nrow(prop.trans) <= 2) {
    message("Setting robust to FALSE for eBayes for less than 3 cell types")
    robust <- FALSE
  }
  fit <- lmFit(prop.trans, design)
  fit.cont <- contrasts.fit(fit, contrasts = contrasts)
  fit.cont <- eBayes(fit.cont, robust = robust, trend = trend)
  if (length(contrasts) == 2) {
    fit.prop <- lmFit(prop, design)
    z <- apply(fit.prop$coefficients, 1, function(x) x^contrasts)
    RR <- apply(z, 2, prod)
  }
  else {
    new.des <- design[, contrasts != 0]
    fit.prop <- lmFit(prop, new.des)
    new.cont <- contrasts[contrasts != 0]
    z <- apply(fit.prop$coefficients, 1, function(x) x^new.cont)
    RR <- apply(z, 2, prod)
  }
  fdr <- p.adjust(fit.cont$p.value[, 1], method = "BH")
  out <- data.frame(PropMean = fit.prop$coefficients, PropRatio = RR, 
                    Tstatistic = fit.cont$t[, 1], P.Value = fit.cont$p.value[, 
                                                                             1], FDR = fdr)
  if (sort) {
    o <- order(out$P.Value)
    return(out[o, ])
  }
  else {
    return(out)
  }
}
propeller.anova = function (prop.list = prop.list, design = design, coef = coef, 
                            robust = robust, trend = trend, sort = sort) 
{
  prop.trans <- prop.list$TransformedProps
  prop <- prop.list$Proportions
  fit.prop <- lmFit(prop, design[, coef])
  design[, 1] <- 1
  colnames(design)[1] <- "Int"
  fit <- lmFit(prop.trans, design)
  fit <- eBayes(fit[, coef[-1]], robust = robust, trend = trend)
  p.value <- fit$F.p.value
  fdr <- p.adjust(fit$F.p.value, method = "BH")
  out <- data.frame(PropMean = fit.prop$coefficients, Fstatistic = fit$F, 
                    P.Value = p.value, FDR = fdr)
  if (sort) {
    o <- order(out$P.Value)
    out[o, ]
  }
  else out
}
