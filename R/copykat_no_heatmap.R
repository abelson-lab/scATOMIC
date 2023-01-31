#' copy_kat_no_heatmap
#'
#' @param rawmat matrix of counts
#' @param id.type set S if gene names
#' @param cell.line set to no for tissues
#' @param ngene.chr minimum nchrm, for my tool i use 0 as i want every cell included
#' @param LOW.DR minimal population fractions of genes for smoothing.
#' @param UP.DR minimal population fractions of genes for segmentation.
#' @param win.size minimal window sizes for segmentation.
#' @param norm.cell.names vector of normal cells use the blood and stromal from the classifier
#' @param KS.cut segmentation parameters, input 0 to 1; larger looser criteria.
#' @param sam.name sample name.
#' @param distance default eucledian
#' @param n.cores for parallel computing
#'
#' @return returns a CNV matrix for use in create_summary_matrix
#' @export
#function is code modified from Gao R. 2021 https://github.com/navinlabcode/copykat
copy_kat_no_heatmap <- function (rawmat = rawdata,summary_matrix,  id.type = "S", cell.line = "no",
                                 ngene.chr = 5, LOW.DR = 0.05, UP.DR = 0.1, win.size = 25, KS.cut = 0.1, sam.name = "", distance = "euclidean",
                                 n.cores = (detectCores() - 1))
{
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
    n.cores = 1
  }
  start_time <- Sys.time()
  if(length(which(duplicated(row.names(rawmat))))> 0){
    rawmat <- rawmat[-which(duplicated(row.names(rawmat))),]
  }
  norm.cell.names = row.names(summary_matrix)[which(summary_matrix$layer_1 == "Blood_Cell" |
                                                      summary_matrix$layer_2 %in% c("Endothelial Cells", "Stromal Cell", "Oligodendrocytes"))]
  set.seed(1)
  sample.name <- paste(sam.name, "_copykat_", sep = "")
  print("step1: read and filter data ...")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data",
              sep = ""))
  genes.raw <- apply(rawmat, 2, function(x) (sum(x > 0)))
  if (sum(genes.raw > 200) == 0)
    stop("none cells have more than 200 genes")
  if (sum(genes.raw < 100) > 1) {
    rawmat <- rawmat[, -which(genes.raw < 200)]
    print(paste("filtered out ", sum(genes.raw <= 200), " cells with less than 200 genes; remaining ",
                ncol(rawmat), " cells", sep = ""))
  }
  der <- apply(rawmat, 1, function(x) (sum(x > 0)))/ncol(rawmat)
  if (sum(der > LOW.DR) >= 1) {
    rawmat <- rawmat[which(der > LOW.DR), ]
    print(paste(nrow(rawmat), " genes past LOW.DR filtering",
                sep = ""))
  }
  WNS1 <- "data quality is ok"
  if (nrow(rawmat) < 7000) {
    WNS1 <- "low data quality"
    UP.DR <- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }
  print("step 2: annotations gene coordinates ...")
  anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type)
  anno.mat <- anno.mat[order(anno.mat$abspos, decreasing = FALSE),
  ]
  HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
  toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]),
                                             HLAs))
  if (length(toRev) > 0) {
    anno.mat <- anno.mat[-toRev, ]
  }
  ToRemov2 <- NULL
  for (i in 8:ncol(anno.mat)) {
    cell <- cbind(anno.mat$chromosome_name, anno.mat[, i])
    cell <- cell[cell[, 2] != 0, ]
    if (length(as.numeric(cell)) < 5) {
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    else if (length(rle(cell[, 1])$length) < 23 | min(rle(cell[,
                                                               1])$length) < ngene.chr) {
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    i <- i + 1
  }
  if (length(ToRemov2) == (ncol(anno.mat) - 7))
    stop("all cells are filtered")
  if (length(ToRemov2) > 0) {
    anno.mat <- anno.mat[, -which(colnames(anno.mat) %in%
                                    ToRemov2)]
  }
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat <- log(sqrt(rawmat3) + sqrt(rawmat3 + 1))
  norm.mat <- apply(norm.mat, 2, function(x) (x <- x - mean(x)))
  colnames(norm.mat) <- colnames(rawmat3)
  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c) {
    model <- dlm::dlmModPoly(order = 1, dV = 0.16, dW = 0.001)
    x <- dlm::dlmSmooth(norm.mat[, c], model)$s
    x <- x[2:length(x)]
    x <- x - mean(x)
  }
  test.mc <- parallel::mclapply(1:ncol(norm.mat), dlm.sm, mc.cores = n.cores)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat),
                            byrow = FALSE)
  colnames(norm.mat.smooth) <- colnames(norm.mat)
  print("step 4: measuring baselines ...")
  if (cell.line == "yes") {
    print("running pure cell line mode")
    relt <- baseline.synthetic(norm.mat = norm.mat.smooth,
                               min.cells = 10, n.cores = n.cores)
    norm.mat.relat <- relt$expr.relat
    CL <- relt$cl
    WNS <- "run with cell line mode"
    preN <- NULL
  }
  else if (length(norm.cell.names) > 1) {
    NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in%
                                                    norm.cell.names)])
    print(paste(NNN, " known normal cells found in dataset",
                sep = ""))
    if (NNN == 0)
      stop("known normal cells provided; however none existing in testing dataset")
    print("run with known normal...")
    basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in%
                                             norm.cell.names)], 1, median)
    print("baseline is from known input")
    d <- parallelDist::parDist(t(norm.mat.smooth), threads = n.cores,
                               method = "euclidean")
    km <- 6
    fit <- hclust(d, method = "ward.D2")
    CL <- cutree(fit, km)
    while (!all(table(CL) > 5)) {
      km <- km - 1
      CL <- cutree(fit, k = km)
      if (km == 2) {
        break
      }
    }
    WNS <- "run with known normal"
    preN <- norm.cell.names
    norm.mat.relat <- norm.mat.smooth - basel
  }
  else {
    basa <- baseline.norm.cl(norm.mat.smooth = norm.mat.smooth,
                             min.cells = 5, n.cores = n.cores)
    basel <- basa$basel
    WNS <- basa$WNS
    preN <- basa$preN
    CL <- basa$cl
    if (WNS == "unclassified.prediction") {
      Tc <- colnames(rawmat)[which(as.numeric(apply(rawmat[which(rownames(rawmat) %in%
                                                                   c("PTPRC", "LYZ", "PECAM1")), ], 2, mean)) >
                                     1)]
      length(Tc)
      preN <- intersect(Tc, colnames(norm.mat.smooth))
      if (length(preN) > 5) {
        print("start manual mode")
        WNS <- paste("copykat failed in locating normal cells; manual adjust performed with ",
                     length(preN), " immune cells", sep = "")
        print(WNS)
        basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in%
                                                 preN)], 1, mean)
      }
      else {
        basa <- baseline.GMM(CNA.mat = norm.mat.smooth,
                             max.normal = 5, mu.cut = 0.05, Nfraq.cut = 0.99,
                             RE.before = basa, n.cores = n.cores)
        basel <- basa$basel
        WNS <- basa$WNS
        preN <- basa$preN
      }
    }
    norm.mat.relat <- norm.mat.smooth - basel
  }
  DR2 <- apply(rawmat3, 1, function(x) (sum(x > 0)))/ncol(rawmat3)
  norm.mat.relat <- norm.mat.relat[which(DR2 >= UP.DR), ]
  anno.mat2 <- anno.mat[which(DR2 >= UP.DR), ]
  ToRemov3 <- NULL
  for (i in 8:ncol(anno.mat2)) {
    cell <- cbind(anno.mat2$chromosome_name, anno.mat2[,
                                                       i])
    cell <- cell[cell[, 2] != 0, ]
    if (length(as.numeric(cell)) < 5) {
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    else if (length(rle(cell[, 1])$length) < 23 | min(rle(cell[,
                                                               1])$length) < ngene.chr) {
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    i <- i + 1
  }
  if (length(ToRemov3) == ncol(norm.mat.relat))
    stop("all cells are filtered")
  if (length(ToRemov3) > 0) {
    norm.mat.relat <- norm.mat.relat[, -which(colnames(norm.mat.relat) %in%
                                                ToRemov3)]
  }
  CL <- CL[which(names(CL) %in% colnames(norm.mat.relat))]
  CL <- CL[order(match(names(CL), colnames(norm.mat.relat)))]
  print("step 5: segmentation...")
  results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat, bins = win.size,
                      cut.cor = KS.cut, n.cores = n.cores)
  if (length(results$breaks) < 25) {
    print("too few breakpoints detected; decreased KS.cut to 50%")
    results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat,
                        bins = win.size, cut.cor = 0.5 * KS.cut, n.cores = n.cores)
  }
  if (length(results$breaks) < 25) {
    print("too few breakpoints detected; decreased KS.cut to 75%")
    results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat,
                        bins = win.size, cut.cor = 0.5 * 0.5 * KS.cut, n.cores = n.cores)
  }
  if (length(results$breaks) < 25)
    stop("too few segments; try to decrease KS.cut; or improve data")
  colnames(results$logCNA) <- colnames(norm.mat.relat)
  results.com <- apply(results$logCNA, 2, function(x) (x <- x -
                                                         mean(x)))
  RNA.copycat <- cbind(anno.mat2[, 1:7], results.com)

  print("step 6: convert to genomic bins...")
  Aj <- convert.all.bins.hg20(DNA.mat = DNA.hg20, RNA.mat = RNA.copycat,
                              n.cores = n.cores)
  uber.mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
  print("step 7: adjust baseline ...")
  if (cell.line == "yes") {
    mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])

    if (distance == "euclidean") {
      hcc <- hclust(parallelDist::parDist(t(mat.adj), threads = n.cores,
                                          method = distance), method = "ward.D")
    }
    else {
      hcc <- hclust(as.dist(1 - cor(mat.adj, method = distance)),
                    method = "ward.D")
    }

    end_time <- Sys.time()
    print(end_time - start_time)
    reslts <- list(cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
    names(reslts) <- c("CNAmat", "hclustering")
    return(reslts)
  }
  else {
    if (distance == "euclidean") {
      hcc <- hclust(parallelDist::parDist(t(uber.mat.adj),
                                          threads = n.cores, method = distance), method = "ward.D")
    }
    else {
      hcc <- hclust(as.dist(1 - cor(uber.mat.adj, method = distance)),
                    method = "ward.D")
    }
    hc.umap <- cutree(hcc, 2)
    names(hc.umap) <- colnames(results.com)
    cl.ID <- NULL
    for (i in 1:max(hc.umap)) {
      cli <- names(hc.umap)[which(hc.umap == i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i <- i + 1
    }
    com.pred <- names(hc.umap)
    com.pred[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
    com.pred[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "nondiploid"
    names(com.pred) <- names(hc.umap)
    results.com.rat <- uber.mat.adj - apply(uber.mat.adj[,
                                                         which(com.pred == "diploid")], 1, mean)
    results.com.rat <- apply(results.com.rat, 2, function(x) (x <- x -
                                                                mean(x)))
    results.com.rat.norm <- results.com.rat[, which(com.pred ==
                                                      "diploid")]
    dim(results.com.rat.norm)
    cf.h <- apply(results.com.rat.norm, 1, sd)
    base <- apply(results.com.rat.norm, 1, mean)
    adjN <- function(j) {
      a <- results.com.rat[, j]
      a[abs(a - base) <= 0.25 * cf.h] <- mean(a)
      a
    }
    mc.adjN <- parallel::mclapply(1:ncol(results.com.rat),
                                  adjN, mc.cores = n.cores)
    adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat),
                          byrow = FALSE)
    colnames(adj.results) <- colnames(results.com.rat)
    rang <- 0.5 * (max(adj.results) - min(adj.results))
    mat.adj <- adj.results/rang
    print("step 8: final prediction ...")
    if (distance == "euclidean") {
      hcc <- hclust(parallelDist::parDist(t(mat.adj), threads = n.cores,
                                          method = distance), method = "ward.D")
    }
    else {
      hcc <- hclust(as.dist(1 - cor(mat.adj, method = distance)),
                    method = "ward.D")
    }
    hc.umap <- cutree(hcc, 2)
    names(hc.umap) <- colnames(results.com)

    cl.ID <- NULL
    for (i in 1:max(hc.umap)) {
      cli <- names(hc.umap)[which(hc.umap == i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i <- i + 1
    }
    com.preN <- names(hc.umap)
    com.preN[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
    com.preN[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "aneuploid"
    names(com.preN) <- names(hc.umap)
    if (WNS == "unclassified.prediction") {
      com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
      com.preN[which(com.preN == "nondiploid")] <- "c2:aneuploid:low.conf"
    }
    print("step 9: saving results...")
    res <- cbind(names(com.preN), com.preN)
    colnames(res) <- c("cell.names", "copykat.pred")


    end_time <- Sys.time()
    print(end_time - start_time)

    return(res)
  }
}
