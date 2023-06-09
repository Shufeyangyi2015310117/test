#' wpca.
#'
#' @description
#' wpca weighted PCA
#'
#' @details weighted PCA
#' @param X a normalized gene expression matrix.
#' @param q the number of dimension for low-dimensional embedding
#' @param weighted a logical value specifying whether or not to be weighted
#' @return a list.
#' @export
wpca <- function(X, q, weighted=TRUE){
  if(!is.matrix(X)) stop("wpca: X must be a matrix!")
  if(q< 1) stop("wpca: q must be a positive integer!")
  X <- scale(X, scale=F) # centralize
  out <- wpcaCpp(X, q, weighted)
  return(out)
}

#' marker_list_to_mat.
#'
#' @description
#' marker_list_to_mat, a function that transform list to cell type gene matrix.
#'
#' @details marker_list_to_mat, a function that transform list to cell type gene matrix.
#' @param marker_list a marker list.
#' @param include_other is logical value specifying whether or not include unknown
#' @return a list contains.
#' @export
marker_list_to_mat <- function(marker_list, include_other = TRUE) {
  cell_types <- names(marker_list)

  if(is.null(cell_types)) {
    warning("Marker list has no cell type names - replacing with generics")
    cell_types <- paste0("cell_type_", seq_along(marker_list))
    names(marker_list) <- cell_types
  }

  genes <- sort(unique(unlist(marker_list)))
  genes <- genes[nchar(genes) > 0]

  n_cell_types <- length(cell_types)
  n_genes <- length(genes)

  mat <- matrix(0, nrow = n_cell_types, ncol = n_genes)
  colnames(mat) <- genes
  rownames(mat) <- cell_types

  for(cell_type in names(marker_list)) {
    mat[cell_type,] <- genes %in% marker_list[[cell_type]]
  }

  if(include_other) {
    mat <- rbind(mat, `other` = 0)
  }

  mat <- t(mat) # Make it gene  by cell type

  mat
}



#' SpatialAnno.
#'
#' @description
#' SpatialAnno, an efficient and accurate annotation method for spatial transcriptomics datasets, with capability of effectively leveraging a large number of non-marker genes with “qualitative” information about marker genes, without using a reference dataset.
#'
#' @details SpatialAnno, an efficient and accurate annotation method for spatial transcriptomics datasets, with capability of effectively leveraging a large number of non-marker genes with “qualitative” information about marker genes, without using a reference dataset.
#' @param X a normalized matrix. The row is for spots/cells, and the column is for genes.
#' @param Adj_sp  a integer vector used in SC.MEB. The default is 2:10
#' @param markers  a list that contain cell types and corresponding cell-type specific markers.
#' @param initial  character specifying which method to generate initial value. In general, We recommend to use the method of SCINA for high-dimensional non-markers and the method of scSorter for low-dimensional non-marker. The default is SCINA.
#' @param Unknown  logical value specifying whether or not include unknown
#' @param xi_int  a initial value of smooth parameter
#' @param xi_grid a pre-defined vector of smooth parameter
#' @param mode a string specifying which mode to use. It has four modes: "full", "annotation", "spatial+annotation" and "annotation+factor". The "full" mode means SpatialAnno runs with all three components, marker component, spatial component, and factor component. The "annotation" mode means SpatialAnno runs with only marker component. The "spatial+annotation" means SpatialAnno runs with marker and spatial components. The "annotation+factor" mode means SpatialAnno runs with marker and non-marker components. The default is "full".
#' @param q the dimension of embedding generated by SpatialAnno. The default is 15
#' @return a list.
#' @export
SpatialAnno <- function(X, Adj_sp, markers, initial = c("SCINA", "scSorter"), Unknown = TRUE, xi_int = 1.5, xi_grid=seq(0.1, 2.5, by=0.2), mode = "full", q = 15){


  rho <- marker_list_to_mat(markers, Unknown)
  df_all = as.data.frame(matrix(0,0,2))
  colnames(df_all) = c("Marker", "Type")
  for (i in 1:(dim(rho)[2]-1)){
    print(rownames(rho)[rho[,i]==1])
    df = as.data.frame(rownames(rho)[rho[,i]==1])
    colnames(df) = "Marker"
    df$Type = colnames(rho)[i]
    df_all = rbind(df_all, df)
  }
  anno = df_all[,c(2,1)]

  # --------------------------------------------------------------------
  # Unknown is a swift controlling whether novel celltype is allowed.
  Unknown = TRUE
  anno_processed = scSorter:::design_matrix_builder(anno, weight=2)
  dat <- scSorter:::data_preprocess(t(X), anno_processed)
  if (Unknown == TRUE) dat$designmat$Unknown = 0
  rho <- as.matrix(dat$designmat)
  m <- nrow(rho)
  X_m <- t(dat$dat[1:m, ])
  X_u <- t(dat$dat[-(1:m), ])
  n <- nrow(X_m)
  K <- ncol(rho)


  if (initial == "SCINA"){

    results = adjSCINA(t(X), markers, max_iter = 100, convergence_n = 10, rm_overlap = 0, allow_unknown=TRUE)

    bet_min <- numeric(length(results$theta))
    bet_max <- numeric(length(results$theta))
    alpha_int <- numeric(m)
    bet_int <- matrix(0, nrow = m, ncol = K)
    for(i in 1:length(results$theta)){
      bet_max[i] <- max(results$theta[[i]]$mean[,1] - results$theta[[i]]$mean[,2])
      bet_min[i] <- min(results$theta[[i]]$mean[,1] - results$theta[[i]]$mean[,2])

      bet_int[rho[,i]!=0, i] <- results$theta[[i]]$mean[, 1] - results$theta[[i]]$mean[, 2]
      alpha_int[rho[,i]!=0] <- results$theta[[i]]$mean[, 2]
    }
    lfc <- median(bet_min)


    ################################################
    ## SCINA initial value
    colnames(rho)[K] = "Unknown"
    y_hat <- results$cell_labels;  y_hat[y_hat == "unknown"] = "Unknown"

    if (length(unique(y_hat))<K){
      idx = match(names(markers), unique(y_hat))
      idx2 = which(is.na(idx)==TRUE)
      if (length(idx2) > 0){
        y_hat[1:length(names(markers)[idx2])]=names(markers)[idx2]
      }
      if (is.na(match("Unknown", unique(y_hat)))){
        y_hat[length(names(markers)[idx2])+1]="Unknown"
      }
    }

    y_int <- match(y_hat, colnames(rho))
    R_int <- matrix(0, nrow = n, ncol = K)
    for(i in 1:nrow(X_m)) R_int[i, match(y_hat, colnames(rho))[i]] <- 1
    #all(apply(R_int, 1, which.max) == match(rts$Pred_Type,colnames(rho)))
    Pi_u_int <- colMeans(R_int)
    #alpha_int <- fit_jig$alpha; bet_int <- fit_jig$bet; sigma_int <- fit_jig$sigma
    #idx = match(rownames(rho), toupper(rownames(fit$mle_params$beta)))
    #bet_int <- fit$mle_params$delta[idx, ] * rho
    #alpha_int <- fit$mle_params$beta[idx,1]
    sigma_int <- update_sigma(X_m, rho, R_int, alpha_int, bet_int)


    princ <- wpca(X_u, q=q, weighted=TRUE)
    Lam_u_int <- princ$Lam_vec
    W_u_int <- princ$loadings
    hZ <- princ$PCs
    n_c <- colSums(R_int)
    Mu_u_int <- t(sapply(1:K, function(k) 1/n_c[k] * colSums(R_int[, k] * hZ)))
    q <- ncol(Mu_u_int)
    Sgm_u_int <- init_Sgm(R_int, hZ, matrix(0, ncol=q, nrow=q), Mu_u_int, FALSE)

  }else if(initial == "scSorter"){

    library(scSorter)
    dat <- t(X)
    rts <- scSorter(dat, anno)
    lfc = 0

    y_hat <- rts$Pred_Type;

    if (length(unique(y_hat))<K){
      idx = match(names(markers), unique(y_hat))
      idx2 = which(is.na(idx)==TRUE)
      if (length(idx2) > 0){
        y_hat[1:length(names(markers)[idx2])]=names(markers)[idx2]
      }
      if (is.na(match("Unknown", unique(y_hat)))){
        y_hat[length(names(markers)[idx2])+1]="Unknown"
      }
    }

    # -----------------------------------------------------------------
    # initial values with unknown
    y_int <- match(y_hat, colnames(rho))
    R_int <- matrix(0, nrow = n, ncol = K)
    for(i in 1:n) R_int[i, match(y_hat, colnames(rho))[i]] <- 1
    idx_colsum_R_int <- which(colSums(R_int)==0)
    if (length(idx_colsum_R_int) != 0)  {
      for(i in 1:length(idx_colsum_R_int)) {
        R_int[i, idx_colsum_R_int[i]] = 1
      }
    }
    Pi_u_int <- colMeans(R_int)
    bet_int <- rts$Pred_param[1:m, 1:(K-1)] * rho[,(K-1)]
    bet_int <- cbind(bet_int, matrix(0, m, 1))
    alpha_int <- apply((rts$Pred_param[1:m, 1:(K-1)] -  bet_int[,1:(K-1)]), 1, max)
    sigma_int <- update_sigma(X_m, rho, R_int, alpha_int, bet_int)


    princ <- wpca(X_u, q=15, weighted=TRUE)
    Lam_u_int <- princ$Lam_vec
    W_u_int <- princ$loadings
    hZ <- princ$PCs
    n_c <- colSums(R_int)
    Mu_u_int <- t(sapply(1:K, function(k) 1/n_c[k] * colSums(R_int[, k] * hZ)))
    q <- ncol(Mu_u_int)
    Sgm_u_int <- init_Sgm(R_int, hZ, matrix(0, ncol=q, nrow=q), Mu_u_int, FALSE)
    # -----------------------------------------------------------------

  }else{
    print("You specify the wrong way of generating initial value")
    break
  }

  if (mode == "full"){
    fit <- icmem(X_m, X_u, Adj_sp, rho, lfc,
                 y_int, Pi_u_int*0, xi_int, xi_grid,
                 alpha_int, bet_int, sigma_int,
                 Mu_u_int, W_u_int, Sgm_u_int, Lam_u_int,
                 300, 10, 1e-6, TRUE,
                 FALSE, FALSE)
    return(fit)
  }else if (mode == "annotation"){
    xi_int = 0
    xi_grid = 0
    fit <- icmem_anno(X_m, X_u, Adj_sp, rho, lfc,
                 y_int, Pi_u_int*0, xi_int, xi_grid,
                 alpha_int, bet_int, sigma_int,
                 Mu_u_int, W_u_int, Sgm_u_int, Lam_u_int,
                 300, 10, 1e-6, TRUE,
                 FALSE, FALSE)
    return(fit)
  }else if(mode == "spatial+annotation"){
    fit <- icmem_spatialanno(X_m, X_u, Adj_sp, rho, lfc,
                      y_int, Pi_u_int*0, xi_int, xi_grid,
                      alpha_int, bet_int, sigma_int,
                      Mu_u_int, W_u_int, Sgm_u_int, Lam_u_int,
                      300, 10, 1e-6, TRUE,
                      FALSE, FALSE)
    return(fit)
  }else if(mode == "annotation+factor"){
    xi_int = 0
    xi_grid = 0
    fit <- icmem_annofactor(X_m, X_u, Adj_sp, rho, lfc,
                               y_int, Pi_u_int*0, xi_int, xi_grid,
                               alpha_int, bet_int, sigma_int,
                               Mu_u_int, W_u_int, Sgm_u_int, Lam_u_int,
                               300, 10, 1e-6, TRUE,
                               FALSE, FALSE)
    return(fit)
  }else{
    print("Please specify the correct mode of SpatialAnno")
  }



}

