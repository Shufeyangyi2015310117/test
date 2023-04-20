#' evaluate.
#'
#' @description
#' calculate the mF1 between true label and predicated label
#'
#' @details calculate the mF1 between true label and predicated label
#' @param true_labe true label.
#' @param pred_lab predicated label
#' @param Indices a index set
#' @return a list.
#' @export
evaluate <- function(true_lab, pred_lab, Indices = NULL){

  #true_lab <- unlist(read.csv(TrueLabelsPath))
  #pred_lab <- unlist(read.csv(PredLabelsPath))

  if (! is.null(Indices)){
    true_lab <- true_lab[Indices]
    pred_lab <- pred_lab[Indices]
  }

  unique_true <- unlist(unique(true_lab))
  unique_pred <- unlist(unique(pred_lab))

  unique_all <- unique(c(unique_true,unique_pred))
  conf <- table(true_lab,pred_lab)
  pop_size <- rowSums(conf)

  pred_lab = gsub('Node..','Node',pred_lab)

  conf_F1 <- table(true_lab,pred_lab,exclude = c('unassigned','Unassigned', 'rand','Node','ambiguous','Unknown'))

  F1 <- vector()
  sum_acc <- 0

  for (i in c(1:length(unique_true))){
    findLabel = colnames(conf_F1) == row.names(conf_F1)[i]
    if(sum(findLabel)){
      prec <- conf_F1[i,findLabel] / colSums(conf_F1)[findLabel]
      rec <- conf_F1[i,findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0){
        F1[i] = 0
      } else{
        F1[i] <- (2*prec*rec) / (prec + rec)
      }
      sum_acc <- sum_acc + conf_F1[i,findLabel]
    } else {
      F1[i] = 0
    }
  }

  pop_size <- pop_size[pop_size > 0]

  names(F1) <- names(pop_size)

  med_F1 <- median(F1)

  total <- length(pred_lab)
  num_unlab <- sum(pred_lab == 'unassigned') + sum(pred_lab == 'Unassigned') + sum(pred_lab == 'rand') + sum(pred_lab == 'Unknown') + sum(pred_lab == 'unknown') + sum(pred_lab == 'Node') + sum(pred_lab == 'ambiguous')
  per_unlab <- num_unlab / total

  acc <- sum_acc/sum(conf_F1)

  result <- list(Conf = conf, MedF1 = med_F1, F1 = F1, Acc = acc, PercUnl = per_unlab, PopSize = pop_size)

  return(result)
}

