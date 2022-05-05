#' scATOMICTree
#'
#' @param predictions_list list of layered predictions, from run_scATOMIC()
#' @param summary_matrix summary matrix of results, from create_summary_matrix()
#' @param interactive_mode whether to return collapsibleTree chart, if F flow chart returned
#' @param collapsed if interactive mode is T, whether to collapse nodes or open all nodes
#' @param save_results logical if you want to save the tree as an html file
#' @param save_dir directory to save html file, we recommend working directory
#' @param project_name prefix name of project for save file, saved as project_name_annotation_results_tree.html
#' @param height height in pixels (optional, defaults to automatic sizing)
#' @param width width in pixels (optional, defaults to automatic sizing)
#' @param fontSize font size of the label text in pixels
#' @param linkLength length of the horizontal links that connect nodes in pixels
#' @param ncell_size logical whether to scale size of dots in interactive mode by number of cells (set to false if using relatively homogeneous dataset)
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' #load data
#' lung_cancer_demo_data <- demo_lung_data
#' # Quality control filtering
#' pct_mt <- colSums(lung_cancer_demo_data[grep("^MT-", row.names(lung_cancer_demo_data)),])/colSums(lung_cancer_demo_data) * 100
#' nFeatureRNA <- apply(lung_cancer_demo_data, 2,function(x){length(which(x != 0))})
#' lung_cancer_demo_data <- lung_cancer_demo_data[, names(which(pct_mt < 25))]
#' lung_cancer_demo_data <- lung_cancer_demo_data[, intersect(names(which(nFeatureRNA > 500)), colnames(lung_cancer_demo_data))]
#' cell_predictions <- run_scATOMIC(lung_cancer_demo_data)
#' #create results matrix - no CNVs
#' results_matrix <- create_summary_matrix(prediction_list = cell_predictions, use_CNVs = F, modify_results = T, mc.cores = 1, raw_counts = lung_cancer_demo_data, min_prop = 0.5 )
#' #visualize on tree
#' #interactive mode
#' tree_results_lung <- scATOMICTree(predictions_list = cell_predictions, summary_matrix = results_lung,
#' interactive_mode = T, collapsed = T, save_results = F,height = 700, width = 1000)
#' #non interactive mode
#' tree_results_non_interactive <- scATOMICTree(predictions_list = cell_predictions, summary_matrix = results_lung,
#' interactive_mode = F, save_results = F)
#' }
scATOMICTree <- function(predictions_list, summary_matrix, interactive_mode = T, collapsed = T, save_results = T, save_dir = getwd(), project_name, width =NULL, height= NULL,fontSize = 10,linkLength = 150, ncell_size = T){
  layer_median_scores <- list()
  layer_IQRs <- list()
  for(i in 1:length(predictions_list)){
    layer_pred <- predictions_list[[i]]
    layer <- names(predictions_list)[i]
    predicted_cells <- levels(as.factor(layer_pred$predicted_tissue_with_cutoff))
    per_cell_type_median_score <- c()
    per_cell_type_median_IQR <- c()
    #get IQR for each number and take max
    for(m in 1:length(predicted_cells)){
      cell_use <- predicted_cells[m]
      numeric_column_index <- colnames(select_if(layer_pred, is.numeric))
      layer_pred_cell_use <- layer_pred[which(layer_pred$predicted_tissue_with_cutoff == cell_use), numeric_column_index]
      per_cell_type_median_score[m] <- max(unlist(lapply(layer_pred_cell_use, median), use.names = T))[1]
      names(per_cell_type_median_score)[m] <- cell_use
      per_cell_type_median_IQR[m] <- unlist(lapply(layer_pred_cell_use, IQR), use.names = T)[which(unlist(lapply(layer_pred_cell_use, median), use.names = T) == per_cell_type_median_score[m])][1]
      names(per_cell_type_median_IQR)[m] <- cell_use
    }
    layer_median_scores[[i]] <- per_cell_type_median_score
    layer_IQRs[[i]] <- per_cell_type_median_IQR
    names(layer_median_scores)[i] <- layer
    names(layer_IQRs)[i] <- layer
  }
  #check if CNV mode was run
   if(length(grep("pan_cancer_cluster", colnames(summary_matrix))) > 0) {
    normal_index <- which(summary_matrix$layer_6 != summary_matrix$scATOMIC_pred &
                            summary_matrix$scATOMIC_pred == "Normal Tissue Cell")

    summary_matrix[normal_index, c("layer_1", "layer_2","layer_3","layer_4", "layer_5","layer_6", "scATOMIC_pred")] <- "Normal Tissue Cell"
    unknown_cancer_index <- which(summary_matrix$layer_6 != summary_matrix$scATOMIC_pred &
                                    summary_matrix$scATOMIC_pred == "Unknown Cancer Cell")

    summary_matrix[unknown_cancer_index, c("layer_1", "layer_2","layer_3","layer_4", "layer_5","layer_6", "scATOMIC_pred")] <- "Unknown Cancer Cell"

  }else{
    unknown_cancer_index <- which(summary_matrix$layer_6 != summary_matrix$scATOMIC_pred &
                                    summary_matrix$scATOMIC_pred == "Unknown Cancer Cell")

    summary_matrix[unknown_cancer_index, c("layer_1", "layer_2","layer_3","layer_4", "layer_5","layer_6", "scATOMIC_pred")] <- "Unknown Cancer Cell"

  }
  unknown_index <- which(summary_matrix$scATOMIC_pred == "Unknown Cells")
  summary_matrix[unknown_index, c("layer_1", "layer_2","layer_3","layer_4", "layer_5")] <- "unclassified_any_cell"
  summary_matrix[unknown_index, c("layer_6", "scATOMIC_pred")] <- "Any Cell"


  cancer_converted_index <- which(summary_matrix$layer_6 != summary_matrix$scATOMIC_pred &
                                    summary_matrix$scATOMIC_pred != "Normal Tissue Cell")
  cancer_type <- levels(as.factor(summary_matrix$scATOMIC_pred[grep("Cancer Cell|oma$", summary_matrix$scATOMIC_pred)]))
  if(length(cancer_type) ==1){
    summary_matrix[cancer_converted_index,c("layer_1", "layer_2","layer_3","layer_4", "layer_5", "layer_6") ]<-
      summary_matrix[grep("Cancer Cell|oma$", summary_matrix$scATOMIC_pred)[1],c("layer_1", "layer_2","layer_3","layer_4", "layer_5","layer_6")  ]
  }
  collapsed_results <- summary_matrix[,c( "layer_1","layer_2","layer_3","layer_4","layer_5","layer_6",
                                          "scATOMIC_pred")]

  collapsed_results <- ddply(collapsed_results,.(layer_1,layer_2,layer_3,layer_4,layer_5,layer_6,scATOMIC_pred),nrow)
  colnames(collapsed_results)[8] <- "Number of Cells"

  median_scores <- unlist(layer_median_scores,use.names = T)
  index_name_change <- which(!names(median_scores) %in% c("layer_4_CD4_CD8.CD8+ T cell", "layer_4_CD8_NK.CD8+ T cell"))

  names(median_scores)[index_name_change] <- paste0(substr(names(median_scores)[index_name_change], 1,7),
                                                    "__", unlist(lapply(strsplit(names(median_scores)[index_name_change], "\\."), function(x){x[2]}), use.names = F))
  names(median_scores)[-index_name_change] <- gsub("\\.", "__", names(median_scores)[-index_name_change] )
  names(median_scores) <- gsub("Cancer$", "Cancer Cell",names(median_scores))
  names(median_scores) <- gsub("Tissue_Cell_Normal_or_Cancer Cell", "Tissue_Cell_Normal_or_Cancer",names(median_scores))


  IQR_scores <- unlist(layer_IQRs,use.names = F)
  names(IQR_scores) <- names(median_scores)

  collapsed_results[which(collapsed_results$layer_3 == "Soft Tissue or Neuro Cancer Cell"),"layer_3"] <- "Skin_Lung_Soft Tissue or Neuro Cancer Cell"
  collapsed_results_na <- collapsed_results

  collapsed_results_na$layer_0 <- "Any cell"
  #add score IQRS

  collapsed_results_na <- collapsed_results_na[,c( "layer_0","layer_1","layer_2","layer_3","layer_4","layer_5","layer_6",
                                                   "scATOMIC_pred","Number of Cells" )]
  for(i in 1:nrow(collapsed_results_na)){
    index_dup <- which(duplicated(as.character(collapsed_results_na[i,])))
    #index_dup[length(index_dup)] <- index_dup[1]-1

    collapsed_results_na[i,index_dup] <- NA

  }


  for(i in 2:7){
    cell_types <- levels(as.factor(collapsed_results_na[,i]))
    layer_use <- colnames(collapsed_results_na)[i]
    if(length(cell_types)>0){
      for(p in 1:length(cell_types)){
        if(cell_types[p] %in%c("Normal Tissue Cell", "Unknown Cancer Cell")){

        } else if(cell_types[p] != "CD8+ T cell"){
          if(cell_types[p] == "Skin_Lung_Soft Tissue or Neuro Cancer Cell"){
            median_score_cell <- median_scores[paste0(layer_use, "__Soft Tissue or Neuro Cancer Cell")]
            IQR_score_cell <- IQR_scores[paste0(layer_use, "__Soft Tissue or Neuro Cancer Cell")]

            collapsed_results_na[which(collapsed_results_na[,i] == cell_types[p]), i] <- paste0(cell_types[p], "\n", "median score ", round(median_score_cell, digits = 3), " IQR ", round(IQR_score_cell, digits = 3))
          }else{
            median_score_cell <- median_scores[paste0(layer_use, "__", cell_types[p])]
            IQR_score_cell <- IQR_scores[paste0(layer_use, "__", cell_types[p])]

            collapsed_results_na[which(collapsed_results_na[,i] == cell_types[p]), i] <- paste0(cell_types[p], "\n", "median score ", round(median_score_cell, digits = 3), " IQR ", round(IQR_score_cell, digits = 3))
          }
        } else{
          index_CD8 <- which(collapsed_results_na[,i] == cell_types[p])
          for(d in index_CD8){
            if(grepl("^CD4 or CD8 T cell",collapsed_results_na[d,i-1])){
              median_score_cell <- median_scores["layer_4_CD4_CD8__CD8+ T cell"]
              IQR_score_cell <- IQR_scores["layer_4_CD4_CD8__CD8+ T cell"]
              collapsed_results_na[d, i] <- paste0(cell_types[p], "\n", "median score ", round(median_score_cell, digits = 3), " IQR ", round(IQR_score_cell, digits = 3))
            } else{
              median_score_cell <- median_scores["layer_4_CD8_NK__CD8+ T cell"]
              IQR_score_cell <- IQR_scores["layer_4_CD8_NK__CD8+ T cell"]
              collapsed_results_na[d, i] <- paste0(cell_types[p], "\n", "median score ", round(median_score_cell, digits = 3), " IQR ", round(IQR_score_cell, digits = 3))

            }
          }



        }
      }
    }
    index_na <- grep("score NA", collapsed_results_na[,i])
    if(length(index_na) > 0){
      collapsed_results_na[index_na,i] <- NA
    }
  }


  #add number of cells...
  for(i in 1:8){
    if (i == 1){
      number_of_cells <- sum(collapsed_results_na[,"Number of Cells"])
      collapsed_results_na[, i] <- paste0("Any cell", "\n", number_of_cells, " cells")

    } else{
      cell_types <- levels(as.factor(collapsed_results_na[,i]))
      if(length(cell_types)>0){
        for(p in 1:length(cell_types)){
          if(grepl("^CD8+ T cell\n", cell_types[p])){
            number_of_cells <- sum(collapsed_results_na[which(collapsed_results_na[,i] == cell_types[p]),"Number of Cells"])
            collapsed_results_na[which(collapsed_results_na[,i] == cell_types[p]), i] <- paste0(cell_types[p], "\n", number_of_cells, " cells")
          } else{
            collapsed_results_na[which(collapsed_results_na[,i] == cell_types[p]), i] <- paste0(collapsed_results_na[which(collapsed_results_na[,i] == cell_types[p]), i],
                                                                                                "\n", sum(collapsed_results_na[which(collapsed_results_na[,i] == cell_types[p]), "Number of Cells"]), " cells")
          }
        }
      }
    }
  }



  collapsed_results_na$layer_1 <- gsub("Blood_Cell","Blood cell",collapsed_results_na$layer_1)
  collapsed_results_na$layer_1 <- gsub("Tissue_Cell_Normal_or_Cancer","Non Blood cell",collapsed_results_na$layer_1)
  collapsed_results_na$layer_1 <- gsub("unclassified_any_cell","Any cell",collapsed_results_na$layer_1)
  collapsed_results_na$layer_2 <- gsub("unclassified_any_cell","Any cell",collapsed_results_na$layer_2)
  collapsed_results_na$layer_3 <- gsub("unclassified_any_cell","Any cell",collapsed_results_na$layer_3)
  collapsed_results_na$layer_4 <- gsub("unclassified_any_cell","Any cell",collapsed_results_na$layer_4)
  collapsed_results_na$layer_5 <- gsub("unclassified_any_cell","Any cell",collapsed_results_na$layer_5)
  collapsed_results_na$layer_6 <- gsub("unclassified_any_cell","Any cell",collapsed_results_na$layer_6)

  collapsed_results_na$scATOMIC_pred <- gsub("Any Cell","Any cell",collapsed_results_na$scATOMIC_pred)
  collapsed_results_na$layer_6 <- gsub("Any Cell","Any cell",collapsed_results_na$layer_6)

  collapsed_results_na$layer_2 <- gsub("B_cell_score","B cell or Plasmablast",collapsed_results_na$layer_2)
  collapsed_results_na$layer_2 <- gsub("macrophage_or_dendritic_cell","Macrophage or DC",collapsed_results_na$layer_2)
  collapsed_results_na$layer_2 <- gsub("T_or_NK_lymphocyte","T or NK cell",collapsed_results_na$layer_2)
  collapsed_results_na$layer_2 <- gsub("unclassified_blood_cell","Blood cell",collapsed_results_na$layer_2)
  collapsed_results_na$layer_2 <- gsub("unclassified_normal_or_cancer_tissue","Non Blood cell",collapsed_results_na$layer_2)
  collapsed_results_na$layer_3 <- gsub("unclassified_blood_cell","Blood cell",collapsed_results_na$layer_3)
  collapsed_results_na$layer_3 <- gsub("unclassified_normal_or_cancer_tissue","Non Blood cell",collapsed_results_na$layer_3)
  collapsed_results_na$layer_4 <- gsub("unclassified_blood_cell","Blood cell",collapsed_results_na$layer_4)
  collapsed_results_na$layer_4 <- gsub("unclassified_normal_or_cancer_tissue","Non Blood cell",collapsed_results_na$layer_4)
  collapsed_results_na$layer_5 <- gsub("unclassified_blood_cell","Blood cell",collapsed_results_na$layer_5)
  collapsed_results_na$layer_5 <- gsub("unclassified_normal_or_cancer_tissue","Non Blood cell",collapsed_results_na$layer_5)
  collapsed_results_na$layer_6 <- gsub("unclassified_blood_cell","Blood cell",collapsed_results_na$layer_6)
  collapsed_results_na$layer_6 <- gsub("unclassified_normal_or_cancer_tissue","Non Blood cell",collapsed_results_na$layer_6)
  collapsed_results_na$scATOMIC_pred <- gsub("Blood Cell","Blood cell",collapsed_results_na$scATOMIC_pred)
  collapsed_results_na$scATOMIC_pred <- gsub("^T or NK Cell","T or NK cell",collapsed_results_na$scATOMIC_pred)
  collapsed_results_na$layer_6 <- gsub("Blood Cell","Blood cell",collapsed_results_na$layer_6)
  collapsed_results_na$layer_6 <- gsub("^T or NK Cell","T or NK cell",collapsed_results_na$layer_6)
  collapsed_results_na$layer_3 <- gsub("unclassified_B_cell_or_plasmablast","B cell or Plasmablast",collapsed_results_na$layer_3)
  collapsed_results_na$layer_3 <- gsub("unclassified_T_or_NK_cell","T or NK cell",collapsed_results_na$layer_3)
  collapsed_results_na$layer_4 <- gsub("unclassified_B_cell_or_plasmablast","B cell or Plasmablast",collapsed_results_na$layer_4)
  collapsed_results_na$layer_4 <- gsub("unclassified_T_or_NK_cell","T or NK cell",collapsed_results_na$layer_4)
  collapsed_results_na$layer_5 <- gsub("unclassified_B_cell_or_plasmablast","B cell or Plasmablast",collapsed_results_na$layer_5)
  collapsed_results_na$layer_5 <- gsub("unclassified_T_or_NK_cell","T or NK cell",collapsed_results_na$layer_5)
  collapsed_results_na$layer_6 <- gsub("unclassified_B_cell_or_plasmablast","B cell or Plasmablast",collapsed_results_na$layer_6)
  collapsed_results_na$layer_6 <- gsub("unclassified_T_or_NK_cell","T or NK cell",collapsed_results_na$layer_6)
  collapsed_results_na$layer_4 <- gsub("unclassified_B_cell_or_plasmablast","B cell or Plasmablast",collapsed_results_na$layer_4)
  collapsed_results_na$layer_4 <- gsub("CD8 T or NK cell","NK or CD8 T cell",collapsed_results_na$layer_4)
  collapsed_results_na$layer_5 <- gsub("CD8 T or NK cell","NK or CD8 T cell",collapsed_results_na$layer_5)
  collapsed_results_na$layer_6 <- gsub("CD8 T or NK cell","NK or CD8 T cell",collapsed_results_na$layer_6)
  collapsed_results_na$scATOMIC_pred <- gsub("CD8 T or NK cell","NK or CD8 T cell",collapsed_results_na$scATOMIC_pred)

  collapsed_results_na$scATOMIC_pred <- gsub("Unknown/Normal Tissue Cells","Normal Tissue Cell",collapsed_results_na$scATOMIC_pred)
  collapsed_results_na$layer_6 <- gsub("Unknown/Normal Tissue Cells","Normal Tissue Cell",collapsed_results_na$layer_6)



  collapsed_results_na$pathString <- paste(collapsed_results_na$layer_0, collapsed_results_na$layer_1,
                                           collapsed_results_na$layer_2,
                                           collapsed_results_na$layer_3,
                                           collapsed_results_na$layer_4,
                                           collapsed_results_na$layer_5,
                                           collapsed_results_na$layer_6,
                                           collapsed_results_na$scATOMIC_pred, sep = "layer_sep")

  collapsed_results_na$pathString <- gsub("layer_sepNA", "", collapsed_results_na$pathString)
  collapsed_results_na$pathString <- gsub("NAlayer_sep", "", collapsed_results_na$pathString)






  mydf <- collapsed_results_na[,c("pathString", "Number of Cells")]
  colnames(mydf)[2] <- "Number_of_Cells"




  mydf <- mydf %>%  dplyr::mutate(tree_level = stringr::str_count(string = pathString, pattern = "layer_sep") + 1,
                                  tree_group = stringr::str_replace(string = pathString, pattern = "layer_sep", replacement = ""),
                                  node_type = "decision_node")

  mytree <- data.tree::as.Node(mydf, pathDelimiter = "layer_sep")



  mytree_old <- mytree
  #collapsibleTree( mytree, linkLength = 300, tooltip = T,fillByLevel = T, fontSize = 20,   collapsed = T)
  #collapsibleTree(mytree, root = "Any cell 2997 cells",linkLength = 200)




  for(i in 1:8){
    if(i == 1){
      layer <- colnames(collapsed_results_na)[i]
      cell_types <- levels(as.factor(collapsed_results_na[,i]))
      parent <- NA
      child <-  "All cells"
      Number_of_Cells <- unlist(lapply(strsplit(cell_types, "\n| cells"), function(x){x[2]}), use.names = F)
      Median_Score <- NA
      IQR <- NA
      long_form_mat <- data.frame(parent, child, Number_of_Cells, Median_Score, IQR, layer)
    } else{
      cell_types <- levels(as.factor(collapsed_results_na[,i]))
      if(length(cell_types)>0){
        for(m in 1:length(cell_types)){
          layer <- colnames(collapsed_results_na)[i]
          parent <- unique(collapsed_results_na[which(collapsed_results_na[,i] == cell_types[m]), i-1])
          parent <- unlist(lapply(strsplit(parent, "\n"), function(x){x[1]}), use.names = F)
          child <-  unlist(lapply(strsplit(cell_types[m], "\n"), function(x){x[1]}), use.names = F)
          if(grepl("Unknown Cancer Cell", child)){
            Number_of_Cells <- collapsed_results_na[grep("Unknown Cancer Cell",collapsed_results_na$layer_1),"Number of Cells"]
            Median_Score <- NA
            IQR <- NA
          } else{
            Number_of_Cells <- unlist(lapply(strsplit(cell_types[m], "\n| cells"), function(x){x[length(x)]}), use.names = F)
            Median_Score <- unlist(lapply(strsplit(cell_types[m], "median score | IQR"), function(x){x[2]}), use.names = F)
            IQR <- unlist(lapply(strsplit(cell_types[m], "IQR|\n"), function(x){x[3]}), use.names = F)
          }
          long_form_mat_addition <- data.frame(parent, child, Number_of_Cells, Median_Score, IQR,layer)
          long_form_mat <- rbind(long_form_mat, long_form_mat_addition)
        }
      }
    }
  }
  long_form_mat$parent <- gsub("Any cell", "All cells",long_form_mat$parent )
  long_form_mat[which(is.na(long_form_mat$Median_Score)), "Median_Score"] <- ""
  long_form_mat[which(is.na(long_form_mat$IQR)), "IQR"] <- ""
  layers_long <- levels(as.factor(long_form_mat$layer))
  for(i in 1:length(layers_long)){
    layer_used <- layers_long[i]
    if(layer_used == "layer_0"){

      children <- levels(as.factor(long_form_mat$child[which(long_form_mat$layer == layer_used)]))
      long_form_mat[,"string_1"] <- children

    } else{
      children <- levels(as.factor(long_form_mat$child[which(long_form_mat$layer == layer_used)]))
      long_form_mat[,paste0("string_", i)] <- "undetermined"
      index_children <- which(long_form_mat$child %in% children | long_form_mat$parent %in% children)
      if(length(index_children)>0){
        for(m in index_children){
          if(long_form_mat[m, paste0("string_", (i-1))] == long_form_mat[m, "parent"] ){
            long_form_mat[m,paste0("string_", i)] <- long_form_mat[m,"child"]
          } else if (long_form_mat[m, paste0("string_", (i-1))] != long_form_mat[m, "child"]){
            long_form_mat[m,paste0("string_", i)] <- long_form_mat[m,"parent"]
          }
        }
      }
    }
  }


  for(i in which(long_form_mat$layer == "layer_3")){
    index_true_string_2 <- which(long_form_mat$child == long_form_mat[i,"string_3"] &long_form_mat$parent != long_form_mat[i,"string_3"])
    long_form_mat[i,"string_2"] <- long_form_mat[index_true_string_2, "parent"]
  }
  for(i in which(long_form_mat$layer == "layer_4")){
    index_true_string_3 <- which(long_form_mat$child == long_form_mat[i,"string_4"] & long_form_mat$parent != long_form_mat[i,"string_4"])
    long_form_mat[i,"string_3"] <- long_form_mat[index_true_string_3, "parent"]
    index_true_string_2 <- which(long_form_mat$child == long_form_mat[i,"string_3"] &long_form_mat$parent != long_form_mat[i,"string_3"])
    long_form_mat[i,"string_2"] <- long_form_mat[index_true_string_2, "parent"]
  }
  for(i in which(long_form_mat$layer == "layer_5")){
    index_true_string_4 <- which(long_form_mat$child == long_form_mat[i,"string_5"] &long_form_mat$parent != long_form_mat[i,"string_5"])
    long_form_mat[i,"string_4"] <- long_form_mat[index_true_string_4, "parent"]
    index_true_string_3 <- which(long_form_mat$child == long_form_mat[i,"string_4"] & long_form_mat$parent != long_form_mat[i,"string_4"])
    long_form_mat[i,"string_3"] <- long_form_mat[index_true_string_3, "parent"]
    index_true_string_2 <- which(long_form_mat$child == long_form_mat[i,"string_3"] &long_form_mat$parent != long_form_mat[i,"string_3"])
    long_form_mat[i,"string_2"] <- long_form_mat[index_true_string_2, "parent"]
  }
  for(i in which(long_form_mat$layer == "layer_6")){
    index_true_string_5 <- which(long_form_mat$child == long_form_mat[i,"string_6"] &long_form_mat$parent != long_form_mat[i,"string_6"])
    long_form_mat[i,"string_5"] <- long_form_mat[index_true_string_5, "parent"]
    index_true_string_4 <- which(long_form_mat$child == long_form_mat[i,"string_5"] &long_form_mat$parent != long_form_mat[i,"string_5"])
    long_form_mat[i,"string_4"] <- long_form_mat[index_true_string_4, "parent"]
    index_true_string_3 <- which(long_form_mat$child == long_form_mat[i,"string_4"] & long_form_mat$parent != long_form_mat[i,"string_4"])
    long_form_mat[i,"string_3"] <- long_form_mat[index_true_string_3, "parent"]
    index_true_string_2 <- which(long_form_mat$child == long_form_mat[i,"string_3"] &long_form_mat$parent != long_form_mat[i,"string_3"])
    long_form_mat[i,"string_2"] <- long_form_mat[index_true_string_2, "parent"]
  }
  #fix_endothelial cells

  long_form_mat <- long_form_mat[,c(colnames(long_form_mat)[1:6],sort(colnames(long_form_mat)[grep("^string", colnames(long_form_mat))]))]

  for(i in 1:nrow(long_form_mat)){
    index_un <- which(long_form_mat[i,] != "undetermined")
    if(length(index_un)/(index_un[length(index_un)] - index_un[1] + 1) == 1){
      long_form_mat[i,-index_un] <- NA
    }
    dup <- which(duplicated(na.omit(as.character(long_form_mat[i,grep("string", colnames(long_form_mat))]))) == T)
    if(length(dup) > 0){
      long_form_mat[i,paste0("string_", dup)] <- paste0("Unclassified ", long_form_mat[i,paste0("string_", dup)])
    }
  }



  index_string <- sort(colnames(long_form_mat)[grep("^string", colnames(long_form_mat))])
  for(i in 1:length(index_string)){
    if(i == 1){
      long_form_mat$pathString <- long_form_mat$string_1
    } else{
      long_form_mat$pathString <- paste(long_form_mat$pathString,long_form_mat[,index_string[i]], sep = "layer_sep")
    }
  }

  long_form_mat$pathString <- gsub("layer_sepNA", "", long_form_mat$pathString)
  long_form_mat$pathString <- gsub("NAlayer_sep", "", long_form_mat$pathString)
  long_form_mat$Metrics <- paste0(paste0("<br>",long_form_mat$Number_of_Cells, " cells"),
                                  paste0("<br>Median Score:",long_form_mat$Median_Score),
                                  paste0("<br>IQR:",long_form_mat$IQR))
  final_cell_call <- unlist(lapply(strsplit(long_form_mat$pathString, "layer_sep"), function(x)x[length(x)]), use.names = F)
  mydf <- long_form_mat[,c("pathString", "Metrics", "Number_of_Cells")]
  mydf$final_cell_call <- final_cell_call


  mydf <- mydf %>%  dplyr::mutate(tree_level = stringr::str_count(string = pathString, pattern = "layer_sep") + 1,
                                  tree_group = stringr::str_replace(string = pathString, pattern = "layer_sep", replacement = ""),
                                  node_type = "decision_node")
  mydf$Number_of_Cells <- as.numeric(mydf$Number_of_Cells)
  index_no_cells <- which(is.na(mydf$Number_of_Cells))
  if(length(index_no_cells)>0){
    mydf <- mydf[-index_no_cells,]
  }
  mytree <- data.tree::as.Node(mydf, pathDelimiter = "layer_sep")

  nodes <- mytree$Get("name")
  color_pal <- character(length = length(nodes))
  names(color_pal) <- nodes
  color_pal_use <- scales::hue_pal()(length(nodes))
  for(i in 1:length(nodes)){
    if(nodes[i] == "All cells"){
      color_pal[which(names(color_pal) == nodes[i])] <- "#FFFFFF"
    }else{
      color_pal[which(names(color_pal) == nodes[i])] <- color_pal_use[i]
    }
  }
  mytree$Set(fill = color_pal)
  if(interactive_mode == T){
    if(ncell_size == T){
      nodeSize = "Number_of_Cells"
    } else{
      nodeSize = NULL

    }
    if(collapsed == T){
      tree_output <- collapsibleTree::collapsibleTree( mytree, linkLength = linkLength, tooltip = T,fillByLevel = T, fontSize = fontSize,   collapsed = T, attribute =c("Metrics"), nodeSize = nodeSize, fill = "fill",
                                                       width =width, height= height )
    } else{
      tree_output <- collapsibleTree::collapsibleTree( mytree, linkLength = linkLength, tooltip = T,fillByLevel = T, fontSize = fontSize,   collapsed = F, attribute =c("Metrics"), nodeSize = nodeSize, fill = "fill",
                                                       width =width, height= height)

    }


  } else if(interactive_mode == F){
    my_tree_non_interactive <- mytree



    data.tree::SetGraphStyle(my_tree_non_interactive, rankdir = "TB")
    data.tree::SetEdgeStyle(my_tree_non_interactive, arrowhead = "vee", color = "blue", penwidth = 2)
    #per default, Node style attributes will be inherited:




    nodes_layer_1 <- names(which(lapply(my_tree_non_interactive, class ) == "list"))
    for(layer_1 in 1:length(nodes_layer_1)){
      data.tree::SetNodeStyle(my_tree_non_interactive[[nodes_layer_1[layer_1]]],
                              color = "black",
                              label = paste0(my_tree_non_interactive[[nodes_layer_1[layer_1]]]$name,gsub("<br>", "\n", my_tree_non_interactive[[nodes_layer_1[layer_1]]]$Metrics)),
                              penwidth = 3,
                              fontcolor = "black", fontsize = 45, inherit = T)
      if(length(my_tree_non_interactive[[nodes_layer_1[layer_1]]]$children)>0){
        nodes_layer_2 <- names(which(lapply(my_tree_non_interactive[[nodes_layer_1[layer_1]]], class ) == "list"))
        for(layer_2 in 1:length(nodes_layer_2)){
          data.tree::SetNodeStyle(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]],
                                  color = "black",
                                  label = paste0(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]]$name, gsub("<br>", "\n", my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]]$Metrics)),
                                  penwidth = 3,
                                  fontcolor = "black", fontsize = 45, inherit = T)
          if(length(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]]$children)>0){
            nodes_layer_3 <- names(which(lapply(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]], class ) == "list"))
            for(layer_3 in 1:length(nodes_layer_3)){
              data.tree::SetNodeStyle(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]],
                                      color = "black",
                                      label = paste0(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]]$name,gsub("<br>", "\n", my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]]$Metrics)),
                                      penwidth = 3,
                                      fontcolor = "black", fontsize = 45, inherit = T)

              if(length(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]]$children)>0){
                nodes_layer_4 <- names(which(lapply(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]], class ) == "list"))
                for(layer_4 in 1:length(nodes_layer_4)){
                  data.tree::SetNodeStyle(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]],
                                          color = "black",
                                          label = paste0(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]]$name,gsub("<br>", "\n", my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]]$Metrics)),
                                          penwidth = 3,
                                          fontcolor = "black", fontsize = 45, inherit = T)
                  if(length(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]]$children)>0){
                    nodes_layer_5 <- names(which(lapply(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]], class ) == "list"))
                    for(layer_5 in 1:length(nodes_layer_5)){
                      data.tree::SetNodeStyle(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]][[nodes_layer_5[layer_5]]],
                                              color = "black",
                                              label = paste0(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]][[nodes_layer_5[layer_5]]]$name,gsub("<br>", "\n", my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]][[nodes_layer_5[layer_5]]]$Metrics)),
                                              penwidth = 3,
                                              fontcolor = "black", fontsize = 45, inherit = T)
                      if(length(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]][[nodes_layer_5[layer_5]]]$children)>0){
                        nodes_layer_6 <- names(which(lapply(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]][[nodes_layer_5[layer_5]]], class ) == "list"))
                        for(layer_6 in 1:length(nodes_layer_6)){
                          data.tree::SetNodeStyle(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]][[nodes_layer_5[layer_5]]][[nodes_layer_6[layer_6]]],
                                                  color = "black",
                                                  label = paste0(my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]][[nodes_layer_5[layer_5]]][[nodes_layer_6[layer_6]]]$name,gsub("<br>", "\n", my_tree_non_interactive[[nodes_layer_1[layer_1]]][[nodes_layer_2[layer_2]]][[nodes_layer_3[layer_3]]][[nodes_layer_4[layer_4]]][[nodes_layer_5[layer_5]]][[nodes_layer_6[layer_6]]]$Metrics)),
                                                  penwidth = 3,
                                                  fontcolor = "black", fontsize = 45, inherit = T)
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }


    data.tree::SetNodeStyle(my_tree_non_interactive, style = "filled,rounded", shape = "box", fillcolor = color_pal,
                            fontname = "helvetica", tooltip = data.tree::GetDefaultTooltip, fontcolor = "black", fontsize = 45,label = paste0(my_tree_non_interactive$root$name,gsub("<br>", "\n", my_tree_non_interactive$root$Metrics)),
    )
    data.tree::SetGraphStyle(my_tree_non_interactive, rankdir = "TB")
    tree_output <- plot(my_tree_non_interactive)
  }

  if(save_results == T){
    htmlwidgets::saveWidget(tree_output,file = paste0(save_dir,"/", project_name, "_annotation_results_tree.html"), libdir = NULL)
  }
  return(tree_output)
}

