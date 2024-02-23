#' create_summary_matrix_pan_can_v3
#'
#' @param prediction_list output from run_classifier()
#' @param use_CNVs logical whether to use CNVs inference to validate cancer cluster - note takes long time
#' @param modify_results logical whether to assume one cancer in the sample and reassign minority cancer calls as Unknown cells
#' @param mc.cores number of cores to use in parallel analysis
#' @param raw_counts gene by cell count matrix
#' @param min_prop minimum % of major cancer class for modifications
#' @param breast_mode run breast subclassfication
#' @param fine_grained_T run T cell subclassfication
#' @param confidence_cutoff logical to use confidence cutoffs
#' @param pan_cancer logical to use pan cancer features or cancer specific features
#' @param cancer_confidence assesing whether a cancer prediction should be interpreted as confident. by default we use layer specific cutoff derived from our paper that optimally split true positives and true negatives. Set to a numeric from 0-1 to set a hard cutoff for all layers.
#' @param normal_tissue Set to true if you know the sample is only normal cells, will convert all cancer annotations to normal tissue cell
#' @param low_res_mode Set to true to convert all high res subtypes into low res cell types
#' @param known_cancer_type converts all cancer cell type annotations to provided string, recommended if one knows with 100% certainty what the cancer type is.
#'
#' @return dataframe of cell names and final prediction
#' @export
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
#' summary_masterrix <- create_summary_matrix(prediction_list = cell_predictions, use_CNVs = F, modify_results = T, mc.cores = 1, raw_counts = lung_cancer_demo_data, min_prop = 0.5 )
#' #create results matrix with CNV corrections
#' results_lung_CNV <- create_summary_matrix(prediction_list = cell_predictions, use_CNVs = T, modify_results = T, mc.cores = 6, raw_counts = lung_cancer_demo_data, min_prop = 0.5 )
#' }
create_summary_matrix <- function(raw_counts, prediction_list, use_CNVs = FALSE, modify_results = TRUE, mc.cores = 1,
                                  min_prop = 0.5, breast_mode = F, fine_grained_T = T, confidence_cutoff = T, pan_cancer = F,
                                  cancer_confidence = "default", normal_tissue = F, low_res_mode = FALSE,
                                  known_cancer_type = NULL){
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
  }
  if(confidence_cutoff == T){
    summary_master <- data.frame(row.names(prediction_list[["layer_1"]]),stringsAsFactors = F)
    colnames(summary_master) <- c("cell_names")
    layer_1 <- prediction_list[["layer_1"]]$predicted_tissue_with_cutoff

    summary_master <-  cbind(summary_master, layer_1)
    layer_2 <- c()
    median_score_class_layer_1 <- c()
    for (i in 1:nrow(summary_master)){
      if(summary_master[i, "layer_1"] %in% c("HSPC","B cell","CD4+ T cell", "CD8+ T cell", "macrophage_or_dendritic_cell","Mast cell",
                                             "Natural killer cell","Blood_Cell")){

        median_score_class_layer_1[i] <- median(as.numeric(prediction_list[["layer_1"]][which(prediction_list[["layer_1"]]$predicted_tissue_with_cutoff == "Blood_Cell"),"blood_score"]))
        layer_2[i] <- as.character(prediction_list[["layer_2_blood"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

      } else if(summary_master[i, "layer_1"] %in%
                c( "Bile Duct Cancer","Bladder Cancer",  "Bone Cancer",     "Brain Cancer",    "Breast Cancer",
                   "Colon/Colorectal Cancer","Endometrial/Uterine Cancer","Esophageal Cancer",
                   "Gallbladder Cancer", "Gastric Cancer", "Glial Cells",  "Kidney Cancer", "Liver Cancer", "Lung Cancer",
                   "Neuroblastoma",  "Oligodendrocytes","Ovarian Cancer",
                   "Pancreatic Cancer",
                   "Prostate Cancer", "Skin Cancer",  "Endothelial Cells", "Fibroblasts", "Myofibroblasts",  "Smooth Muscle Cells",
                   "Sarcoma", "Tissue_Cell_Normal_or_Cancer")){
        layer_2[i] <- as.character(prediction_list[["layer_2_non_blood"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_1[i] <- median(as.numeric(prediction_list[["layer_1"]][which(prediction_list[["layer_1"]]$predicted_tissue_with_cutoff == "Tissue_Cell_Normal_or_Cancer"),"cancer_normal_stromal_score"]))


      } else {
        layer_2[i] <- as.character(layer_1[i])
        median_score_class_layer_1[i] <- NA
      }
    }
    summary_master <-  cbind(summary_master, layer_2)
    median_score_class_layer_2 <- c()

    layer_3 <- c()
    for (i in 1:nrow(summary_master)){
      if(summary_master[i, "layer_2"] %in% c("CD4+ T cell", "CD8+ T cell", "Natural killer cell", "T_or_NK_lymphocyte")){
        layer_3[i] <- as.character(prediction_list[["layer_3_TNK"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_blood"]][which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff == "T_or_NK_lymphocyte"),"non_B_lymphocyte_score"]))


      } else if(summary_master[i, "layer_2"] %in% c("macrophage_or_dendritic_cell", "macrophage_DC_score")){
        layer_3[i] <- as.character(prediction_list[["layer_3_myeloid"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_blood"]][which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff == "macrophage_or_dendritic_cell"),"macrophage_DC_score"]))

      } else if(summary_master[i, "layer_2"] %in% c("B cell", "B_cell_score", "B cell or Plasmablast")){
        layer_3[i] <- as.character(prediction_list[["layer_3_BCell"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_blood"]][which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff == "B cell or Plasmablast"),"B_cell_score"]))


      } else if(summary_master[i, "layer_2"] %in% c(
        "Bile Duct Cancer","Bladder Cancer",  "Bone Cancer","Brain Cancer","Breast Cancer",
        "Colon/Colorectal Cancer","Endometrial/Uterine Cancer","Esophageal Cancer",
        "Gallbladder Cancer", "Gastric Cancer",   "Kidney Cancer", "Liver Cancer", "Lung Cancer",
        "Neuroblastoma",  "Ovarian Cancer",
        "Pancreatic Cancer",
        "Prostate Cancer", "Skin Cancer",
        "Sarcoma", "Non Stromal Cell")){
        layer_3[i] <- as.character(prediction_list[["layer_3_non_stromal"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_non_blood"]][which(prediction_list[["layer_2_non_blood"]]$predicted_tissue_with_cutoff == "Non Stromal Cell"),"non_stromal_score"]))

      } else if(summary_master[i, "layer_2"] %in%
                c("Stromal Cell","Fibroblasts", "Myofibroblasts","Smooth Muscle Cells","Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts",
                  "fibroblasts_score")){
        layer_3[i] <- as.character(prediction_list[["layer_3_stromal"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_non_blood"]][which(prediction_list[["layer_2_non_blood"]]$predicted_tissue_with_cutoff == "Stromal Cell"),"stromal_score"]))
      } else if(summary_master[i, "layer_2"] == "Cell low quality"){
        summary_master[i, "layer_2"] <- summary_master[i, "layer_1"]
        layer_3[i] <- summary_master[i, "layer_1"]
        median_score_class_layer_2[i] <- NA

      } else {
        layer_3[i] <- as.character(layer_2[i])
        if(layer_2[i] %in% c("Mast cell")){
          median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_blood"]][which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff == "Mast cell"),"mast_score"]))
        } else if(layer_2[i] %in% c("HSPC", "HSPC_score")){
          median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_blood"]][which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff %in% c("HSPC", "HSPC_score")),"HSPC_score"]))
        } else if(layer_2[i] %in% c("Glial Cells")){
          median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_non_blood"]][which(prediction_list[["layer_2_non_blood"]]$predicted_tissue_with_cutoff %in% c("Glial Cells")),"Glial Cells"]))
        } else if(layer_2[i] %in% c("Oligodendrocytes")){
          median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_non_blood"]][which(prediction_list[["layer_2_non_blood"]]$predicted_tissue_with_cutoff %in% c("Oligodendrocytes")),"Oligodendrocytes"]))
        } else{
          median_score_class_layer_2[i] <- median_score_class_layer_1[i]
        }
      }
    }
    summary_master <-  cbind(summary_master, layer_3)
    layer_4 <- c()
    median_score_class_layer_3 <- c()

    for (i in 1:nrow(summary_master)){
      if(summary_master[i, "layer_3"]  %in% c("CD4+ T cell" ,"CD4 or CD8 T cell")){
        layer_4[i] <- as.character(prediction_list[["layer_4_CD4_CD8"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_TNK"]][which(prediction_list[["layer_3_TNK"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell" ,"CD4 or CD8 T cell")),"CD4_CD8_score"]))


      } else if (summary_master[i, "layer_3"] %in% c("Natural killer cell", "NK or CD8 T cell")) {
        layer_4[i] <- as.character(prediction_list[["layer_4_CD8_NK"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_TNK"]][which(prediction_list[["layer_3_TNK"]]$predicted_tissue_with_cutoff %in% c("Natural killer cell", "NK or CD8 T cell")),"NK_CD8_score"]))

      } else if (summary_master[i, "layer_3"] == "Dendritic Cell") {
        layer_4[i] <- as.character(prediction_list[["layer_4_dendritic"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_myeloid"]][which(prediction_list[["layer_3_myeloid"]]$predicted_tissue_with_cutoff %in% c("Dendritic Cell")),"DC_score"]))

      } else if (summary_master[i, "layer_3"] == "Macrophage or Monocyte") {
        layer_4[i] <- as.character(prediction_list[["layer_4_macrophage"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_myeloid"]][which(prediction_list[["layer_3_myeloid"]]$predicted_tissue_with_cutoff %in% c("Macrophage or Monocyte")),"macrophage_score"]))

      } else if (summary_master[i, "layer_3"] %in%c(
        "Breast Cancer","Endometrial/Uterine Cancer", "Lung Cancer", "Ovarian Cancer",
        "Prostate Cancer", "Kidney Cancer",  "Non GI Epithelial Cell")) {
        layer_4[i] <- as.character(prediction_list[["layer_4_non_GI"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_non_stromal"]][which(prediction_list[["layer_3_non_stromal"]]$predicted_tissue_with_cutoff %in% c("Breast Cancer","Endometrial/Uterine Cancer", "Lung Cancer", "Ovarian Cancer",
                                                                                                                                                                                       "Prostate Cancer", "Kidney Cancer",  "Non GI Epithelial Cell")),"non_GI_score"]))

      } else if (summary_master[i, "layer_3"] %in% c("Bile Duct Cancer","Bladder Cancer",
                                                     "Colon/Colorectal Cancer","Esophageal Cancer",
                                                     "Gallbladder Cancer", "Gastric Cancer", "Liver Cancer",
                                                     "Pancreatic Cancer","GI Epithelial Cell")) {
        layer_4[i] <- as.character(prediction_list[["layer_4_GI"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_non_stromal"]][which(prediction_list[["layer_3_non_stromal"]]$predicted_tissue_with_cutoff %in% c("Bile Duct Cancer","Bladder Cancer",
                                                                                                                                                                                       "Colon/Colorectal Cancer","Esophageal Cancer",
                                                                                                                                                                                       "Gallbladder Cancer", "Gastric Cancer", "Liver Cancer",
                                                                                                                                                                                       "Pancreatic Cancer","GI Epithelial Cell")),"GI_score"]))

      } else if (summary_master[i, "layer_3"] %in% c("Bone Cancer","Brain Cancer",
                                                     "Neuroblastoma", "Skin Cancer",
                                                     "Sarcoma", "Soft Tissue or Neuro Cancer Cell")) {
        layer_4[i] <- as.character(prediction_list[["layer_4_soft_tissue_neuro"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_non_stromal"]][which(prediction_list[["layer_3_non_stromal"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer","Brain Cancer",
                                                                                                                                                                                       "Neuroblastoma", "Skin Cancer",
                                                                                                                                                                                       "Sarcoma", "Soft Tissue or Neuro Cancer Cell")),"soft_tissue_neuro_score"]))



      } else if(summary_master[i, "layer_3"] %in% c("Cancer Associated Fibroblasts")){
        layer_4[i] <- as.character(prediction_list[["layer_4_CAF"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_stromal"]][which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff %in% c("Cancer Associated Fibroblasts")),"Cancer Associated Fibroblasts"]))

      } else if(summary_master[i, "layer_3"] == "Cell low quality"){
        summary_master[i, "layer_3"] <- summary_master[i, "layer_2"]
        layer_4[i] <- summary_master[i, "layer_2"]
        median_score_class_layer_3[i] <- NA

      } else {
        layer_4[i] <- as.character(layer_3[i])
        if(layer_3[i] %in% c("B Cell")){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_BCell"]][which(prediction_list[["layer_3_BCell"]]$predicted_tissue_with_cutoff == "B Cell"),"B_cell_score"]))
        } else if(layer_3[i] %in% c("Plasmablast")){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_BCell"]][which(prediction_list[["layer_3_BCell"]]$predicted_tissue_with_cutoff == "Plasmablast"),"plasmablast_score"]))
        } else if(layer_3[i] %in% c("Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts" )){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_stromal"]][which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff %in% c("Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts")),"cancer_associated_fibroblast_score"]))
        } else if(layer_3[i] %in% c("Smooth Muscle Cells" )){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_stromal"]][which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff %in% c("Smooth Muscle Cells")),"Smooth Muscle Cells"]))
        } else if(layer_3[i] %in% c("Fibroblasts", "Myofibroblasts")){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_stromal"]][which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff %in% c("Fibroblasts", "Myofibroblasts")),"fibroblasts_score"]))
        } else if(layer_3[i] %in% c("Endothelial Cells")){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_stromal"]][which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff %in% c("Endothelial Cells")),"Endothelial Cells"]))
        } else{
          median_score_class_layer_3[i] <- median_score_class_layer_2[i]
        }


      }
    }
    summary_master <-  cbind(summary_master, layer_4)
    median_score_class_layer_4 <- c()
    layer_5 <- c()
    for (i in 1:nrow(summary_master)){
      if(summary_master[i, "layer_4"] %in%c(
        "Endometrial/Uterine Cancer", "Ovarian Cancer", "Kidney Cancer", "Ovarian/Endometrial/Kidney Cell")){
        layer_5[i] <- as.character(prediction_list[["layer_5_ov_endo_kid"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_non_GI"]][which(prediction_list[["layer_4_non_GI"]]$predicted_tissue_with_cutoff %in% c("Endometrial/Uterine Cancer", "Ovarian Cancer", "Kidney Cancer", "Ovarian/Endometrial/Kidney Cell")),"ov_endo_kid_score"]))

      } else if(summary_master[i, "layer_4"] %in%c(
        "Lung Cancer") & summary_master[i, "layer_3"] %in%c(
          "Soft Tissue or Neuro Cancer Cell")){
        layer_5[i] <- "Lung Cancer"
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_soft_tissue_neuro"]][which(prediction_list[["layer_4_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Lung Cancer")),"Lung Cancer"]))

      } else if(summary_master[i, "layer_4"] %in%c(
        "Breast Cancer","Lung Cancer",
        "Prostate Cancer",  "Breast/Lung/Prostate Cell")){
        layer_5[i] <- as.character(prediction_list[["layer_5_breast_lung_prostate"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_non_GI"]][which(prediction_list[["layer_4_non_GI"]]$predicted_tissue_with_cutoff %in% c("Breast Cancer","Lung Cancer",
                                                                                                                                                                             "Prostate Cancer",  "Breast/Lung/Prostate Cell")),"breast_lung_prostate_score"]))

      } else if(summary_master[i, "layer_4"] %in% c("Bile Duct Cancer","Bladder Cancer",
                                                    "Gallbladder Cancer",  "Liver Cancer",
                                                    "Pancreatic Cancer", "Billiary Cell")){
        layer_5[i] <- as.character(prediction_list[["layer_5_biliary"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_GI"]][which(prediction_list[["layer_4_GI"]]$predicted_tissue_with_cutoff %in% c("Bile Duct Cancer","Bladder Cancer",
                                                                                                                                                                     "Gallbladder Cancer",  "Liver Cancer",
                                                                                                                                                                     "Pancreatic Cancer", "Billiary Cell")),"billiary_score"]))


      } else if(summary_master[i, "layer_4"] %in% c("Colon/Colorectal Cancer","Esophageal Cancer",
                                                    "Gastric Cancer","Colorectal/Esophageal/Gastric Cell")){
        layer_5[i] <- as.character(prediction_list[["layer_5_digestive"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_GI"]][which(prediction_list[["layer_4_GI"]]$predicted_tissue_with_cutoff %in% c("Colon/Colorectal Cancer","Esophageal Cancer",
                                                                                                                                                                     "Gastric Cancer","Colorectal/Esophageal/Gastric Cell")),"digestive_score"]))
      } else if(summary_master[i, "layer_4"] %in% c("Bone Cancer","Brain Cancer",
                                                    "Neuroblastoma",
                                                    "Sarcoma", "Soft Tissue or Neuro Cancer Cell", "Soft Tissue or Neuro Cancer Cell")){
        layer_5[i] <- as.character(prediction_list[["layer_5_soft_tissue_neuro"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_soft_tissue_neuro"]][which(prediction_list[["layer_4_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer","Brain Cancer",
                                                                                                                                                                                                   "Neuroblastoma",
                                                                                                                                                                                                   "Sarcoma", "Soft Tissue or Neuro Cancer Cell", "Soft Tissue or Neuro Cancer Cell")),"soft_tissue_neuro_score"]))

      } else if(summary_master[i, "layer_4"] %in% c("cDC")){
        layer_5[i] <- as.character(prediction_list[["layer_5_cDC"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_dendritic"]][which(prediction_list[["layer_4_dendritic"]]$predicted_tissue_with_cutoff %in% c("cDC")),"cDC"]))

      } else if(summary_master[i, "layer_4"] %in% c("Macrophage")){
        layer_5[i] <- as.character(prediction_list[["layer_5_macrophage"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_macrophage"]][which(prediction_list[["layer_4_macrophage"]]$predicted_tissue_with_cutoff %in% c("Macrophage")),"Macrophage"]))

      } else if(summary_master[i, "layer_4"] %in% c("Monocyte")){
        layer_5[i] <- as.character(prediction_list[["layer_5_monocyte"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_macrophage"]][which(prediction_list[["layer_4_macrophage"]]$predicted_tissue_with_cutoff %in% c("Monocyte")),"Monocyte"]))

      } else if(summary_master[i, "layer_4"] %in% c("CD4+ T cell")){
        if(fine_grained_T == T){
          layer_5[i] <- as.character(prediction_list[["layer_5_CD4"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell")),"CD4+ T cell"]))

        } else{
          layer_5[i] <- as.character(layer_4[i])
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell")),"CD4+ T cell"]))
        }
      } else if(summary_master[i, "layer_4"] %in% c("CD8+ T cell")){
        if(fine_grained_T == T){
          layer_5[i] <- as.character(prediction_list[["layer_5_CD8"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

          if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "CD4 or CD8 T cell"){
            median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
          } else if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "NK or CD8 T cell"){
            median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD8_NK"]][which(prediction_list[["layer_4_CD8_NK"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
          }
        } else{
          layer_5[i] <- as.character(layer_4[i])
          if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "CD4 or CD8 T cell"){
            median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
          } else if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "NK or CD8 T cell"){
            median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD8_NK"]][which(prediction_list[["layer_4_CD8_NK"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
          }
        }
      } else if(summary_master[i, "layer_4"] == "Cell low quality"){
        summary_master[i, "layer_4"] <- summary_master[i, "layer_3"]
        median_score_class_layer_4[i] <- NA
        layer_5[i] <- summary_master[i, "layer_3"]
      } else {
        layer_5[i] <- as.character(layer_4[i])
        if(layer_4[i] %in% c("Skin Cancer")){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_soft_tissue_neuro"]][which(prediction_list[["layer_4_soft_tissue_neuro"]]$predicted_tissue_with_cutoff == "Skin Cancer"),"Skin Cancer"]))
        } else if(layer_4[i] %in% c("cDC")){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_dendritic"]][which(prediction_list[["layer_4_dendritic"]]$predicted_tissue_with_cutoff == "cDC"),"cDC"]))
        } else if(layer_4[i] %in% c("pDC")){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_dendritic"]][which(prediction_list[["layer_4_dendritic"]]$predicted_tissue_with_cutoff == "pDC"),"pDC"]))
        } else if(layer_4[i] %in% c("ASDC")){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_dendritic"]][which(prediction_list[["layer_4_dendritic"]]$predicted_tissue_with_cutoff == "ASDC"),"ASDC"]))
        } else if(layer_4[i] %in% c("Macrophage" )){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_macrophage"]][which(prediction_list[["layer_4_macrophage"]]$predicted_tissue_with_cutoff == "Macrophage"),"Macrophage"]))
        } else if(layer_4[i] %in% c("Monocyte" )){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_macrophage"]][which(prediction_list[["layer_4_macrophage"]]$predicted_tissue_with_cutoff == "Monocyte"),"Monocyte"]))
        } else if(layer_4[i] %in% c("CD4+ T cell") ){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell")),"CD4+ T cell"]))
        } else if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "CD4 or CD8 T cell"){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
        } else if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "NK or CD8 T cell"){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD8_NK"]][which(prediction_list[["layer_4_CD8_NK"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
        } else if(layer_4[i] %in% c("Natural killer cell")){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD8_NK"]][which(prediction_list[["layer_4_CD8_NK"]]$predicted_tissue_with_cutoff %in% c("Natural killer cell")),"Natural killer cell"]))
        } else{
          median_score_class_layer_4[i] <- median_score_class_layer_3[i]
        }
      }
    }


    summary_master <-  cbind(summary_master, layer_5)

    layer_6 <- c()
    median_score_class_layer_5 <- c()
    for (i in 1:nrow(summary_master)){
      if(summary_master[i, "layer_5"] %in% c("Brain Cancer",
                                             "Neuroblastoma", "Brain/Neuroblastoma Cancer Cell")){
        layer_6[i] <- as.character(prediction_list[["layer_6_brain_nbm"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_soft_tissue_neuro"]][which(prediction_list[["layer_5_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Brain Cancer",
                                                                                                                                                                                                   "Neuroblastoma", "Brain/Neuroblastoma Cancer Cell")),"brain_nbm_score"]))

      } else if(summary_master[i, "layer_5"] %in% c("Bone Cancer",
                                                    "Sarcoma", "Soft Tissue Cancer Cell")){
        layer_6[i] <- as.character(prediction_list[["layer_6_soft_tissue"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_soft_tissue_neuro"]][which(prediction_list[["layer_5_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer",
                                                                                                                                                                                                   "Sarcoma", "Soft Tissue Cancer Cell")),"soft_tissue_score"]))


      } else if(summary_master[i, "layer_5"] %in% c("Breast Cancer Cell")){
        if(breast_mode == T){
          layer_6[i] <- as.character(prediction_list[["layer_6_breast"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_breast_lung_prostate"]][which(prediction_list[["layer_5_breast_lung_prostate"]]$predicted_tissue_with_cutoff %in% c("Breast Cancer Cell")),"Breast Cancer"]))

        } else{
          layer_6[i] <- as.character(layer_5[i])
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_breast_lung_prostate"]][which(prediction_list[["layer_5_breast_lung_prostate"]]$predicted_tissue_with_cutoff %in% c("Breast Cancer Cell")),"Breast Cancer"]))
        }
      } else if(summary_master[i, "layer_5"] == "Cell low quality"){
        summary_master[i, "layer_5"] <- summary_master[i, "layer_4"]
        layer_6[i] <- summary_master[i, "layer_4"]
        median_score_class_layer_5[i] <- NA
      } else {
        layer_6[i] <- as.character(layer_5[i])
        if(layer_5[i] == "Lung Cancer Cell" & layer_4[i] == "Breast/Lung/Prostate Cell"){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_breast_lung_prostate"]][which(prediction_list[["layer_5_breast_lung_prostate"]]$predicted_tissue_with_cutoff %in% c("Lung Cancer Cell")),"Lung Cancer"]))
        } else if(layer_5[i] == "Prostate Cancer Cell" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_breast_lung_prostate"]][which(prediction_list[["layer_5_breast_lung_prostate"]]$predicted_tissue_with_cutoff %in% c("Prostate Cancer Cell")),"Prostate Cancer"]))
        } else if(layer_5[i] == "Ovarian Cancer Cell" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_ov_endo_kid"]][which(prediction_list[["layer_5_ov_endo_kid"]]$predicted_tissue_with_cutoff %in% c("Ovarian Cancer Cell")),"Ovarian Cancer"]))
        } else if(layer_5[i] == "Kidney Cancer Cell" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_ov_endo_kid"]][which(prediction_list[["layer_5_ov_endo_kid"]]$predicted_tissue_with_cutoff %in% c("Kidney Cancer Cell")),"Kidney Cancer"]))
        } else if(layer_5[i] == "Endometrial Cancer Cell" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_ov_endo_kid"]][which(prediction_list[["layer_5_ov_endo_kid"]]$predicted_tissue_with_cutoff %in% c("Endometrial Cancer Cell")),"Endometrial/Uterine Cancer"]))
        } else if(layer_5[i] == "Pancreatic Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_biliary"]][which(prediction_list[["layer_5_biliary"]]$predicted_tissue_with_cutoff %in% c("Pancreatic Cancer")),"Pancreatic Cancer"]))
        } else if(layer_5[i] == "Liver Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_biliary"]][which(prediction_list[["layer_5_biliary"]]$predicted_tissue_with_cutoff %in% c("Liver Cancer")),"Liver Cancer"]))
        } else if(layer_5[i] == "Gallbladder Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_biliary"]][which(prediction_list[["layer_5_biliary"]]$predicted_tissue_with_cutoff %in% c("Gallbladder Cancer")),"Gallbladder Cancer"]))
        } else if(layer_5[i] == "Bladder Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_biliary"]][which(prediction_list[["layer_5_biliary"]]$predicted_tissue_with_cutoff %in% c("Bladder Cancer")),"Bladder Cancer"]))
        } else if(layer_5[i] == "Bile Duct Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_biliary"]][which(prediction_list[["layer_5_biliary"]]$predicted_tissue_with_cutoff %in% c("Bile Duct Cancer")),"Bile Duct Cancer"]))
        } else if(layer_5[i] == "Colon/Colorectal Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_digestive"]][which(prediction_list[["layer_5_digestive"]]$predicted_tissue_with_cutoff %in% c("Colon/Colorectal Cancer")),"Colon/Colorectal Cancer"]))
        } else if(layer_5[i] == "Esophageal Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_digestive"]][which(prediction_list[["layer_5_digestive"]]$predicted_tissue_with_cutoff %in% c("Esophageal Cancer")),"Esophageal Cancer"]))
        } else if(layer_5[i] == "Gastric Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_digestive"]][which(prediction_list[["layer_5_digestive"]]$predicted_tissue_with_cutoff %in% c("Gastric Cancer")),"Gastric Cancer"]))
        } else{
          median_score_class_layer_5[i] <- median_score_class_layer_4[i]
        }
      }
    }

    summary_master <-  cbind(summary_master, layer_6)
    median_score_class_layer_6 <- c()
    for (i in 1:nrow(summary_master)){
      if(layer_6[i] == "TNBC Breast Cancer Cell"){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_breast"]][which(prediction_list[["layer_6_breast"]]$predicted_tissue_with_cutoff %in% c("TNBC Breast Cancer Cell")),"TNBC"]))
      } else if(layer_6[i] == "Her2+ Breast Cancer Cell"){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_breast"]][which(prediction_list[["layer_6_breast"]]$predicted_tissue_with_cutoff %in% c("Her2+ Breast Cancer Cell")),"HER2+"]))
      } else if(layer_6[i] == "ER+ Breast Cancer Cell"){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_breast"]][which(prediction_list[["layer_6_breast"]]$predicted_tissue_with_cutoff %in% c("ER+ Breast Cancer Cell")),"ER+"]))
      } else if(layer_6[i] == "Brain Cancer" ){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_brain_nbm"]][which(prediction_list[["layer_6_brain_nbm"]]$predicted_tissue_with_cutoff %in% c("Brain Cancer")),"Brain Cancer"]))
      } else if(layer_6[i] == "Neuroblastoma" ){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_brain_nbm"]][which(prediction_list[["layer_6_brain_nbm"]]$predicted_tissue_with_cutoff %in% c("Neuroblastoma")),"Neuroblastoma"]))
      } else if(layer_6[i] == "Bone Cancer" ){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_soft_tissue"]][which(prediction_list[["layer_6_soft_tissue"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer")),"Bone Cancer"]))
      } else if(layer_6[i] == "Sarcoma" ){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_soft_tissue"]][which(prediction_list[["layer_6_soft_tissue"]]$predicted_tissue_with_cutoff %in% c("Sarcoma")),"Sarcoma"]))
      } else{
        median_score_class_layer_6[i] <- median_score_class_layer_5[i]
      }
    }
    summary_master <- cbind(summary_master, median_score_class_layer_1,median_score_class_layer_2,median_score_class_layer_3,median_score_class_layer_4,median_score_class_layer_5,median_score_class_layer_6  )
    row.names(summary_master) <- summary_master$cell_names

    summary_master$layer_3 <- gsub("Cancer$", "Cancer Cell", summary_master$layer_3)
    summary_master$layer_4 <- gsub("Cancer$", "Cancer Cell", summary_master$layer_4)
    summary_master$layer_5 <- gsub("Cancer$", "Cancer Cell", summary_master$layer_5)

    summary_master$layer_6 <- gsub("Cancer$", "Cancer Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_B_cell_or_plasmablast", "B cell or Plasmablast", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_blood_cell", "Blood Cell", summary_master$layer_6)

    summary_master$layer_6 <- gsub("unclassified_any_cell", "Any Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_normal_or_cancer_tissue", "Non Blood Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_T_or_NK_cell", "T or NK Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_macrophage_or_DC", "Macrophage or Dendritic Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_non_GI_epithelial_cell", "Epithelial Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_GI_epithelial_cell", "GI Tract Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("^NF$", "Normal Fibroblast", summary_master$layer_6)

    if(modify_results == TRUE){
      scATOMIC_pred <- summary_master$layer_6
      summary_master <- cbind(summary_master, scATOMIC_pred)
      if(use_CNVs == TRUE){
        seurat_object <- Seurat::CreateSeuratObject(raw_counts)
        seurat_object <- Seurat::NormalizeData(seurat_object)
        seurat_object <- Seurat::FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
        all.genes <- rownames(seurat_object)
        seurat_object <- Seurat::ScaleData(seurat_object, features = all.genes)
        seurat_object <- Seurat::RunPCA(seurat_object, features = Seurat::VariableFeatures(object = seurat_object))
        seurat_object <- Seurat::FindNeighbors(seurat_object, dims = 1:50)
        seurat_object <- Seurat::FindClusters(seurat_object, resolution = 0.5)
        summary_master <- cbind(summary_master, seurat_object@meta.data[row.names(summary_master), "seurat_clusters"])
        colnames(summary_master)[ncol(summary_master)] <- "seurat_clusters"

        copy_kat_res <- scATOMIC::copy_kat_no_heatmap(rawmat = raw_counts, summary_matrix = summary_master,
                                                                id.type = "S", cell.line = "no",
                                                                ngene.chr = 5, LOW.DR = 0.05, UP.DR = 0.1,
                                                                win.size = 25, KS.cut = 0.1, sam.name = "", distance = "euclidean",
                                                                n.cores = mc.cores)
        copy_kat_res <- as.data.frame(copy_kat_res[,"copykat.pred"], row.names = row.names(copy_kat_res))
        colnames(copy_kat_res) <- "CNV_status"
        summary_master <- merge(summary_master,copy_kat_res,by="row.names",all.x=TRUE)
        row.names(summary_master) <- summary_master$Row.names
        summary_master <- summary_master[,-1]
      }

      s.genes <- Seurat::cc.genes$s.genes
      g2m.genes <- Seurat::cc.genes$g2m.genes
      seurat_object <- CreateSeuratObject(raw_counts, meta.data = summary_master)
      seurat_object <- Seurat::NormalizeData(seurat_object, verbose = F)
      seurat_object <- Seurat::FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = F)
      all.genes <- rownames(seurat_object)
      #try regressing cell cycle genes
      seurat_object <- Seurat::CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
      seurat_object <- Seurat::ScaleData(seurat_object, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(object = seurat_object))
      if(ncol(seurat_object) < 50){
        seurat_object <- Seurat::RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), verbose = F, npcs = (ncol(seurat_object) - 1))

      } else{
        seurat_object <- Seurat::RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), verbose = F)
      }
      seurat_object <- Seurat::FindNeighbors(seurat_object, dims = 1:10, verbose = F)
      seurat_object <- Seurat::FindClusters(seurat_object, resolution = 0.2, verbose = F)

      final_classes <- model_layer_1$classes
      final_classes <- gsub("Cancer$", "Cancer Cell", final_classes)
      cancer_classes <- final_classes[grep("Cancer Cell|oma$|Breast/Lung/Prostate|Ovarian/Endometrial/Kidney|Colorectal Cancer Cell|Endometrial Cancer Cell|Biliary/Hepatic Cancer Cell", final_classes)]
      cancer_classes <- c(cancer_classes,"Breast/Lung/Prostate", "Ovarian/Endometrial/Kidney", "Endometrial Cancer Cell",
                          "Biliary/Hepatic Cancer Cell", "Brain/Neuroblastoma Cancer Cell", "Digestive Tract Cancer Cell",
                          "Soft Tissue Cancer Cell", "Soft Tissue or Neuro Cancer Cell", "Unclassified Soft Tissue or Neuro Cancer Cell","ER+ Breast Cancer Cell",
                          "HER2+ Breast Cancer Cell","Her2+ Breast Cancer Cell","TNBC Breast Cancer Cell")
      predicted_cancer <- seurat_object@meta.data[which(seurat_object@meta.data$layer_6 %in% cancer_classes),]
      if(nrow(predicted_cancer) > 0){
        frequency <- as.data.frame(table(predicted_cancer$scATOMIC_pred), stringsAsFactors = F)
        if(breast_mode == T){
          frequency <- as.data.frame(table(predicted_cancer$layer_5), stringsAsFactors = F)

        }
        proportion <- frequency$Freq/sum(frequency$Freq)
        frequency <- cbind(frequency, proportion)
        max_proportion <- max(frequency$proportion)
        major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
        if(length(major_cancer) > 1){
          major_cancer <- major_cancer[1]
        }
        #if there is more than 1 seurat cluster for cancer we want to split into normal and cancer
        if(major_cancer == "Bile Duct Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "CHOL") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "CHOL")),"Gene"]
        }else if(major_cancer == "Bladder Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "BLCA") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "BLCA")),"Gene"]
        } else if(major_cancer == "Bone Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "SARC") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "SARC")),"Gene"]
        } else if(major_cancer == "Brain Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "GBM", "LGG") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "GBM", "LGG")),"Gene"]
        } else if(major_cancer %in% c( "Breast Cancer Cell","ER+ Breast Cancer Cell","HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell","Her2+ Breast Cancer Cell"  ) ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "BRCA") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "BRCA")),"Gene"]
        } else if(major_cancer == "Colon/Colorectal Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "COAD") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "COAD")),"Gene"]
        } else if(major_cancer %in% c("Endometrial/Uterine Cancer Cell", "Endometrial Cancer Cell" )){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "UCS", "UCEC") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "UCS", "UCEC")),"Gene"]
        } else if(major_cancer == "Gallbladder Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "CHOL") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "CHOL")),"Gene"]
        } else if(major_cancer == "Kidney Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "KICH", "KIRC", "KIRP") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "KICH", "KIRC", "KIRP")),"Gene"]
        } else if(major_cancer == "Liver Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "LIHC") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "LIHC")),"Gene"]
        } else if(major_cancer == "Lung Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "LUAD") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "LUAD")),"Gene"]
        } else if(major_cancer == "Ovarian Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "OV") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "OV")),"Gene"]
        } else if(major_cancer == "Pancreatic Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "PAAD") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "PAAD")),"Gene"]
        } else if(major_cancer == "Prostate Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "PRAD") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "PRAD")),"Gene"]
        } else if(major_cancer == "Sarcoma" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "SARC") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "SARC")),"Gene"]
        } else if(major_cancer == "Skin Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "SKCM") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "SKCM")),"Gene"]
        } else{
          cancer_specific_upreg <- pan_cancer_upreg
          cancer_specific_downreg <- pan_cancer_downreg
        }
        if(pan_cancer == T){
          cancer_specific_upreg <- pan_cancer_upreg
          cancer_specific_downreg <- pan_cancer_downreg
        }
        if(max_proportion < min_prop){
          cancer_specific_upreg <- pan_cancer_upreg
          cancer_specific_downreg <- pan_cancer_downreg
        }
        genes_upreg_shared <- intersect(row.names(raw_counts), cancer_specific_upreg)
        genes_downreg_shared <- intersect(row.names(raw_counts), cancer_specific_downreg)
        index_cancer <- row.names(predicted_cancer)
        count_matrix_up_reg <- raw_counts[c(genes_upreg_shared),index_cancer]
        if(length(index_cancer) == 1){
          count_matrix_up_reg <- matrix(count_matrix_up_reg, nrow = length(count_matrix_up_reg), ncol = length(index_cancer), dimnames = list(names(count_matrix_up_reg), index_cancer))
        }
        frequency_per_row_nozero <- apply(count_matrix_up_reg,1,function(x){length(which(x!=0))})/ncol(count_matrix_up_reg)
        genes_upreg_shared <- names(frequency_per_row_nozero)[which(frequency_per_row_nozero > 0)]
        count_matrix_down_reg <- raw_counts[c(genes_downreg_shared),index_cancer]
        if(length(index_cancer) == 1){
          count_matrix_down_reg <- matrix(count_matrix_down_reg, nrow = length(count_matrix_down_reg), ncol = length(index_cancer), dimnames = list(names(count_matrix_down_reg), index_cancer))
        }
        frequency_per_row_nozero <- apply(count_matrix_down_reg,1,function(x){length(which(x!=0))})/ncol(count_matrix_down_reg)
        genes_downreg_shared <- names(frequency_per_row_nozero)[which(frequency_per_row_nozero > 0)]
        cancer_subset <- seurat_object
        cancer_subset <- Rmagic::magic(cancer_subset, seed = 123)
        #cancer_subset <- Seurat::AddModuleScore(cancer_subset, features = genes_upreg_shared, verbose = F, name = "upreg_genes", assay = "MAGIC_RNA")
        #cancer_subset <- Seurat::AddModuleScore(cancer_subset, features = genes_downreg_shared, verbose = F, name = "downreg_genes", assay = "MAGIC_RNA")
        cancer_subset <- Seurat::AddModuleScore(cancer_subset, features = list(upreg_genes = genes_upreg_shared, downreg_genes = genes_downreg_shared), verbose = F, name = c("upreg_genes", "downreg_genes"), assay = "MAGIC_RNA")

        dist_cancer_signature <- amap::Dist(cancer_subset@meta.data[,c("upreg_genes1", "downreg_genes2")], method = "euclidean", nbproc = mc.cores) # distance matrix
        fit_clusters <- hclust(dist_cancer_signature, method="ward.D2")

        groups <- cutree(fit_clusters, k=2)
        cancer_subset <- Seurat::AddMetaData(cancer_subset, groups, "pan_cancer_cluster")
        index_normal <- unique(c(row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$upreg_genes1 <
                                                                            cancer_subset@meta.data$downreg_genes2)],
                                 row.names(cancer_subset@meta.data)[-which(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes)]))
        index_non_normal <- outersect(row.names(cancer_subset@meta.data), index_normal)
        pct_normal_cluster_1 <- length(intersect(index_normal, row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == 1)]))/length(which(cancer_subset@meta.data$pan_cancer_cluster == 1))

        pct_normal_cluster_2 <- length(intersect(index_normal, row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == 2)]))/length(which(cancer_subset@meta.data$pan_cancer_cluster == 2))
        # DimPlot(cancer_subset, group.by = "pan_cancer_cluster")
        # DimPlot(cancer_subset, group.by = "scATOMIC_pred", label = T) + NoLegend()
        # DimPlot(cancer_subset, cells.highlight =  index_normal)
        if (max(c(pct_normal_cluster_1,pct_normal_cluster_2)) - min(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.1){
          clust_1_cancer <- length(which(cancer_subset@meta.data$pan_cancer_cluster == 1 &
                                           cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes ))/nrow(cancer_subset@meta.data)
          clust_2_cancer <- length(which(cancer_subset@meta.data$pan_cancer_cluster == 2 &
                                           cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes ))/nrow(cancer_subset@meta.data)

          if(clust_1_cancer > clust_2_cancer){
            cancer_clust = "1"

          } else{
            cancer_clust = "2"

          }
          index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes )]

        } else if(pct_normal_cluster_1 > pct_normal_cluster_2){
          cancer_clust = "2"
          if(max(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.5){

            index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal) )]

          } else{
            index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust)]

          }


        } else if (pct_normal_cluster_1 < pct_normal_cluster_2){
          cancer_clust = "1"
          if(max(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.5){
            index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal) )]

          } else{
            index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust)]

          }

        }
        continue = TRUE
        if(length(index_cancer) < 20){
          continue = FALSE
        }
        if(pct_normal_cluster_1 != pct_normal_cluster_2){
          if(min(c(pct_normal_cluster_1,pct_normal_cluster_2)) > 0 & max(pct_normal_cluster_1,pct_normal_cluster_2) > 0.5){
            while(continue == TRUE){
              index_cancer_last <- index_cancer
              cancer_subset <- subset(seurat_object, cells = index_cancer)
              cancer_subset <- Seurat::NormalizeData(cancer_subset, verbose = F)
              cancer_subset <- Seurat::ScaleData(cancer_subset, verbose = F)
              cancer_subset <- Rmagic::magic(cancer_subset, seed = 123)
              cancer_subset <- Seurat::AddModuleScore(cancer_subset, features = list(upreg_genes = genes_upreg_shared, downreg_genes = genes_downreg_shared), verbose = F, name = c("upreg_genes", "downreg_genes"), assay = "MAGIC_RNA")

              dist_cancer_signature <- amap::Dist(cancer_subset@meta.data[,c("upreg_genes1", "downreg_genes2")], method = "euclidean", nbproc = mc.cores) # distance matrix


              #hierarcherical clustering based on these scores
              fit_clusters <- hclust(dist_cancer_signature, method="ward.D2")

              groups <- cutree(fit_clusters, k=2)
              cancer_subset <- Seurat::AddMetaData(cancer_subset, groups, "pan_cancer_cluster")
              # DimPlot(cancer_subset,group.by = "pan_cancer_cluster" )
              index_normal <- unique(c(row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$upreg_genes1 <
                                                                                  cancer_subset@meta.data$downreg_genes2)],
                                       row.names(cancer_subset@meta.data)[-which(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes)]))
              index_non_normal <- outersect(row.names(cancer_subset@meta.data), index_normal)
              # DimPlot(cancer_subset, cells.highlight =  index_normal)

              pct_normal_cluster_1 <- length(intersect(index_normal, row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == 1)]))/length(which(cancer_subset@meta.data$pan_cancer_cluster == 1))

              pct_normal_cluster_2 <- length(intersect(index_normal, row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == 2)]))/length(which(cancer_subset@meta.data$pan_cancer_cluster == 2))
              if(max(pct_normal_cluster_1,pct_normal_cluster_2) > 0.2){
                continue <- TRUE
              } else{
                continue <- FALSE
              }
              if (max(c(pct_normal_cluster_1,pct_normal_cluster_2)) - min(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.25){
                clust_1_cancer <- length(which(cancer_subset@meta.data$pan_cancer_cluster == 1 &
                                                 cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes ))/nrow(cancer_subset@meta.data)
                clust_2_cancer <- length(which(cancer_subset@meta.data$pan_cancer_cluster == 2 &
                                                 cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes ))/nrow(cancer_subset@meta.data)

                if(clust_1_cancer > clust_2_cancer){
                  cancer_clust = "1"

                } else{
                  cancer_clust = "2"

                }
                index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal)]


              } else if (max(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.9 ){
                if(pct_normal_cluster_1 > pct_normal_cluster_2){
                  cancer_clust = "2"
                } else{
                  cancer_clust = "1"

                }
                index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal)]

              }else if(pct_normal_cluster_1 > pct_normal_cluster_2){
                cancer_clust = "2"
                if(max(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.5 ){
                  index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal) )]
                  continue <- TRUE
                } else{
                  index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust)]

                }


              } else if (pct_normal_cluster_1 < pct_normal_cluster_2){
                cancer_clust = "1"
                if( max(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.5 ){
                  index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal) )]
                  continue <- TRUE
                } else{
                  index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust)]

                }

              }
              if(length(index_cancer_last) == length(index_cancer)){
                if(length(levels(as.factor(index_cancer_last == index_cancer)) == TRUE) == 1){
                  if (levels(as.factor(index_cancer_last == index_cancer)) == TRUE){
                    continue = FALSE
                  }
                }
              }


              if(continue == FALSE){
                index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$upreg_genes1 >  cancer_subset@meta.data$downreg_genes2 )]
              }
              if(length(index_cancer) < 20){
                continue = FALSE
              }
            }
            index_normal <- outersect(row.names(seurat_object@meta.data)[which(seurat_object@meta.data$scATOMIC_pred %in% cancer_classes)],
                                      index_cancer)
            seurat_object@meta.data[index_normal,"scATOMIC_pred"] <- "Normal Tissue Cell"
            predicted_cancer_new <- seurat_object@meta.data[which(seurat_object@meta.data$scATOMIC_pred %in%c(cancer_classes,"Breast/Lung/Prostate", "Ovarian/Endometrial/Kidney", "Endometrial Cancer Cell",
                                                                                                              "Biliary/Hepatic Cancer Cell", "Brain/Neuroblastoma Cancer Cell", "Digestive Tract Cancer Cell",
                                                                                                              "Soft Tissue Cancer Cell", "Soft Tissue or Neuro Cancer Cell", "Unclassified Soft Tissue or Neuro Cancer Cell","ER+ Breast Cancer Cell",
                                                                                                              "HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell", "Her2+ Breast Cancer Cell")),]

            if(breast_mode == F){
              frequency <- as.data.frame(table(predicted_cancer_new$scATOMIC_pred), stringsAsFactors = F)
              proportion <- frequency$Freq/sum(frequency$Freq)
              frequency <- cbind(frequency, proportion)
              max_proportion <- max(frequency$proportion)
              major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
              if(length(major_cancer) > 1){
                major_cancer <- major_cancer[1]
              }
              seurat_object@meta.data[row.names(predicted_cancer_new),"scATOMIC_pred"] <- major_cancer
              seurat_object@meta.data[,"pan_cancer_cluster"] <- "Normal"
              seurat_object@meta.data[row.names(predicted_cancer_new),"pan_cancer_cluster"] <- "Cancer"
            } else{
              frequency <- as.data.frame(table(predicted_cancer_new$layer_5), stringsAsFactors = F)
              proportion <- frequency$Freq/sum(frequency$Freq)
              frequency <- cbind(frequency, proportion)
              max_proportion <- max(frequency$proportion)
              major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
              if(length(major_cancer) > 1){
                major_cancer <- major_cancer[1]
              }
              index_non_breast <- row.names(predicted_cancer_new)[-which(predicted_cancer_new$scATOMIC_pred %in% c("ER+ Breast Cancer Cell",
                                                                                                                   "HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell",
                                                                                                                   "Her2+ Breast Cancer Cell", "Breast Cancer Cell") )]
              seurat_object@meta.data[index_non_breast,"scATOMIC_pred"] <- major_cancer
              seurat_object@meta.data[,"pan_cancer_cluster"] <- "Normal"
              seurat_object@meta.data[row.names(predicted_cancer_new),"pan_cancer_cluster"] <- "Cancer"
            }
          }
          else{
            index_cancer <- intersect(index_cancer, outersect(index_cancer, row.names(seurat_object@meta.data)[-which(seurat_object@meta.data$scATOMIC_pred %in% cancer_classes)]))
            index_normal <- outersect(row.names(seurat_object@meta.data)[which(seurat_object@meta.data$scATOMIC_pred %in% cancer_classes)],
                                      index_cancer)
            seurat_object@meta.data[index_normal,"scATOMIC_pred"] <- "Normal Tissue Cell"
            predicted_cancer_new <- seurat_object@meta.data[which(seurat_object@meta.data$scATOMIC_pred %in%c(cancer_classes,"Breast/Lung/Prostate", "Ovarian/Endometrial/Kidney", "Endometrial Cancer Cell",
                                                                                                              "Biliary/Hepatic Cancer Cell", "Brain/Neuroblastoma Cancer Cell", "Digestive Tract Cancer Cell",
                                                                                                              "Soft Tissue Cancer Cell", "Soft Tissue or Neuro Cancer Cell", "Unclassified Soft Tissue or Neuro Cancer Cell","ER+ Breast Cancer Cell",
                                                                                                              "HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell", "Her2+ Breast Cancer Cell")),]

            if(breast_mode == F){
              frequency <- as.data.frame(table(predicted_cancer_new$scATOMIC_pred), stringsAsFactors = F)
              proportion <- frequency$Freq/sum(frequency$Freq)
              frequency <- cbind(frequency, proportion)
              max_proportion <- max(frequency$proportion)
              major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
              if(length(major_cancer) > 1){
                major_cancer <- major_cancer[1]
              }
              seurat_object@meta.data[row.names(predicted_cancer_new),"scATOMIC_pred"] <- major_cancer
              seurat_object@meta.data[,"pan_cancer_cluster"] <- "Normal"
              seurat_object@meta.data[row.names(predicted_cancer_new),"pan_cancer_cluster"] <- "Cancer"
            } else{
              frequency <- as.data.frame(table(predicted_cancer_new$layer_5), stringsAsFactors = F)
              proportion <- frequency$Freq/sum(frequency$Freq)
              frequency <- cbind(frequency, proportion)
              max_proportion <- max(frequency$proportion)
              major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
              if(length(major_cancer) > 1){
                major_cancer <- major_cancer[1]
              }
              index_non_breast <- row.names(predicted_cancer_new)[-which(predicted_cancer_new$scATOMIC_pred %in% c("ER+ Breast Cancer Cell",
                                                                                                                   "HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell",
                                                                                                                   "Her2+ Breast Cancer Cell", "Breast Cancer Cell") )]
              seurat_object@meta.data[index_non_breast,"scATOMIC_pred"] <- major_cancer
              seurat_object@meta.data[,"pan_cancer_cluster"] <- "Normal"
              seurat_object@meta.data[row.names(predicted_cancer_new),"pan_cancer_cluster"] <- "Cancer"
            }

          }
        } else{
          index_cancer <- intersect(index_cancer, outersect(index_cancer, row.names(seurat_object@meta.data)[-which(seurat_object@meta.data$scATOMIC_pred %in% cancer_classes)]))
          index_normal <- outersect(row.names(seurat_object@meta.data)[which(seurat_object@meta.data$scATOMIC_pred %in% cancer_classes)],
                                    index_cancer)
          predicted_cancer_new <- seurat_object@meta.data[which(seurat_object@meta.data$scATOMIC_pred %in%c(cancer_classes,"Breast/Lung/Prostate", "Ovarian/Endometrial/Kidney", "Endometrial Cancer Cell",
                                                                                                            "Biliary/Hepatic Cancer Cell", "Brain/Neuroblastoma Cancer Cell", "Digestive Tract Cancer Cell",
                                                                                                            "Soft Tissue Cancer Cell", "Soft Tissue or Neuro Cancer Cell", "Unclassified Soft Tissue or Neuro Cancer Cell","ER+ Breast Cancer Cell",
                                                                                                            "HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell")),]
          frequency <- as.data.frame(table(predicted_cancer_new$scATOMIC_pred), stringsAsFactors = F)
          proportion <- frequency$Freq/sum(frequency$Freq)
          frequency <- cbind(frequency, proportion)
          max_proportion <- max(frequency$proportion)
          major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
          if(length(major_cancer) > 1){
            major_cancer <- major_cancer[1]
          }
          seurat_object@meta.data[row.names(predicted_cancer_new),"scATOMIC_pred"] <- major_cancer
          seurat_object@meta.data[,"pan_cancer_cluster"] <- "Normal"
          seurat_object@meta.data[row.names(predicted_cancer_new),"pan_cancer_cluster"] <- "Cancer"

          summary_master <- seurat_object@meta.data
        }

        summary_master <- seurat_object@meta.data

      }
      summary_master <- seurat_object@meta.data
      summary_master$classification_confidence <- "confident"
      if(nrow(predicted_cancer) > 0){
        median_score_major_cancer <- min(c(summary_master[which(summary_master$layer_6 == major_cancer), "median_score_class_layer_4"],
                                           summary_master[which(summary_master$layer_6 == major_cancer), "median_score_class_layer_5"],
                                           summary_master[which(summary_master$layer_6 == major_cancer), "median_score_class_layer_6"]))
        if(cancer_confidence == "default"){
          index_non_confident <- row.names(summary_master)[which((summary_master$layer_1 %in% c("Tissue_Cell_Normal_or_Cancer") & summary_master$median_score_class_layer_2 < 0.55 ) |
                                                                   (summary_master$scATOMIC_pred == "Brain Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.6 | summary_master$median_score_class_layer_5 < 0.62 |summary_master$median_score_class_layer_6 < 0.8 ) )|
                                                                   (summary_master$scATOMIC_pred == "Breast Cancer Cell" & (summary_master$median_score_class_layer_5 < 0.55) )|
                                                                   (summary_master$scATOMIC_pred == "Colorectal Cancer Cell" & (summary_master$median_score_class_layer_4 < 0.7 |summary_master$median_score_class_layer_5 < 0.68  ) )|
                                                                   (summary_master$scATOMIC_pred == "Colon/Colorectal Cancer Cell" & (summary_master$median_score_class_layer_4 < 0.7 |summary_master$median_score_class_layer_5 < 0.68  ) )|
                                                                   (summary_master$scATOMIC_pred == "Gastric Cancer Cell" & (summary_master$median_score_class_layer_4 < 0.7 |summary_master$median_score_class_layer_5 < 0.68  ) )|
                                                                   (summary_master$scATOMIC_pred == "Endometrial Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.6|summary_master$median_score_class_layer_4 < 0.7|summary_master$median_score_class_layer_5 < 0.5   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Endometrial/Uterine Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.6|summary_master$median_score_class_layer_4 < 0.7|summary_master$median_score_class_layer_5 < 0.5   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Kidney Cancer Cell" & ( summary_master$median_score_class_layer_4 < 0.7) )|
                                                                   (summary_master$scATOMIC_pred == "Liver Cancer Cell" & (summary_master$median_score_class_layer_5 < 0.69   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Lung Cancer Cell" & (summary_master$median_score_class_layer_4 < 0.59 |summary_master$median_score_class_layer_5 < 0.6    ) ) |
                                                                   (summary_master$scATOMIC_pred == "Neuroblastoma Cell" & (summary_master$median_score_class_layer_6 < 0.55   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Ovarian Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.6|summary_master$median_score_class_layer_4 < 0.7|summary_master$median_score_class_layer_5 < 0.5   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Pancreatic Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.59 |  summary_master$median_score_class_layer_4 < 0.55   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Bile Duct Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.59 |  summary_master$median_score_class_layer_4 < 0.55   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Skin Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.6   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Prostate Cancer Cell" & (summary_master$median_score_class_layer_5 < 0.68 )))]

          if(length(index_non_confident) > 0 ){
            summary_master[index_non_confident, "classification_confidence"] <-"low_confidence"
          }
          if(length(index_non_confident) > 0 ){
            preds_layer_6 <- summary_master[index_non_confident, "layer_6"]
            if(major_cancer %in% preds_layer_6){
              warning("Level of confidence for cancer typing is low")
            }
          }



          index_new_major_cancer <- which(summary_master$pan_cancer_cluster == "Cancer" & summary_master$classification_confidence == "confident")
          if(length(index_new_major_cancer) > 0){
            major_cancer_fixed <- names(sort(table(summary_master[index_new_major_cancer, "layer_6"]), decreasing = T))[1]
            summary_master[which(summary_master$pan_cancer_cluster == "Cancer"), "scATOMIC_pred"] <- major_cancer_fixed
          }
          index_esophageal <- row.names(summary_master)[which(summary_master$scATOMIC_pred == "Esophageal Cancer Cell")]
          if(length(index_esophageal) > 0){
            warning("Cancer was called esophageal, but may be a different squamous cell cancer type")
            summary_master[index_esophageal, "classification_confidence"] <-"low_confidence"
          }



          index_na_confident <- row.names(summary_master)[which(summary_master$scATOMIC_pred %in%  c("Bile Duct Cancer Cell",
                                                                                                     "Bone Cancer Cell",
                                                                                                     "Endometrial/Uterine Cancer Cell",
                                                                                                     "Endometrial Cancer Cell",
                                                                                                     "Gallbladder Cancer Cell",
                                                                                                     "Gastric Cancer Cell")) ]
          if(length(index_na_confident) > 0){
            warning("Cancer was not externally validated, no confidence assigned")
            summary_master[index_na_confident, "classification_confidence"] <- NA
          }

        } else if(is.numeric(cancer_confidence) & median_score_major_cancer < cancer_confidence){
          warning(paste0("Level of confidence for cancer typing = ", median_score_major_cancer, "\nThis classfication is not reliable"))
          summary_master[which(summary_master$scATOMIC_pred == major_cancer), "classification_confidence"] <-"low_confidence"
        }

      }
      if(cancer_confidence == "default"){
        index_non_confident_blood_stromal <- row.names(summary_master)[which(
          (summary_master$scATOMIC_pred == "B cell" & summary_master$median_score_class_layer_3 < 0.7) |
            (summary_master$layer_4 == "CD4+ T cell" & (summary_master$median_score_class_layer_2 < 0.6 | summary_master$median_score_class_layer_3 < 0.9 |summary_master$median_score_class_layer_4 < 0.9 ) )|
            (summary_master$layer_4 == "CD8+ T cell" & (summary_master$median_score_class_layer_2 < 0.6 | summary_master$median_score_class_layer_3 < 0.9 |summary_master$median_score_class_layer_4 < 0.6 ) )|
            (summary_master$scATOMIC_pred == "cDC" & (summary_master$median_score_class_layer_2 < 0.9 ) )|
            (summary_master$scATOMIC_pred == "Macrophage" & (summary_master$median_score_class_layer_2 < 0.9 | summary_master$median_score_class_layer_3 < 0.66|summary_master$median_score_class_layer_4 < 0.75 ) ) |
            (summary_master$scATOMIC_pred == "Natural killer cell" & ( summary_master$median_score_class_layer_3 < 0.9) )|
            (summary_master$scATOMIC_pred == "pDC" & (summary_master$median_score_class_layer_2 < 0.9   ) ) |
            (summary_master$scATOMIC_pred == "Plasmablast" & (summary_master$median_score_class_layer_2 < 0.8 | summary_master$median_score_class_layer_3 < 0.9 )
            ))]
        summary_master[index_non_confident_blood_stromal, "classification_confidence"] <- "low_confidence"

      }



    }
    if(normal_tissue == T){
      index_potential_cancer <- which(summary_master$scATOMIC_pred %in% c("Bile Duct Cancer Cell","Bladder Cancer Cell",
                                                                          "Bone Cancer Cell","Brain Cancer Cell"    ,
                                                                          "Breast Cancer Cell","Colon/Colorectal Cancer Cell"   ,
                                                                           "Endometrial/Uterine Cancer Cell","Esophageal Cancer Cell"   ,
                                                                           "Gallbladder Cancer Cell","Gastric Cancer Cell"     ,
                                                                           "Kidney Cancer Cell","Liver Cancer Cell"   ,
                                                                           "Lung Cancer Cell","Neuroblastoma"       ,
                                                                           "Ovarian Cancer Cell","Pancreatic Cancer Cell"  ,
                                                                          "Prostate Cancer Cell","Sarcoma"     ,
                                                                          "Skin Cancer Cell","Breast/Lung/Prostate"    ,
                                                                           "Ovarian/Endometrial/Kidney","Endometrial Cancer Cell"    ,
                                                                          "Biliary/Hepatic Cancer Cell","Brain/Neuroblastoma Cancer Cell",
                                                                          "Digestive Tract Cancer Cell","Soft Tissue Cancer Cell" ,
                                                                           "Soft Tissue or Neuro Cancer Cell","Unclassified Soft Tissue or Neuro Cancer Cell",
                                                                           "ER+ Breast Cancer Cell","HER2+ Breast Cancer Cell" ,
                                                                           "Her2+ Breast Cancer Cell","TNBC Breast Cancer Cell", "Epithelial Cell","GI Tract Cell",
                                                                          "Non Blood Cell","Non Stromal Cell"))
      if(length(index_potential_cancer) > 0){
        summary_master[index_potential_cancer, "scATOMIC_pred"] <- "Normal Tissue Cell"
      }
      index_CAF <- which(summary_master$layer_3 %in% c("Cancer Associated Fibroblasts"))
      if(length(index_CAF) > 0){
        summary_master[index_CAF, "scATOMIC_pred"] <- "Normal Fibroblast"
      }
    }
    if(low_res_mode == T){
      #index CD8 T cells
      index_CD8 <- which(summary_master$layer_4 %in% c("CD8+ T cell"))
      if(length(index_CD8) > 0){
        summary_master[index_CD8, "scATOMIC_pred"] <- "CD8+ T cell"
      }
      #index CD4 T cells
      index_CD4 <- which(summary_master$layer_4 %in% c("CD4+ T cell"))
      if(length(index_CD4) > 0){
        summary_master[index_CD4, "scATOMIC_pred"] <- "CD4+ T cell"
      }
      #index macrophages
      index_macro <- which(summary_master$layer_4 %in% c("Macrophage"))
      if(length(index_macro) > 0){
        summary_master[index_macro, "scATOMIC_pred"] <- "Macrophage"
      }
      #index fibroblasts
      index_CAF <- which(summary_master$layer_3 %in% c("Cancer Associated Fibroblasts") & summary_master$scATOMIC_pred != "Normal Fibroblast")
      if(length(index_CAF) > 0){
        summary_master[index_CAF, "scATOMIC_pred"] <- "Cancer Associated Fibroblasts"
      }
      #index cDCs
      index_cDC <- which(summary_master$layer_4 %in% c("cDC"))
      if(length(index_cDC) > 0){
        summary_master[index_cDC, "scATOMIC_pred"] <- "cDC"
      }
    }
    percent_confident <- round(length(which(summary_master$classification_confidence == "confident"))/nrow(summary_master), digits = 2)
    if(percent_confident == 1){
      print(paste0("Sample classification confidence = 1.00"))
    } else{
      print(paste0("Sample classification confidence = ", percent_confident))
    }
    if(percent_confident < 0.75){
      warning("Sample classification confidence is low")
    }
    if(length(known_cancer_type) == 1){
      index_cancer_predicted <- which(summary_master$scATOMIC_pred %in% c("Bile Duct Cancer Cell", "Bladder Cancer Cell", "Bone Cancer Cell",
                                                                   "Brain Cancer Cell", "Breast Cancer Cell", "Colon/Colorectal Cancer Cell",
                                                                   "Endometrial/Uterine Cancer Cell", "Esophageal Cancer Cell", "Gallbladder Cancer Cell",
                                                                   "Gastric Cancer Cell", "Kidney Cancer Cell", "Liver Cancer Cell", "Lung Cancer Cell",
                                                                   "Neuroblastoma", "Ovarian Cancer Cell", "Pancreatic Cancer Cell", "Prostate Cancer Cell",
                                                                   "Sarcoma", "Skin Cancer Cell", "Breast/Lung/Prostate", "Ovarian/Endometrial/Kidney",
                                                                   "Endometrial Cancer Cell", "Biliary/Hepatic Cancer Cell", "Brain/Neuroblastoma Cancer Cell",
                                                                   "Digestive Tract Cancer Cell", "Soft Tissue Cancer Cell", "Soft Tissue or Neuro Cancer Cell",
                                                                   "Unclassified Soft Tissue or Neuro Cancer Cell", "ER+ Breast Cancer Cell", "HER2+ Breast Cancer Cell",
                                                                   "Her2+ Breast Cancer Cell", "TNBC Breast Cancer Cell", "Epithelial Cell", "GI Tract Cell"))
      summary_master[index_cancer_predicted, "scATOMIC_pred"] <- known_cancer_type
    }
    return(summary_master)


  }
  else{
    summary_master <- data.frame(row.names(prediction_list[["layer_1"]]),stringsAsFactors = F)
    colnames(summary_master) <- c("cell_names")
    layer_1 <- prediction_list[["layer_1"]]$predicted_tissue_with_cutoff

    summary_master <-  cbind(summary_master, layer_1)
    layer_2 <- c()
    median_score_class_layer_1 <- c()
    for (i in 1:nrow(summary_master)){
      if(summary_master[i, "layer_1"] %in% c("HSPC","B cell","CD4+ T cell", "CD8+ T cell", "macrophage_or_dendritic_cell","Mast cell",
                                             "Natural killer cell","Blood_Cell")){

        median_score_class_layer_1[i] <- median(as.numeric(prediction_list[["layer_1"]][which(prediction_list[["layer_1"]]$predicted_tissue_with_cutoff == "Blood_Cell"),"blood_score"]))
        layer_2[i] <- as.character(prediction_list[["layer_2_blood"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

      } else if(summary_master[i, "layer_1"] %in%
                c( "Bile Duct Cancer","Bladder Cancer",  "Bone Cancer",     "Brain Cancer",    "Breast Cancer",
                   "Colon/Colorectal Cancer","Endometrial/Uterine Cancer","Esophageal Cancer",
                   "Gallbladder Cancer", "Gastric Cancer", "Glial Cells",  "Kidney Cancer", "Liver Cancer", "Lung Cancer",
                   "Neuroblastoma",  "Oligodendrocytes","Ovarian Cancer",
                   "Pancreatic Cancer",
                   "Prostate Cancer", "Skin Cancer",  "Endothelial Cells", "Fibroblasts", "Myofibroblasts",  "Smooth Muscle Cells",
                   "Sarcoma", "Tissue_Cell_Normal_or_Cancer")){
        layer_2[i] <- as.character(prediction_list[["layer_2_non_blood"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_1[i] <- median(as.numeric(prediction_list[["layer_1"]][which(prediction_list[["layer_1"]]$predicted_tissue_with_cutoff == "Tissue_Cell_Normal_or_Cancer"),"cancer_normal_stromal_score"]))


      } else {
        layer_2[i] <- as.character(layer_1[i])
        median_score_class_layer_1[i] <- NA
      }
    }
    summary_master <-  cbind(summary_master, layer_2)
    median_score_class_layer_2 <- c()

    layer_3 <- c()
    for (i in 1:nrow(summary_master)){
      if(summary_master[i, "layer_2"] %in% c("CD4+ T cell", "CD8+ T cell", "Natural killer cell", "T_or_NK_lymphocyte")){
        layer_3[i] <- as.character(prediction_list[["layer_3_TNK"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_blood"]][which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff == "T_or_NK_lymphocyte"),"non_B_lymphocyte_score"]))


      } else if(summary_master[i, "layer_2"] %in% c("macrophage_or_dendritic_cell", "macrophage_DC_score")){
        layer_3[i] <- as.character(prediction_list[["layer_3_myeloid"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_blood"]][which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff == "macrophage_or_dendritic_cell"),"macrophage_DC_score"]))

      } else if(summary_master[i, "layer_2"] %in% c("B cell", "B_cell_score", "B cell or Plasmablast")){
        layer_3[i] <- as.character(prediction_list[["layer_3_BCell"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_blood"]][which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff == "B cell or Plasmablast"),"B_cell_score"]))


      } else if(summary_master[i, "layer_2"] %in% c(
        "Bile Duct Cancer","Bladder Cancer",  "Bone Cancer","Brain Cancer","Breast Cancer",
        "Colon/Colorectal Cancer","Endometrial/Uterine Cancer","Esophageal Cancer",
        "Gallbladder Cancer", "Gastric Cancer",   "Kidney Cancer", "Liver Cancer", "Lung Cancer",
        "Neuroblastoma",  "Ovarian Cancer",
        "Pancreatic Cancer",
        "Prostate Cancer", "Skin Cancer",
        "Sarcoma", "Non Stromal Cell")){
        layer_3[i] <- as.character(prediction_list[["layer_3_non_stromal"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_non_blood"]][which(prediction_list[["layer_2_non_blood"]]$predicted_tissue_with_cutoff == "Non Stromal Cell"),"non_stromal_score"]))

      } else if(summary_master[i, "layer_2"] %in%
                c("Stromal Cell","Fibroblasts", "Myofibroblasts","Smooth Muscle Cells","Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts",
                  "fibroblasts_score")){
        layer_3[i] <- as.character(prediction_list[["layer_3_stromal"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_non_blood"]][which(prediction_list[["layer_2_non_blood"]]$predicted_tissue_with_cutoff == "Stromal Cell"),"stromal_score"]))
      } else if(summary_master[i, "layer_2"] == "Cell low quality"){
        summary_master[i, "layer_2"] <- summary_master[i, "layer_1"]
        layer_3[i] <- summary_master[i, "layer_1"]
        median_score_class_layer_2[i] <- NA

      } else {
        layer_3[i] <- as.character(layer_2[i])
        if(layer_2[i] %in% c("Mast cell")){
          median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_blood"]][which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff == "Mast cell"),"mast_score"]))
        } else if(layer_2[i] %in% c("HSPC", "HSPC_score")){
          median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_blood"]][which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff %in% c("HSPC", "HSPC_score")),"HSPC_score"]))
        } else if(layer_2[i] %in% c("Glial Cells")){
          median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_non_blood"]][which(prediction_list[["layer_2_non_blood"]]$predicted_tissue_with_cutoff %in% c("Glial Cells")),"Glial Cells"]))
        } else if(layer_2[i] %in% c("Oligodendrocytes")){
          median_score_class_layer_2[i] <- median(as.numeric(prediction_list[["layer_2_non_blood"]][which(prediction_list[["layer_2_non_blood"]]$predicted_tissue_with_cutoff %in% c("Oligodendrocytes")),"Oligodendrocytes"]))
        } else{
          median_score_class_layer_2[i] <- median_score_class_layer_1[i]
        }
      }
    }
    summary_master <-  cbind(summary_master, layer_3)
    layer_4 <- c()
    median_score_class_layer_3 <- c()

    for (i in 1:nrow(summary_master)){
      if(summary_master[i, "layer_3"]  %in% c("CD4+ T cell" ,"CD4 or CD8 T cell")){
        layer_4[i] <- as.character(prediction_list[["layer_4_CD4_CD8"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_TNK"]][which(prediction_list[["layer_3_TNK"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell" ,"CD4 or CD8 T cell")),"CD4_CD8_score"]))


      } else if (summary_master[i, "layer_3"] %in% c("Natural killer cell", "NK or CD8 T cell")) {
        layer_4[i] <- as.character(prediction_list[["layer_4_CD8_NK"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_TNK"]][which(prediction_list[["layer_3_TNK"]]$predicted_tissue_with_cutoff %in% c("Natural killer cell", "NK or CD8 T cell")),"NK_CD8_score"]))

      }else if (summary_master[i, "layer_3"] == "Dendritic Cell") {
        layer_4[i] <- as.character(prediction_list[["layer_4_dendritic"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_myeloid"]][which(prediction_list[["layer_3_myeloid"]]$predicted_tissue_with_cutoff %in% c("Dendritic Cell")),"DC_score"]))

      }else if (summary_master[i, "layer_3"] == "Macrophage or Monocyte") {
        layer_4[i] <- as.character(prediction_list[["layer_4_macrophage"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_myeloid"]][which(prediction_list[["layer_3_myeloid"]]$predicted_tissue_with_cutoff %in% c("Macrophage or Monocyte")),"macrophage_score"]))

      } else if (summary_master[i, "layer_3"] %in%c(
        "Breast Cancer","Endometrial/Uterine Cancer", "Lung Cancer", "Ovarian Cancer",
        "Prostate Cancer", "Kidney Cancer",  "Non GI Epithelial Cell")) {
        layer_4[i] <- as.character(prediction_list[["layer_4_non_GI"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_non_stromal"]][which(prediction_list[["layer_3_non_stromal"]]$predicted_tissue_with_cutoff %in% c("Breast Cancer","Endometrial/Uterine Cancer", "Lung Cancer", "Ovarian Cancer",
                                                                                                                                                                                       "Prostate Cancer", "Kidney Cancer",  "Non GI Epithelial Cell")),"non_GI_score"]))

      } else if (summary_master[i, "layer_3"] %in% c("Bile Duct Cancer","Bladder Cancer",
                                                     "Colon/Colorectal Cancer","Esophageal Cancer",
                                                     "Gallbladder Cancer", "Gastric Cancer", "Liver Cancer",
                                                     "Pancreatic Cancer","GI Epithelial Cell")) {
        layer_4[i] <- as.character(prediction_list[["layer_4_GI"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_non_stromal"]][which(prediction_list[["layer_3_non_stromal"]]$predicted_tissue_with_cutoff %in% c("Bile Duct Cancer","Bladder Cancer",
                                                                                                                                                                                       "Colon/Colorectal Cancer","Esophageal Cancer",
                                                                                                                                                                                       "Gallbladder Cancer", "Gastric Cancer", "Liver Cancer",
                                                                                                                                                                                       "Pancreatic Cancer","GI Epithelial Cell")),"GI_score"]))

      } else if (summary_master[i, "layer_3"] %in% c("Bone Cancer","Brain Cancer",
                                                     "Neuroblastoma", "Skin Cancer",
                                                     "Sarcoma", "Soft Tissue or Neuro Cancer Cell")) {
        layer_4[i] <- as.character(prediction_list[["layer_4_soft_tissue_neuro"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_non_stromal"]][which(prediction_list[["layer_3_non_stromal"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer","Brain Cancer",
                                                                                                                                                                                       "Neuroblastoma", "Skin Cancer",
                                                                                                                                                                                       "Sarcoma", "Soft Tissue or Neuro Cancer Cell")),"soft_tissue_neuro_score"]))



      } else if(summary_master[i, "layer_3"] %in% c("Cancer Associated Fibroblasts")){
        layer_4[i] <- as.character(prediction_list[["layer_4_CAF"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

        median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_stromal"]][which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff %in% c("Cancer Associated Fibroblasts")),"Cancer Associated Fibroblasts"]))

      } else if(summary_master[i, "layer_3"] == "Cell low quality"){
        summary_master[i, "layer_3"] <- summary_master[i, "layer_2"]
        layer_4[i] <- summary_master[i, "layer_2"]
        median_score_class_layer_3[i] <- NA

      } else {
        layer_4[i] <- as.character(layer_3[i])
        if(layer_3[i] %in% c("B Cell")){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_BCell"]][which(prediction_list[["layer_3_BCell"]]$predicted_tissue_with_cutoff == "B Cell"),"B_cell_score"]))
        } else if(layer_3[i] %in% c("Plasmablast")){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_BCell"]][which(prediction_list[["layer_3_BCell"]]$predicted_tissue_with_cutoff == "Plasmablast"),"plasmablast_score"]))
        } else if(layer_3[i] %in% c("Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts" )){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_stromal"]][which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff %in% c("Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts")),"cancer_associated_fibroblast_score"]))
        } else if(layer_3[i] %in% c("Smooth Muscle Cells" )){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_stromal"]][which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff %in% c("Smooth Muscle Cells")),"Smooth Muscle Cells"]))
        } else if(layer_3[i] %in% c("Fibroblasts", "Myofibroblasts")){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_stromal"]][which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff %in% c("Fibroblasts", "Myofibroblasts")),"fibroblasts_score"]))
        } else if(layer_3[i] %in% c("Endothelial Cells")){
          median_score_class_layer_3[i] <- median(as.numeric(prediction_list[["layer_3_stromal"]][which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff %in% c("Endothelial Cells")),"Endothelial Cells"]))
        } else{
          median_score_class_layer_3[i] <- median_score_class_layer_2[i]
        }


      }
    }
    summary_master <-  cbind(summary_master, layer_4)
    median_score_class_layer_4 <- c()
    layer_5 <- c()
    for (i in 1:nrow(summary_master)){
      if(summary_master[i, "layer_4"] %in%c(
        "Endometrial/Uterine Cancer", "Ovarian Cancer", "Kidney Cancer", "Ovarian/Endometrial/Kidney Cell")){
        layer_5[i] <- as.character(prediction_list[["layer_5_ov_endo_kid"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_non_GI"]][which(prediction_list[["layer_4_non_GI"]]$predicted_tissue_with_cutoff %in% c("Endometrial/Uterine Cancer", "Ovarian Cancer", "Kidney Cancer", "Ovarian/Endometrial/Kidney Cell")),"ov_endo_kid_score"]))

      } else if(summary_master[i, "layer_4"] %in%c(
        "Lung Cancer") & summary_master[i, "layer_3"] %in%c(
          "Soft Tissue or Neuro Cancer Cell")){
        layer_5[i] <- "Lung Cancer"
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_soft_tissue_neuro"]][which(prediction_list[["layer_4_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Lung Cancer")),"Lung Cancer"]))

      } else if(summary_master[i, "layer_4"] %in%c(
        "Breast Cancer","Lung Cancer",
        "Prostate Cancer",  "Breast/Lung/Prostate Cell")){
        layer_5[i] <- as.character(prediction_list[["layer_5_breast_lung_prostate"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_non_GI"]][which(prediction_list[["layer_4_non_GI"]]$predicted_tissue_with_cutoff %in% c("Breast Cancer","Lung Cancer",
                                                                                                                                                                             "Prostate Cancer",  "Breast/Lung/Prostate Cell")),"breast_lung_prostate_score"]))

      } else if(summary_master[i, "layer_4"] %in% c("Bile Duct Cancer","Bladder Cancer",
                                                    "Gallbladder Cancer",  "Liver Cancer",
                                                    "Pancreatic Cancer", "Billiary Cell")){
        layer_5[i] <- as.character(prediction_list[["layer_5_biliary"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_GI"]][which(prediction_list[["layer_4_GI"]]$predicted_tissue_with_cutoff %in% c("Bile Duct Cancer","Bladder Cancer",
                                                                                                                                                                     "Gallbladder Cancer",  "Liver Cancer",
                                                                                                                                                                     "Pancreatic Cancer", "Billiary Cell")),"billiary_score"]))


      } else if(summary_master[i, "layer_4"] %in% c("Colon/Colorectal Cancer","Esophageal Cancer",
                                                    "Gastric Cancer","Colorectal/Esophageal/Gastric Cell")){
        layer_5[i] <- as.character(prediction_list[["layer_5_digestive"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_GI"]][which(prediction_list[["layer_4_GI"]]$predicted_tissue_with_cutoff %in% c("Colon/Colorectal Cancer","Esophageal Cancer",
                                                                                                                                                                     "Gastric Cancer","Colorectal/Esophageal/Gastric Cell")),"digestive_score"]))
      } else if(summary_master[i, "layer_4"] %in% c("Bone Cancer","Brain Cancer",
                                                    "Neuroblastoma",
                                                    "Sarcoma", "Soft Tissue or Neuro Cancer Cell", "Soft Tissue or Neuro Cancer Cell")){
        layer_5[i] <- as.character(prediction_list[["layer_5_soft_tissue_neuro"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_soft_tissue_neuro"]][which(prediction_list[["layer_4_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer","Brain Cancer",
                                                                                                                                                                                                   "Neuroblastoma",
                                                                                                                                                                                                   "Sarcoma", "Soft Tissue or Neuro Cancer Cell", "Soft Tissue or Neuro Cancer Cell")),"soft_tissue_neuro_score"]))

      }  else if(summary_master[i, "layer_4"] %in% c("cDC")){
        layer_5[i] <- as.character(prediction_list[["layer_5_cDC"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_dendritic"]][which(prediction_list[["layer_4_dendritic"]]$predicted_tissue_with_cutoff %in% c("cDC")),"cDC"]))

      } else if(summary_master[i, "layer_4"] %in% c("Macrophage")){
        layer_5[i] <- as.character(prediction_list[["layer_5_macrophage"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_macrophage"]][which(prediction_list[["layer_4_macrophage"]]$predicted_tissue_with_cutoff %in% c("Macrophage")),"Macrophage"]))

      } else if(summary_master[i, "layer_4"] %in% c("Monocyte")){
        layer_5[i] <- as.character(prediction_list[["layer_5_monocyte"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

        median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_macrophage"]][which(prediction_list[["layer_4_macrophage"]]$predicted_tissue_with_cutoff %in% c("Monocyte")),"Monocyte"]))

      } else if(summary_master[i, "layer_4"] %in% c("CD4+ T cell")){
        if(fine_grained_T == T){
          layer_5[i] <- as.character(prediction_list[["layer_5_CD4"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell")),"CD4+ T cell"]))

        } else{
          layer_5[i] <- as.character(layer_4[i])
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell")),"CD4+ T cell"]))
        }
      } else if(summary_master[i, "layer_4"] %in% c("CD8+ T cell")){
        if(fine_grained_T == T){
          layer_5[i] <- as.character(prediction_list[["layer_5_CD8"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

          if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "CD4 or CD8 T cell"){
            median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
          } else if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "NK or CD8 T cell"){
            median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD8_NK"]][which(prediction_list[["layer_4_CD8_NK"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
          }
        } else{
          layer_5[i] <- as.character(layer_4[i])
          if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "CD4 or CD8 T cell"){
            median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
          } else if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "NK or CD8 T cell"){
            median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD8_NK"]][which(prediction_list[["layer_4_CD8_NK"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
          }
        }
      } else if(summary_master[i, "layer_4"] == "Cell low quality"){
        summary_master[i, "layer_4"] <- summary_master[i, "layer_3"]
        median_score_class_layer_4[i] <- NA
        layer_5[i] <- summary_master[i, "layer_3"]
      } else {
        layer_5[i] <- as.character(layer_4[i])
        if(layer_4[i] %in% c("Skin Cancer")){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_soft_tissue_neuro"]][which(prediction_list[["layer_4_soft_tissue_neuro"]]$predicted_tissue_with_cutoff == "Skin Cancer"),"Skin Cancer"]))
        } else if(layer_4[i] %in% c("cDC")){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_dendritic"]][which(prediction_list[["layer_4_dendritic"]]$predicted_tissue_with_cutoff == "cDC"),"cDC"]))
        } else if(layer_4[i] %in% c("pDC")){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_dendritic"]][which(prediction_list[["layer_4_dendritic"]]$predicted_tissue_with_cutoff == "pDC"),"pDC"]))
        } else if(layer_4[i] %in% c("ASDC")){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_dendritic"]][which(prediction_list[["layer_4_dendritic"]]$predicted_tissue_with_cutoff == "ASDC"),"ASDC"]))
        } else if(layer_4[i] %in% c("Macrophage" )){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_macrophage"]][which(prediction_list[["layer_4_macrophage"]]$predicted_tissue_with_cutoff == "Macrophage"),"Macrophage"]))
        } else if(layer_4[i] %in% c("Monocyte" )){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_macrophage"]][which(prediction_list[["layer_4_macrophage"]]$predicted_tissue_with_cutoff == "Monocyte"),"Monocyte"]))
        } else if(layer_4[i] %in% c("CD4+ T cell") ){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell")),"CD4+ T cell"]))
        } else if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "CD4 or CD8 T cell"){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD4_CD8"]][which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
        } else if(layer_4[i] %in% c("CD8+ T cell") & layer_3[i] == "NK or CD8 T cell"){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD8_NK"]][which(prediction_list[["layer_4_CD8_NK"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell")),"CD8+ T cell"]))
        } else if(layer_4[i] %in% c("Natural killer cell")){
          median_score_class_layer_4[i] <- median(as.numeric(prediction_list[["layer_4_CD8_NK"]][which(prediction_list[["layer_4_CD8_NK"]]$predicted_tissue_with_cutoff %in% c("Natural killer cell")),"Natural killer cell"]))
        } else{
          median_score_class_layer_4[i] <- median_score_class_layer_3[i]
        }
      }
    }


    summary_master <-  cbind(summary_master, layer_5)

    layer_6 <- c()
    median_score_class_layer_5 <- c()
    for (i in 1:nrow(summary_master)){
      if(summary_master[i, "layer_5"] %in% c("Brain Cancer",
                                             "Neuroblastoma", "Brain/Neuroblastoma Cancer Cell")){
        layer_6[i] <- as.character(prediction_list[["layer_6_brain_nbm"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_soft_tissue_neuro"]][which(prediction_list[["layer_5_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Brain Cancer",
                                                                                                                                                                                                   "Neuroblastoma", "Brain/Neuroblastoma Cancer Cell")),"brain_nbm_score"]))

      } else if(summary_master[i, "layer_5"] %in% c("Bone Cancer",
                                                    "Sarcoma", "Soft Tissue Cancer Cell")){
        layer_6[i] <- as.character(prediction_list[["layer_6_soft_tissue"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])
        median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_soft_tissue_neuro"]][which(prediction_list[["layer_5_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer",
                                                                                                                                                                                                   "Sarcoma", "Soft Tissue Cancer Cell")),"soft_tissue_score"]))


      } else if(summary_master[i, "layer_5"] %in% c("Breast Cancer Cell")){
        if(breast_mode == T){
          layer_6[i] <- as.character(prediction_list[["layer_6_breast"]][summary_master$cell_names[i],"predicted_tissue_with_cutoff"])

          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_breast_lung_prostate"]][which(prediction_list[["layer_5_breast_lung_prostate"]]$predicted_tissue_with_cutoff %in% c("Breast Cancer Cell")),"Breast Cancer"]))

        } else{
          layer_6[i] <- as.character(layer_5[i])
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_breast_lung_prostate"]][which(prediction_list[["layer_5_breast_lung_prostate"]]$predicted_tissue_with_cutoff %in% c("Breast Cancer Cell")),"Breast Cancer"]))
        }
      } else if(summary_master[i, "layer_5"] == "Cell low quality"){
        summary_master[i, "layer_5"] <- summary_master[i, "layer_4"]
        layer_6[i] <- summary_master[i, "layer_4"]
        median_score_class_layer_5[i] <- NA
      } else {
        layer_6[i] <- as.character(layer_5[i])
        if(layer_5[i] == "Lung Cancer Cell" & layer_4[i] == "Breast/Lung/Prostate Cell"){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_breast_lung_prostate"]][which(prediction_list[["layer_5_breast_lung_prostate"]]$predicted_tissue_with_cutoff %in% c("Lung Cancer Cell")),"Lung Cancer"]))
        } else if(layer_5[i] == "Prostate Cancer Cell" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_breast_lung_prostate"]][which(prediction_list[["layer_5_breast_lung_prostate"]]$predicted_tissue_with_cutoff %in% c("Prostate Cancer Cell")),"Prostate Cancer"]))
        } else if(layer_5[i] == "Ovarian Cancer Cell" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_ov_endo_kid"]][which(prediction_list[["layer_5_ov_endo_kid"]]$predicted_tissue_with_cutoff %in% c("Ovarian Cancer Cell")),"Ovarian Cancer"]))
        } else if(layer_5[i] == "Kidney Cancer Cell" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_ov_endo_kid"]][which(prediction_list[["layer_5_ov_endo_kid"]]$predicted_tissue_with_cutoff %in% c("Kidney Cancer Cell")),"Kidney Cancer"]))
        } else if(layer_5[i] == "Endometrial Cancer Cell" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_ov_endo_kid"]][which(prediction_list[["layer_5_ov_endo_kid"]]$predicted_tissue_with_cutoff %in% c("Endometrial Cancer Cell")),"Endometrial/Uterine Cancer"]))
        } else if(layer_5[i] == "Pancreatic Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_biliary"]][which(prediction_list[["layer_5_biliary"]]$predicted_tissue_with_cutoff %in% c("Pancreatic Cancer")),"Pancreatic Cancer"]))
        } else if(layer_5[i] == "Liver Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_biliary"]][which(prediction_list[["layer_5_biliary"]]$predicted_tissue_with_cutoff %in% c("Liver Cancer")),"Liver Cancer"]))
        } else if(layer_5[i] == "Gallbladder Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_biliary"]][which(prediction_list[["layer_5_biliary"]]$predicted_tissue_with_cutoff %in% c("Gallbladder Cancer")),"Gallbladder Cancer"]))
        } else if(layer_5[i] == "Bladder Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_biliary"]][which(prediction_list[["layer_5_biliary"]]$predicted_tissue_with_cutoff %in% c("Bladder Cancer")),"Bladder Cancer"]))
        } else if(layer_5[i] == "Bile Duct Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_biliary"]][which(prediction_list[["layer_5_biliary"]]$predicted_tissue_with_cutoff %in% c("Bile Duct Cancer")),"Bile Duct Cancer"]))
        } else if(layer_5[i] == "Colon/Colorectal Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_digestive"]][which(prediction_list[["layer_5_digestive"]]$predicted_tissue_with_cutoff %in% c("Colon/Colorectal Cancer")),"Colon/Colorectal Cancer"]))
        } else if(layer_5[i] == "Esophageal Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_digestive"]][which(prediction_list[["layer_5_digestive"]]$predicted_tissue_with_cutoff %in% c("Esophageal Cancer")),"Esophageal Cancer"]))
        } else if(layer_5[i] == "Gastric Cancer" ){
          median_score_class_layer_5[i] <- median(as.numeric(prediction_list[["layer_5_digestive"]][which(prediction_list[["layer_5_digestive"]]$predicted_tissue_with_cutoff %in% c("Gastric Cancer")),"Gastric Cancer"]))
        } else{
          median_score_class_layer_5[i] <- median_score_class_layer_4[i]
        }
      }
    }

    summary_master <-  cbind(summary_master, layer_6)
    median_score_class_layer_6 <- c()
    for (i in 1:nrow(summary_master)){
      if(layer_6[i] == "TNBC Breast Cancer Cell"){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_breast"]][which(prediction_list[["layer_6_breast"]]$predicted_tissue_with_cutoff %in% c("TNBC Breast Cancer Cell")),"TNBC"]))
      } else if(layer_6[i] == "Her2+ Breast Cancer Cell"){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_breast"]][which(prediction_list[["layer_6_breast"]]$predicted_tissue_with_cutoff %in% c("Her2+ Breast Cancer Cell")),"HER2+"]))
      } else if(layer_6[i] == "ER+ Breast Cancer Cell"){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_breast"]][which(prediction_list[["layer_6_breast"]]$predicted_tissue_with_cutoff %in% c("ER+ Breast Cancer Cell")),"ER+"]))
      } else if(layer_6[i] == "Brain Cancer" ){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_brain_nbm"]][which(prediction_list[["layer_6_brain_nbm"]]$predicted_tissue_with_cutoff %in% c("Brain Cancer")),"Brain Cancer"]))
      } else if(layer_6[i] == "Neuroblastoma" ){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_brain_nbm"]][which(prediction_list[["layer_6_brain_nbm"]]$predicted_tissue_with_cutoff %in% c("Neuroblastoma")),"Neuroblastoma"]))
      } else if(layer_6[i] == "Bone Cancer" ){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_soft_tissue"]][which(prediction_list[["layer_6_soft_tissue"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer")),"Bone Cancer"]))
      } else if(layer_6[i] == "Sarcoma" ){
        median_score_class_layer_6[i] <- median(as.numeric(prediction_list[["layer_6_soft_tissue"]][which(prediction_list[["layer_6_soft_tissue"]]$predicted_tissue_with_cutoff %in% c("Sarcoma")),"Sarcoma"]))
      } else{
        median_score_class_layer_6[i] <- median_score_class_layer_5[i]
      }
    }
    summary_master <- cbind(summary_master, median_score_class_layer_1,median_score_class_layer_2,median_score_class_layer_3,median_score_class_layer_4,median_score_class_layer_5,median_score_class_layer_6  )

    row.names(summary_master) <- summary_master$cell_names
    summary_master$layer_4 <- gsub("Cancer$", "Cancer Cell", summary_master$layer_4)
    summary_master$layer_5 <- gsub("Cancer$", "Cancer Cell", summary_master$layer_5)
    summary_master$layer_6 <- gsub("Cancer$", "Cancer Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_B_cell_or_plasmablast", "B cell or Plasmablast", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_blood_cell", "Blood Cell", summary_master$layer_6)

    summary_master$layer_6 <- gsub("unclassified_any_cell", "Any Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_normal_or_cancer_tissue", "Non Blood Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_T_or_NK_cell", "T or NK Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_macrophage_or_DC", "Macrophage or Dendritic Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_non_GI_epithelial_cell", "Epithelial Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("unclassified_GI_epithelial_cell", "GI Tract Cell", summary_master$layer_6)
    summary_master$layer_6 <- gsub("^NF$", "Normal Fibroblast", summary_master$layer_6)
    if(modify_results == TRUE){
      scATOMIC_pred <- summary_master$layer_6
      summary_master <- cbind(summary_master, scATOMIC_pred)
      if(use_CNVs == TRUE){
        seurat_object <- Seurat::CreateSeuratObject(raw_counts)
        seurat_object <- Seurat::NormalizeData(seurat_object)
        seurat_object <- Seurat::FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
        all.genes <- rownames(seurat_object)
        seurat_object <- Seurat::ScaleData(seurat_object, features = all.genes)
        seurat_object <- Seurat::RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
        seurat_object <- Seurat::FindNeighbors(seurat_object, dims = 1:50)
        seurat_object <- Seurat::FindClusters(seurat_object, resolution = 0.5)
        summary_master <- cbind(summary_master, seurat_object@meta.data[row.names(summary_master), "seurat_clusters"])
        colnames(summary_master)[ncol(summary_master)] <- "seurat_clusters"

        copy_kat_res <- scATOMIC::copy_kat_no_heatmap(rawmat = raw_counts, summary_matrix = summary_master,
                                                                id.type = "S", cell.line = "no",
                                                                ngene.chr = 5, LOW.DR = 0.05, UP.DR = 0.1,
                                                                win.size = 25, KS.cut = 0.1, sam.name = "", distance = "euclidean",
                                                                n.cores = mc.cores)
        copy_kat_res <- as.data.frame(copy_kat_res[,"copykat.pred"], row.names = row.names(copy_kat_res))
        colnames(copy_kat_res) <- "CNV_status"
        summary_master <- merge(summary_master,copy_kat_res,by="row.names",all.x=TRUE)
        row.names(summary_master) <- summary_master$Row.names
        summary_master <- summary_master[,-1]
      }
      s.genes <- Seurat::cc.genes$s.genes
      g2m.genes <- Seurat::cc.genes$g2m.genes
      seurat_object <- CreateSeuratObject(raw_counts, meta.data = summary_master)
      seurat_object <- Seurat::NormalizeData(seurat_object, verbose = F)
      seurat_object <- Seurat::FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = F)
      all.genes <- rownames(seurat_object)
      #try regressing cell cycle genes
      seurat_object <- Seurat::CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
      seurat_object <- Seurat::ScaleData(seurat_object, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(object = seurat_object))
      if(ncol(seurat_object) < 50){
        seurat_object <- Seurat::RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), verbose = F, npcs = (ncol(seurat_object) - 1))

      } else{
        seurat_object <- Seurat::RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), verbose = F)
      }
      seurat_object <- Seurat::FindNeighbors(seurat_object, dims = 1:10, verbose = F)
      seurat_object <- Seurat::FindClusters(seurat_object, resolution = 0.2, verbose = F)

      final_classes <- model_layer_1$classes
      final_classes <- gsub("Cancer$", "Cancer Cell", final_classes)
      cancer_classes <- final_classes[grep("Cancer Cell|oma$|Breast/Lung/Prostate|Ovarian/Endometrial/Kidney|Colorectal Cancer Cell|Endometrial Cancer Cell|Biliary/Hepatic Cancer Cell", final_classes)]
      cancer_classes <- c(cancer_classes,"Breast/Lung/Prostate", "Ovarian/Endometrial/Kidney", "Endometrial Cancer Cell",
                          "Biliary/Hepatic Cancer Cell", "Brain/Neuroblastoma Cancer Cell", "Digestive Tract Cancer Cell",
                          "Soft Tissue Cancer Cell", "Soft Tissue or Neuro Cancer Cell", "Unclassified Soft Tissue or Neuro Cancer Cell","ER+ Breast Cancer Cell",
                          "HER2+ Breast Cancer Cell","Her2+ Breast Cancer Cell","TNBC Breast Cancer Cell")
      predicted_cancer <- seurat_object@meta.data[which(seurat_object@meta.data$layer_6 %in% cancer_classes),]
      if(nrow(predicted_cancer) > 0){
        frequency <- as.data.frame(table(predicted_cancer$scATOMIC_pred), stringsAsFactors = F)
        if(breast_mode == T){
          frequency <- as.data.frame(table(predicted_cancer$layer_5), stringsAsFactors = F)

        }
        proportion <- frequency$Freq/sum(frequency$Freq)
        frequency <- cbind(frequency, proportion)
        max_proportion <- max(frequency$proportion)
        major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
        if(length(major_cancer) > 1){
          major_cancer <- major_cancer[1]
        }

        #if there is more than 1 seurat cluster for cancer we want to split into normal and cancer
        if(major_cancer == "Bile Duct Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "CHOL") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "CHOL")),"Gene"]
        }else if(major_cancer == "Bladder Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "BLCA") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "BLCA")),"Gene"]
        } else if(major_cancer == "Bone Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "SARC") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "SARC")),"Gene"]
        } else if(major_cancer == "Brain Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "GBM", "LGG") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "GBM", "LGG")),"Gene"]
        } else if(major_cancer %in% c( "Breast Cancer Cell","ER+ Breast Cancer Cell","HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell","Her2+ Breast Cancer Cell"  ) ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "BRCA") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "BRCA")),"Gene"]
        } else if(major_cancer == "Colon/Colorectal Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "COAD") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "COAD")),"Gene"]
        } else if(major_cancer %in% c("Endometrial/Uterine Cancer Cell", "Endometrial Cancer Cell" )){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "UCS", "UCEC") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "UCS", "UCEC")),"Gene"]
        } else if(major_cancer == "Gallbladder Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "CHOL") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "CHOL")),"Gene"]
        } else if(major_cancer == "Kidney Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "KICH", "KIRC", "KIRP") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "KICH", "KIRC", "KIRP")),"Gene"]
        } else if(major_cancer == "Liver Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "LIHC") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "LIHC")),"Gene"]
        } else if(major_cancer == "Lung Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "LUAD") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "LUAD")),"Gene"]
        } else if(major_cancer == "Ovarian Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "OV") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "OV")),"Gene"]
        } else if(major_cancer == "Pancreatic Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "PAAD") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "PAAD")),"Gene"]
        } else if(major_cancer == "Prostate Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "PRAD") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "PRAD")),"Gene"]
        } else if(major_cancer == "Sarcoma" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "SARC") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "SARC")),"Gene"]
        } else if(major_cancer == "Skin Cancer Cell" ){
          cancer_specific_upreg <- upreg_list[which(upreg_list$cancer_type %in% c( "SKCM") ),"Gene"]
          cancer_specific_downreg <- downreg_list[which(downreg_list$cancer_type %in% c( "SKCM")),"Gene"]
        } else{
          cancer_specific_upreg <- pan_cancer_upreg
          cancer_specific_downreg <- pan_cancer_downreg
        }
        if(pan_cancer == T){
          cancer_specific_upreg <- pan_cancer_upreg
          cancer_specific_downreg <- pan_cancer_downreg
        }
        if(max_proportion < min_prop){
          cancer_specific_upreg <- pan_cancer_upreg
          cancer_specific_downreg <- pan_cancer_downreg
        }
        genes_upreg_shared <- intersect(row.names(raw_counts), cancer_specific_upreg)
        genes_downreg_shared <- intersect(row.names(raw_counts), cancer_specific_downreg)
        index_cancer <- row.names(predicted_cancer)
        count_matrix_up_reg <- raw_counts[c(genes_upreg_shared),index_cancer]
        if(length(index_cancer) == 1){
          count_matrix_up_reg <- matrix(count_matrix_up_reg, nrow = length(count_matrix_up_reg), ncol = length(index_cancer), dimnames = list(names(count_matrix_up_reg), index_cancer))
        }
        frequency_per_row_nozero <- apply(count_matrix_up_reg,1,function(x){length(which(x!=0))})/ncol(count_matrix_up_reg)
        genes_upreg_shared <- names(frequency_per_row_nozero)[which(frequency_per_row_nozero > 0)]
        count_matrix_down_reg <- raw_counts[c(genes_downreg_shared),index_cancer]
        if(length(index_cancer) == 1){
          count_matrix_down_reg <- matrix(count_matrix_down_reg, nrow = length(count_matrix_down_reg), ncol = length(index_cancer), dimnames = list(names(count_matrix_down_reg), index_cancer))
        }
        frequency_per_row_nozero <- apply(count_matrix_down_reg,1,function(x){length(which(x!=0))})/ncol(count_matrix_down_reg)
        genes_downreg_shared <- names(frequency_per_row_nozero)[which(frequency_per_row_nozero > 0)]
        cancer_subset <- seurat_object
        cancer_subset <- magic(cancer_subset, seed = 123)
        #cancer_subset <- Seurat::AddModuleScore(cancer_subset, features = genes_upreg_shared, verbose = F, name = "upreg_genes", assay = "MAGIC_RNA")
        #cancer_subset <- Seurat::AddModuleScore(cancer_subset, features = genes_downreg_shared, verbose = F, name = "downreg_genes", assay = "MAGIC_RNA")
        cancer_subset <- Seurat::AddModuleScore(cancer_subset, features = list(upreg_genes = genes_upreg_shared, downreg_genes = genes_downreg_shared), verbose = F, name = c("upreg_genes", "downreg_genes"), assay = "MAGIC_RNA")

        dist_cancer_signature <- amap::Dist(cancer_subset@meta.data[,c("upreg_genes1", "downreg_genes2")], method = "euclidean", nbproc = mc.cores) # distance matrix
        fit_clusters <- hclust(dist_cancer_signature, method="ward.D2")

        groups <- cutree(fit_clusters, k=2)
        cancer_subset <- Seurat::AddMetaData(cancer_subset, groups, "pan_cancer_cluster")
        index_normal <- unique(c(row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$upreg_genes1 <
                                                                            cancer_subset@meta.data$downreg_genes2)],
                                 row.names(cancer_subset@meta.data)[-which(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes)]))
        index_non_normal <- outersect(row.names(cancer_subset@meta.data), index_normal)
        pct_normal_cluster_1 <- length(intersect(index_normal, row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == 1)]))/length(which(cancer_subset@meta.data$pan_cancer_cluster == 1))

        pct_normal_cluster_2 <- length(intersect(index_normal, row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == 2)]))/length(which(cancer_subset@meta.data$pan_cancer_cluster == 2))
        # DimPlot(cancer_subset, group.by = "pan_cancer_cluster")
        # DimPlot(cancer_subset, group.by = "scATOMIC_pred", label = T) + NoLegend()
        # DimPlot(cancer_subset, cells.highlight =  index_normal)

        if (max(c(pct_normal_cluster_1,pct_normal_cluster_2)) - min(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.1){
          clust_1_cancer <- length(which(cancer_subset@meta.data$pan_cancer_cluster == 1 &
                                           cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes ))/nrow(cancer_subset@meta.data)
          clust_2_cancer <- length(which(cancer_subset@meta.data$pan_cancer_cluster == 2 &
                                           cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes ))/nrow(cancer_subset@meta.data)

          if(clust_1_cancer > clust_2_cancer){
            cancer_clust = "1"

          } else{
            cancer_clust = "2"

          }
          index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes )]

        } else if(pct_normal_cluster_1 > pct_normal_cluster_2){
          cancer_clust = "2"
          if(max(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.5){

            index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal) )]

          } else{
            index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust)]

          }


        } else if (pct_normal_cluster_1 < pct_normal_cluster_2){
          cancer_clust = "1"
          if(max(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.5){
            index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal) )]

          } else{
            index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust)]

          }

        }
        continue = TRUE
        if(length(index_cancer) < 20){
          continue = FALSE
        }
        if(pct_normal_cluster_1 != pct_normal_cluster_2){
          if(min(c(pct_normal_cluster_1,pct_normal_cluster_2)) > 0 & max(pct_normal_cluster_1,pct_normal_cluster_2) > 0.5){
            while(continue == TRUE){
              index_cancer_last <- index_cancer
              cancer_subset <- subset(seurat_object, cells = index_cancer)
              cancer_subset <- Seurat::NormalizeData(cancer_subset, verbose = F)
              cancer_subset <- Seurat::ScaleData(cancer_subset, verbose = F)
              cancer_subset <- magic(cancer_subset, seed = 123)
              cancer_subset <- Seurat::AddModuleScore(cancer_subset, features = list(upreg_genes = genes_upreg_shared, downreg_genes = genes_downreg_shared), verbose = F, name = c("upreg_genes", "downreg_genes"), assay = "MAGIC_RNA")

              dist_cancer_signature <- amap::Dist(cancer_subset@meta.data[,c("upreg_genes1", "downreg_genes2")], method = "euclidean", nbproc = mc.cores) # distance matrix


              #hierarcherical clustering based on these scores
              fit_clusters <- hclust(dist_cancer_signature, method="ward.D2")

              groups <- cutree(fit_clusters, k=2)
              cancer_subset <- Seurat::AddMetaData(cancer_subset, groups, "pan_cancer_cluster")
              # DimPlot(cancer_subset,group.by = "pan_cancer_cluster" )
              index_normal <- unique(c(row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$upreg_genes1 <
                                                                                  cancer_subset@meta.data$downreg_genes2)],
                                       row.names(cancer_subset@meta.data)[-which(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes)]))
              index_non_normal <- outersect(row.names(cancer_subset@meta.data), index_normal)
              #DimPlot(cancer_subset, cells.highlight =  index_normal)

              pct_normal_cluster_1 <- length(intersect(index_normal, row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == 1)]))/length(which(cancer_subset@meta.data$pan_cancer_cluster == 1))

              pct_normal_cluster_2 <- length(intersect(index_normal, row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == 2)]))/length(which(cancer_subset@meta.data$pan_cancer_cluster == 2))
              if(max(pct_normal_cluster_1,pct_normal_cluster_2) > 0.2){
                continue <- TRUE
              } else{
                continue <- FALSE
              }
              if (max(c(pct_normal_cluster_1,pct_normal_cluster_2)) - min(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.25){
                clust_1_cancer <- length(which(cancer_subset@meta.data$pan_cancer_cluster == 1 &
                                                 cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes ))/nrow(cancer_subset@meta.data)
                clust_2_cancer <- length(which(cancer_subset@meta.data$pan_cancer_cluster == 2 &
                                                 cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes ))/nrow(cancer_subset@meta.data)

                if(clust_1_cancer > clust_2_cancer){
                  cancer_clust = "1"

                } else{
                  cancer_clust = "2"

                }
                index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal)]


              } else if (max(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.9 ){
                if(pct_normal_cluster_1 > pct_normal_cluster_2){
                  cancer_clust = "2"
                } else{
                  cancer_clust = "1"

                }
                index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal)]

              }else if(pct_normal_cluster_1 > pct_normal_cluster_2){
                cancer_clust = "2"
                if(max(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.5 ){
                  index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal) )]
                  continue <- TRUE
                } else{
                  index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust)]

                }


              } else if (pct_normal_cluster_1 < pct_normal_cluster_2){
                cancer_clust = "1"
                if( max(c(pct_normal_cluster_1,pct_normal_cluster_2)) < 0.5 ){
                  index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust |(cancer_subset@meta.data$scATOMIC_pred %in% cancer_classes & row.names(cancer_subset@meta.data) %in% index_non_normal) )]
                  continue <- TRUE
                } else{
                  index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$pan_cancer_cluster == cancer_clust)]

                }

              }

              if(length(index_cancer_last) == length(index_cancer)){
                if(length(levels(as.factor(index_cancer_last == index_cancer)) == TRUE) == 1){
                  if (levels(as.factor(index_cancer_last == index_cancer)) == TRUE){
                    continue = FALSE
                  }
                }
              }

              if(continue == FALSE){
                index_cancer <- row.names(cancer_subset@meta.data)[which(cancer_subset@meta.data$upreg_genes1 >  cancer_subset@meta.data$downreg_genes2 )]
              }
              if(length(index_cancer) < 20){
                continue = FALSE
              }
            }
            index_normal <- outersect(row.names(seurat_object@meta.data)[which(seurat_object@meta.data$scATOMIC_pred %in% cancer_classes)],
                                      index_cancer)
            seurat_object@meta.data[index_normal,"scATOMIC_pred"] <- "Normal Tissue Cell"
            predicted_cancer_new <- seurat_object@meta.data[which(seurat_object@meta.data$scATOMIC_pred %in%c(cancer_classes,"Breast/Lung/Prostate", "Ovarian/Endometrial/Kidney", "Endometrial Cancer Cell",
                                                                                                              "Biliary/Hepatic Cancer Cell", "Brain/Neuroblastoma Cancer Cell", "Digestive Tract Cancer Cell",
                                                                                                              "Soft Tissue Cancer Cell", "Soft Tissue or Neuro Cancer Cell", "Unclassified Soft Tissue or Neuro Cancer Cell","ER+ Breast Cancer Cell",
                                                                                                              "HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell", "Her2+ Breast Cancer Cell")),]

            if(breast_mode == F){
              frequency <- as.data.frame(table(predicted_cancer_new$scATOMIC_pred), stringsAsFactors = F)
              proportion <- frequency$Freq/sum(frequency$Freq)
              frequency <- cbind(frequency, proportion)
              max_proportion <- max(frequency$proportion)
              major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
              if(length(major_cancer) > 1){
                major_cancer <- major_cancer[1]
              }
              seurat_object@meta.data[row.names(predicted_cancer_new),"scATOMIC_pred"] <- major_cancer
              seurat_object@meta.data[,"pan_cancer_cluster"] <- "Normal"
              seurat_object@meta.data[row.names(predicted_cancer_new),"pan_cancer_cluster"] <- "Cancer"
            } else{
              frequency <- as.data.frame(table(predicted_cancer_new$layer_5), stringsAsFactors = F)
              proportion <- frequency$Freq/sum(frequency$Freq)
              frequency <- cbind(frequency, proportion)
              max_proportion <- max(frequency$proportion)
              major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
              if(length(major_cancer) > 1){
                major_cancer <- major_cancer[1]
              }
              index_non_breast <- row.names(predicted_cancer_new)[-which(predicted_cancer_new$scATOMIC_pred %in% c("ER+ Breast Cancer Cell",
                                                                                                                   "HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell",
                                                                                                                   "Her2+ Breast Cancer Cell", "Breast Cancer Cell") )]
              seurat_object@meta.data[index_non_breast,"scATOMIC_pred"] <- major_cancer
              seurat_object@meta.data[,"pan_cancer_cluster"] <- "Normal"
              seurat_object@meta.data[row.names(predicted_cancer_new),"pan_cancer_cluster"] <- "Cancer"
            }
          }
          else{
            index_cancer <- intersect(index_cancer, outersect(index_cancer, row.names(seurat_object@meta.data)[-which(seurat_object@meta.data$scATOMIC_pred %in% cancer_classes)]))
            index_normal <- outersect(row.names(seurat_object@meta.data)[which(seurat_object@meta.data$scATOMIC_pred %in% cancer_classes)],
                                      index_cancer)
            seurat_object@meta.data[index_normal,"scATOMIC_pred"] <- "Normal Tissue Cell"
            predicted_cancer_new <- seurat_object@meta.data[which(seurat_object@meta.data$scATOMIC_pred %in%c(cancer_classes,"Breast/Lung/Prostate", "Ovarian/Endometrial/Kidney", "Endometrial Cancer Cell",
                                                                                                              "Biliary/Hepatic Cancer Cell", "Brain/Neuroblastoma Cancer Cell", "Digestive Tract Cancer Cell",
                                                                                                              "Soft Tissue Cancer Cell", "Soft Tissue or Neuro Cancer Cell", "Unclassified Soft Tissue or Neuro Cancer Cell","ER+ Breast Cancer Cell",
                                                                                                              "HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell", "Her2+ Breast Cancer Cell")),]

            if(breast_mode == F){
              frequency <- as.data.frame(table(predicted_cancer_new$scATOMIC_pred), stringsAsFactors = F)
              proportion <- frequency$Freq/sum(frequency$Freq)
              frequency <- cbind(frequency, proportion)
              max_proportion <- max(frequency$proportion)
              major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
              if(length(major_cancer) > 1){
                major_cancer <- major_cancer[1]
              }
              seurat_object@meta.data[row.names(predicted_cancer_new),"scATOMIC_pred"] <- major_cancer
              seurat_object@meta.data[,"pan_cancer_cluster"] <- "Normal"
              seurat_object@meta.data[row.names(predicted_cancer_new),"pan_cancer_cluster"] <- "Cancer"
            } else{
              frequency <- as.data.frame(table(predicted_cancer_new$layer_5), stringsAsFactors = F)
              proportion <- frequency$Freq/sum(frequency$Freq)
              frequency <- cbind(frequency, proportion)
              max_proportion <- max(frequency$proportion)
              major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
              if(length(major_cancer) > 1){
                major_cancer <- major_cancer[1]
              }
              index_non_breast <- row.names(predicted_cancer_new)[-which(predicted_cancer_new$scATOMIC_pred %in% c("ER+ Breast Cancer Cell",
                                                                                                                   "HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell",
                                                                                                                   "Her2+ Breast Cancer Cell", "Breast Cancer Cell") )]
              seurat_object@meta.data[index_non_breast,"scATOMIC_pred"] <- major_cancer
              seurat_object@meta.data[,"pan_cancer_cluster"] <- "Normal"
              seurat_object@meta.data[row.names(predicted_cancer_new),"pan_cancer_cluster"] <- "Cancer"
            }

          }
        } else{
          index_cancer <- intersect(index_cancer, outersect(index_cancer, row.names(seurat_object@meta.data)[-which(seurat_object@meta.data$scATOMIC_pred %in% cancer_classes)]))
          index_normal <- outersect(row.names(seurat_object@meta.data)[which(seurat_object@meta.data$scATOMIC_pred %in% cancer_classes)],
                                    index_cancer)
          predicted_cancer_new <- seurat_object@meta.data[which(seurat_object@meta.data$scATOMIC_pred %in%c(cancer_classes,"Breast/Lung/Prostate", "Ovarian/Endometrial/Kidney", "Endometrial Cancer Cell",
                                                                                                            "Biliary/Hepatic Cancer Cell", "Brain/Neuroblastoma Cancer Cell", "Digestive Tract Cancer Cell",
                                                                                                            "Soft Tissue Cancer Cell", "Soft Tissue or Neuro Cancer Cell", "Unclassified Soft Tissue or Neuro Cancer Cell","ER+ Breast Cancer Cell",
                                                                                                            "HER2+ Breast Cancer Cell","TNBC Breast Cancer Cell")),]
          frequency <- as.data.frame(table(predicted_cancer_new$scATOMIC_pred), stringsAsFactors = F)
          proportion <- frequency$Freq/sum(frequency$Freq)
          frequency <- cbind(frequency, proportion)
          max_proportion <- max(frequency$proportion)
          major_cancer <- frequency[which(frequency$proportion == max_proportion), "Var1"]
          if(length(major_cancer) > 1){
            major_cancer <- major_cancer[1]
          }
          seurat_object@meta.data[row.names(predicted_cancer_new),"scATOMIC_pred"] <- major_cancer
          seurat_object@meta.data[,"pan_cancer_cluster"] <- "Normal"
          seurat_object@meta.data[row.names(predicted_cancer_new),"pan_cancer_cluster"] <- "Cancer"

          summary_master <- seurat_object@meta.data
        }

        summary_master <- seurat_object@meta.data

      }
      summary_master <- seurat_object@meta.data
      summary_master$classification_confidence <- "confident"

      if(nrow(predicted_cancer) > 0){
        median_score_major_cancer <- min(c(summary_master[which(summary_master$layer_6 == major_cancer), "median_score_class_layer_4"],
                                           summary_master[which(summary_master$layer_6 == major_cancer), "median_score_class_layer_5"],
                                           summary_master[which(summary_master$layer_6 == major_cancer), "median_score_class_layer_6"]))
        summary_master$classification_confidence <- "confident"
        if(cancer_confidence == "default"){
          index_non_confident <- row.names(summary_master)[which((summary_master$layer_1 %in% c("Tissue_Cell_Normal_or_Cancer") & summary_master$median_score_class_layer_2 < 0.55 ) |
                                                                   (summary_master$scATOMIC_pred == "Brain Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.6 | summary_master$median_score_class_layer_5 < 0.62 |summary_master$median_score_class_layer_6 < 0.8 ) )|
                                                                   (summary_master$scATOMIC_pred == "Breast Cancer Cell" & (summary_master$median_score_class_layer_5 < 0.55) )|
                                                                   (summary_master$scATOMIC_pred == "Colorectal Cancer Cell" & (summary_master$median_score_class_layer_4 < 0.7 |summary_master$median_score_class_layer_5 < 0.68  ) )|
                                                                   (summary_master$scATOMIC_pred == "Gastric Cancer Cell" & (summary_master$median_score_class_layer_4 < 0.7 |summary_master$median_score_class_layer_5 < 0.68  ) )|
                                                                   (summary_master$scATOMIC_pred == "Colon/Colorectal Cancer Cell" & (summary_master$median_score_class_layer_4 < 0.7 |summary_master$median_score_class_layer_5 < 0.68  ) )|
                                                                   (summary_master$scATOMIC_pred == "Endometrial Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.6|summary_master$median_score_class_layer_4 < 0.7|summary_master$median_score_class_layer_5 < 0.5   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Endometrial/Uterine Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.6|summary_master$median_score_class_layer_4 < 0.7|summary_master$median_score_class_layer_5 < 0.5   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Kidney Cancer Cell" & ( summary_master$median_score_class_layer_4 < 0.7) )|
                                                                   (summary_master$scATOMIC_pred == "Liver Cancer Cell" & (summary_master$median_score_class_layer_5 < 0.69   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Lung Cancer Cell" & (summary_master$median_score_class_layer_4 < 0.59 |summary_master$median_score_class_layer_5 < 0.6    ) ) |
                                                                   (summary_master$scATOMIC_pred == "Neuroblastoma Cell" & (summary_master$median_score_class_layer_6 < 0.55   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Ovarian Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.6|summary_master$median_score_class_layer_4 < 0.7|summary_master$median_score_class_layer_5 < 0.5   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Pancreatic Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.59 |  summary_master$median_score_class_layer_4 < 0.55   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Bile Duct Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.59 |  summary_master$median_score_class_layer_4 < 0.55   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Skin Cancer Cell" & (summary_master$median_score_class_layer_3 < 0.6   ) ) |
                                                                   (summary_master$scATOMIC_pred == "Prostate Cancer Cell" & (summary_master$median_score_class_layer_5 < 0.68 )))]
          summary_master[index_non_confident, "classification_confidence"] <-"low_confidence"
          index_new_major_cancer <- which(summary_master$pan_cancer_cluster == "Cancer" & summary_master$classification_confidence == "confident")
          if(length(index_new_major_cancer) > 0){
            major_cancer_fixed <- names(sort(table(summary_master[index_new_major_cancer, "layer_6"]), decreasing = T))[1]
            summary_master[which(summary_master$pan_cancer_cluster == "Cancer"), "scATOMIC_pred"] <- major_cancer_fixed
          }
          if(length(index_non_confident) > 0 ){
            summary_master[index_non_confident, "classification_confidence"] <-"low_confidence"
          }
          if(length(index_non_confident) > 0 ){
            preds_layer_6 <- summary_master[index_non_confident, "layer_6"]
            if(major_cancer %in% preds_layer_6){
              warning("Level of confidence for cancer typing is low")
            }
          }

          index_esophageal <- row.names(summary_master)[which(summary_master$scATOMIC_pred == "Esophageal Cancer Cell")]
          if(length(index_esophageal) > 0){
            warning("Cancer was called esophageal, but may be a different squamous cell cancer type")
            summary_master[index_esophageal, "classification_confidence"] <-"low_confidence"
          }
          index_na_confident <- row.names(summary_master)[which(summary_master$scATOMIC_pred %in%  c("Bile Duct Cancer Cell",
                                                                                                     "Bone Cancer Cell",
                                                                                                     "Endometrial/Uterine Cancer Cell",
                                                                                                     "Endometrial Cancer Cell",
                                                                                                     "Gallbladder Cancer Cell",
                                                                                                     "Gastric Cancer Cell")) ]
          if(length(index_na_confident) > 0){
            warning("Cancer was not externally validated, no confidence assigned")
            summary_master[index_na_confident, "classification_confidence"] <- NA
          }
        } else if(is.numeric(cancer_confidence) & median_score_major_cancer < cancer_confidence){
          warning(paste0("Level of confidence for cancer typing = ", median_score_major_cancer, "\nThis classfication is not reliable"))
          summary_master[which(summary_master$scATOMIC_pred == major_cancer), "classification_confidence"] <-"low_confidence"
        }

      }
      if(cancer_confidence == "default"){
        index_non_confident_blood_stromal <- row.names(summary_master)[which(
          (summary_master$scATOMIC_pred == "B cell" & summary_master$median_score_class_layer_3 < 0.7) |
            (summary_master$layer_4 == "CD4+ T cell" & (summary_master$median_score_class_layer_2 < 0.6 | summary_master$median_score_class_layer_3 < 0.9 |summary_master$median_score_class_layer_4 < 0.9 ) )|
            (summary_master$layer_4 == "CD8+ T cell" & (summary_master$median_score_class_layer_2 < 0.6 | summary_master$median_score_class_layer_3 < 0.9 |summary_master$median_score_class_layer_4 < 0.6 ) )|
            (summary_master$scATOMIC_pred == "cDC" & (summary_master$median_score_class_layer_2 < 0.9 ) )|
            (summary_master$scATOMIC_pred == "Macrophage" & (summary_master$median_score_class_layer_2 < 0.9 | summary_master$median_score_class_layer_3 < 0.66|summary_master$median_score_class_layer_4 < 0.75 ) ) |
            (summary_master$scATOMIC_pred == "Natural killer cell" & ( summary_master$median_score_class_layer_3 < 0.9) )|
            (summary_master$scATOMIC_pred == "pDC" & (summary_master$median_score_class_layer_2 < 0.9   ) ) |
            (summary_master$scATOMIC_pred == "Plasmablast" & (summary_master$median_score_class_layer_2 < 0.8 | summary_master$median_score_class_layer_3 < 0.9 )
            ))]
        summary_master[index_non_confident_blood_stromal, "classification_confidence"] <- "low_confidence"

      }



    }
    if(normal_tissue == T){
      index_potential_cancer <- which(summary_master$scATOMIC_pred %in% c("Bile Duct Cancer Cell","Bladder Cancer Cell",
                                                                          "Bone Cancer Cell","Brain Cancer Cell"    ,
                                                                          "Breast Cancer Cell","Colon/Colorectal Cancer Cell"   ,
                                                                          "Endometrial/Uterine Cancer Cell","Esophageal Cancer Cell"   ,
                                                                          "Gallbladder Cancer Cell","Gastric Cancer Cell"     ,
                                                                          "Kidney Cancer Cell","Liver Cancer Cell"   ,
                                                                          "Lung Cancer Cell","Neuroblastoma"       ,
                                                                          "Ovarian Cancer Cell","Pancreatic Cancer Cell"  ,
                                                                          "Prostate Cancer Cell","Sarcoma"     ,
                                                                          "Skin Cancer Cell","Breast/Lung/Prostate"    ,
                                                                          "Ovarian/Endometrial/Kidney","Endometrial Cancer Cell"    ,
                                                                          "Biliary/Hepatic Cancer Cell","Brain/Neuroblastoma Cancer Cell",
                                                                          "Digestive Tract Cancer Cell","Soft Tissue Cancer Cell" ,
                                                                          "Soft Tissue or Neuro Cancer Cell","Unclassified Soft Tissue or Neuro Cancer Cell",
                                                                          "ER+ Breast Cancer Cell","HER2+ Breast Cancer Cell" ,
                                                                          "Her2+ Breast Cancer Cell","TNBC Breast Cancer Cell", "Epithelial Cell","GI Tract Cell",
                                                                          "Non Blood Cell","Non Stromal Cell"))
      if(length(index_potential_cancer) > 0){
        summary_master[index_potential_cancer, "scATOMIC_pred"] <- "Normal Tissue Cell"
      }
      index_CAF <- which(summary_master$layer_3 %in% c("Cancer Associated Fibroblasts"))
      if(length(index_CAF) > 0){
        summary_master[index_CAF, "scATOMIC_pred"] <- "Normal Fibroblast"
      }

    }
    if(low_res_mode == T){
      #index CD8 T cells
      index_CD8 <- which(summary_master$layer_4 %in% c("CD8+ T cell"))
      if(length(index_CD8) > 0){
        summary_master[index_CD8, "scATOMIC_pred"] <- "CD8+ T cell"
      }
      #index CD4 T cells
      index_CD4 <- which(summary_master$layer_4 %in% c("CD4+ T cell"))
      if(length(index_CD4) > 0){
        summary_master[index_CD4, "scATOMIC_pred"] <- "CD4+ T cell"
      }
      #index macrophages
      index_macro <- which(summary_master$layer_4 %in% c("Macrophage"))
      if(length(index_macro) > 0){
        summary_master[index_macro, "scATOMIC_pred"] <- "Macrophage"
      }
      #index fibroblasts
      index_CAF <- which(summary_master$layer_3 %in% c("Cancer Associated Fibroblasts") & summary_master$scATOMIC_pred != "Normal Fibroblast")
      if(length(index_CAF) > 0){
        summary_master[index_CAF, "scATOMIC_pred"] <- "Cancer Associated Fibroblasts"
      }
      #index cDCs
      index_cDC <- which(summary_master$layer_4 %in% c("cDC"))
      if(length(index_cDC) > 0){
        summary_master[index_cDC, "scATOMIC_pred"] <- "cDC"
      }
    }
    percent_confident <- round(length(which(summary_master$classification_confidence == "confident"))/nrow(summary_master), digits = 2)
    if(percent_confident == 1){
      print(paste0("Sample classification confidence = 1.00"))
    } else{
      print(paste0("Sample classification confidence = ", percent_confident))
    }
    if(percent_confident < 0.75){
      warning("Sample classification confidence is low")
    }
    if(length(known_cancer_type) == 1){
      index_cancer_predicted <- which(summary_master$scATOMIC_pred %in% c("Bile Duct Cancer Cell", "Bladder Cancer Cell", "Bone Cancer Cell",
                                                                          "Brain Cancer Cell", "Breast Cancer Cell", "Colon/Colorectal Cancer Cell",
                                                                          "Endometrial/Uterine Cancer Cell", "Esophageal Cancer Cell", "Gallbladder Cancer Cell",
                                                                          "Gastric Cancer Cell", "Kidney Cancer Cell", "Liver Cancer Cell", "Lung Cancer Cell",
                                                                          "Neuroblastoma", "Ovarian Cancer Cell", "Pancreatic Cancer Cell", "Prostate Cancer Cell",
                                                                          "Sarcoma", "Skin Cancer Cell", "Breast/Lung/Prostate", "Ovarian/Endometrial/Kidney",
                                                                          "Endometrial Cancer Cell", "Biliary/Hepatic Cancer Cell", "Brain/Neuroblastoma Cancer Cell",
                                                                          "Digestive Tract Cancer Cell", "Soft Tissue Cancer Cell", "Soft Tissue or Neuro Cancer Cell",
                                                                          "Unclassified Soft Tissue or Neuro Cancer Cell", "ER+ Breast Cancer Cell", "HER2+ Breast Cancer Cell",
                                                                          "Her2+ Breast Cancer Cell", "TNBC Breast Cancer Cell", "Epithelial Cell", "GI Tract Cell"))
      summary_master[index_cancer_predicted, "scATOMIC_pred"] <- known_cancer_type
    }
    return(summary_master)
  }
}





