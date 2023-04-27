#' Run scATOMIC
#'
#' @param rna_counts count matrix
#' @param imputation whether to apply MAGIC imputation - recommended
#' @param ref_based whether to use ref based - currently depreceated
#' @param mc.cores number of cores
#' @param unimodal_nsd number of sd for setting threshold in unimodal distributions
#' @param bimodal_nsd number of sd for setting threshold in bimodal distributions
#' @param breast_mode run breast subclassfication
#' @param confidence_cutoff logical to use confidence cutoffs
#' @param fine_grained_T run fine grained T cell subclassification
#'
#' @return list of layers classified
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
#' }
run_scATOMIC <- function(rna_counts, imputation = TRUE, ref_based = F, mc.cores = 1, unimodal_nsd = 3, bimodal_nsd = 2, breast_mode = F, confidence_cutoff = T, fine_grained_T = T){
  if(.Platform$OS.type == "windows"){
    mc.cores = 1
  }
  if(confidence_cutoff == T){
    prediction_list <- list()
    print("Starting Layer 1")
    normalized_counts <- Rmagic::library.size.normalize(t(as.matrix(rna_counts)))
    normalized_counts <- t(sqrt(normalized_counts))
    prediction_list[["layer_1"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = colnames(rna_counts),
                                                                       layer = "layer_1", imputation = imputation, genes_in_model = top_genes_unlisted_layer_1,
                                                                       model = model_layer_1, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                       bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

    print("Done Layer 1")
    normal_tissue_cancer_predicted <- row.names(prediction_list[["layer_1"]])[which(prediction_list[["layer_1"]]$predicted_tissue_with_cutoff %in%
                                                                                      c( "Bile Duct Cancer","Bladder Cancer",  "Bone Cancer",     "Brain Cancer",    "Breast Cancer",
                                                                                         "Colon/Colorectal Cancer","Endometrial/Uterine Cancer","Esophageal Cancer",
                                                                                         "Gallbladder Cancer", "Gastric Cancer", "Glial Cells",  "Kidney Cancer", "Liver Cancer", "Lung Cancer",
                                                                                         "Neuroblastoma",  "Oligodendrocytes","Ovarian Cancer",
                                                                                         "Pancreatic Cancer", "Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts",
                                                                                         "Prostate Cancer", "Skin Cancer", "Endothelial Cells", "Fibroblasts", "Myofibroblasts",  "Smooth Muscle Cells",
                                                                                         "Sarcoma", "Tissue_Cell_Normal_or_Cancer"))]
    blood_cells_predicted <- row.names(prediction_list[["layer_1"]])[which(prediction_list[["layer_1"]]$predicted_tissue_with_cutoff %in% c("B cell","CD4+ T cell", "CD8+ T cell",  "Mast cell", "macrophage_or_dendritic_cell",
                                                                                                                                            "Natural killer cell", "Blood_Cell"))]
    if (length(normal_tissue_cancer_predicted) > 0){
      print("Starting Layer 2 Non Blood")
      prediction_list[["layer_2_non_blood"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = normal_tissue_cancer_predicted,
                                                                                   layer = "layer_2_non_blood", imputation = imputation, genes_in_model = top_genes_unlisted_layer_2_normal_tissue_cancer,
                                                                                   model = model_layer_2_normal_tissue_cancer, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                   bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
      print("Done Layer 2 Non Blood")
      non_stromal_predicted <- row.names(prediction_list[["layer_2_non_blood"]])[which(prediction_list[["layer_2_non_blood"]]$predicted_tissue_with_cutoff %in% c(
        "Bile Duct Cancer","Bladder Cancer",  "Bone Cancer","Brain Cancer","Breast Cancer",
        "Colon/Colorectal Cancer","Endometrial/Uterine Cancer","Esophageal Cancer",
        "Gallbladder Cancer", "Gastric Cancer",   "Kidney Cancer", "Liver Cancer", "Lung Cancer",
        "Neuroblastoma",  "Ovarian Cancer",
        "Pancreatic Cancer",
        "Prostate Cancer", "Skin Cancer",
        "Sarcoma", "Non Stromal Cell") )]
      normal_stromal_predicted <- row.names(prediction_list[["layer_2_non_blood"]])[which(prediction_list[["layer_2_non_blood"]]$predicted_tissue_with_cutoff %in%
                                                                                            c("Stromal Cell","Endothelial cell","Fibroblasts", "Myofibroblasts","Smooth Muscle Cells",
                                                                                              "Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts",
                                                                                              "fibroblasts_score"))]
      if (length(non_stromal_predicted) > 0){
        print("Starting Layer 3 Non Stromal")
        prediction_list[["layer_3_non_stromal"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = non_stromal_predicted,
                                                                                       layer = "layer_3_non_stromal", imputation = imputation, genes_in_model = top_genes_unlisted_layer_3_non_stromal,
                                                                                       model = model_layer_3_non_stromal, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                       bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
        print("Done Layer 3 Non Stromal")
        non_GI_predicted <- row.names(prediction_list[["layer_3_non_stromal"]])[
          which(prediction_list[["layer_3_non_stromal"]]$predicted_tissue_with_cutoff %in%c(
            "Breast Cancer",
            "Endometrial/Uterine Cancer", "Lung Cancer", "Ovarian Cancer",
            "Prostate Cancer", "Kidney Cancer",  "Non GI Epithelial Cell"))]
        soft_tissue_neuro_cancer_predicted <- row.names(prediction_list[["layer_3_non_stromal"]])[
          which(prediction_list[["layer_3_non_stromal"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer","Brain Cancer",
                                                                                             "Neuroblastoma", "Skin Cancer",
                                                                                             "Sarcoma", "Soft Tissue or Neuro Cancer Cell") )]
        GI_predicted <- row.names(prediction_list[["layer_3_non_stromal"]])[
          which(prediction_list[["layer_3_non_stromal"]]$predicted_tissue_with_cutoff %in% c("Bile Duct Cancer","Bladder Cancer",
                                                                                             "Colon/Colorectal Cancer","Esophageal Cancer",
                                                                                             "Gallbladder Cancer", "Gastric Cancer", "Liver Cancer",
                                                                                             "Pancreatic Cancer","GI Epithelial Cell") )]
        if (length(non_GI_predicted) > 0){
          print("Starting Layer 4 Non GI")
          prediction_list[["layer_4_non_GI"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = non_GI_predicted,
                                                                                    layer = "layer_4_non_GI", imputation = imputation, genes_in_model = top_genes_unlisted_layer_4_non_GI,
                                                                                    model = model_layer_4_non_GI, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                    bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
          print("Done Layer 4 Non GI")
          breast_lung_prostate_predicted <- row.names(prediction_list[["layer_4_non_GI"]] )[
            which(prediction_list[["layer_4_non_GI"]]$predicted_tissue_with_cutoff  %in%c(
              "Breast Cancer","Lung Cancer",
              "Prostate Cancer",  "Breast/Lung/Prostate Cell"))]
          ov_endo_kid_predicted <- row.names(prediction_list[["layer_4_non_GI"]])[
            which(prediction_list[["layer_4_non_GI"]]$predicted_tissue_with_cutoff %in%c(
              "Endometrial/Uterine Cancer", "Ovarian Cancer", "Kidney Cancer", "Ovarian/Endometrial/Kidney Cell"))]
          if (length(breast_lung_prostate_predicted) > 0){
            print("Starting Layer 5 Breast Lung Prostate")
            prediction_list[["layer_5_breast_lung_prostate"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = breast_lung_prostate_predicted,
                                                                                                    layer = "layer_5_breast_lung_prostate", imputation = imputation,
                                                                                                    genes_in_model = top_genes_unlisted_layer_5_breast_lung_prostate,
                                                                                                    model = model_layer_5_breast_lung_prostate,
                                                                                                    ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                    bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
            print("Done Layer 5 Breast Lung Prostate")
            if(breast_mode == T){
              breast_predicted <- row.names(prediction_list[["layer_5_breast_lung_prostate"]])[which(prediction_list[["layer_5_breast_lung_prostate"]]$predicted_tissue_with_cutoff == "Breast Cancer Cell")]
              if (length(breast_predicted) > 0){
                print("Starting Layer 6 Breast")
                prediction_list[["layer_6_breast"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = breast_predicted,
                                                                                          layer = "layer_6_breast", imputation = imputation,
                                                                                          genes_in_model = top_genes_unlisted_layer_6_breast,
                                                                                          model = model_layer_6_breast,
                                                                                          ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                          bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
                print("Done Layer 6 Breast")
              }
            }

          }
          if (length(ov_endo_kid_predicted) > 0){
            print("Starting Layer 5 Ovarian Endometrial Kidney")
            prediction_list[["layer_5_ov_endo_kid"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = ov_endo_kid_predicted,
                                                                                           layer = "layer_5_ov_endo_kid", imputation = imputation,
                                                                                           genes_in_model = top_genes_unlisted_layer_5_ov_endo_kid,
                                                                                           model = model_layer_5_ov_endo_kid,
                                                                                           ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                           bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
            print("Done Layer 5 Ovarian Endometrial Kidney")
          }

        }
        if (length(GI_predicted) > 0){
          print("Starting Layer 4 GI")
          prediction_list[["layer_4_GI"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = GI_predicted,
                                                                                layer = "layer_4_GI", imputation = imputation, genes_in_model = top_genes_unlisted_layer_4_GI,
                                                                                model = model_layer_4_GI, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
          print("Done Layer 4 GI")
          digestive_predicted <- row.names(prediction_list[["layer_4_GI"]])[
            which(prediction_list[["layer_4_GI"]]$predicted_tissue_with_cutoff %in% c("Colon/Colorectal Cancer","Esophageal Cancer",
                                                                                      "Gastric Cancer","Colorectal/Esophageal/Gastric Cell"))]
          biliary_predicted <- row.names(prediction_list[["layer_4_GI"]])[
            which(prediction_list[["layer_4_GI"]]$predicted_tissue_with_cutoff %in% c("Bile Duct Cancer","Bladder Cancer",
                                                                                      "Gallbladder Cancer",  "Liver Cancer",
                                                                                      "Pancreatic Cancer", "Billiary Cell"))]

          if (length(digestive_predicted) > 0){
            print("Starting Layer 5 Digestive")
            prediction_list[["layer_5_digestive"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = digestive_predicted,
                                                                                         layer = "layer_5_digestive", imputation = imputation,
                                                                                         genes_in_model = top_genes_unlisted_layer_5_digestive,
                                                                                         model = model_layer_5_digestive,
                                                                                         ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                         bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
            print("Done Layer 5 Digestive")
          }
          if (length(biliary_predicted) > 0){
            print("Starting Layer 5 Biliary")
            prediction_list[["layer_5_biliary"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = biliary_predicted,
                                                                                       layer = "layer_5_biliary", imputation = imputation,
                                                                                       genes_in_model = top_genes_unlisted_layer_5_biliary,
                                                                                       model = model_layer_5_biliary,
                                                                                       ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                       bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
            print("Done Layer 5 Biliary")
          }
        }
        if (length(soft_tissue_neuro_cancer_predicted) > 0){
          print("Starting Layer 4 Soft Tissue Neuro")
          prediction_list[["layer_4_soft_tissue_neuro"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = soft_tissue_neuro_cancer_predicted,
                                                                                               layer = "layer_4_soft_tissue_neuro", imputation = imputation,
                                                                                               genes_in_model = top_genes_unlisted_layer_4_soft_neuro_cancer,
                                                                                               model = model_layer_4_soft_neuro_cancer,
                                                                                               ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                               bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
          print("Done Layer 4 Soft Tissue Neuro")
          soft_tissue_neuro_cancer_no_lung_skin_predicted <- row.names(prediction_list[["layer_4_soft_tissue_neuro"]])[
            which(prediction_list[["layer_4_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer","Brain Cancer",
                                                                                                     "Neuroblastoma",
                                                                                                     "Sarcoma", "Soft Tissue or Neuro Cancer Cell", "Soft Tissue or Neuro Cancer Cell"))]
          if (length(soft_tissue_neuro_cancer_no_lung_skin_predicted) > 0){
            print("Starting Layer 5 Soft Tissue Neuro")
            prediction_list[["layer_5_soft_tissue_neuro"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = soft_tissue_neuro_cancer_no_lung_skin_predicted,
                                                                                                 layer = "layer_5_soft_tissue_neuro", imputation = imputation,
                                                                                                 genes_in_model = top_genes_unlisted_layer_5_soft_neuro_cancer_no_lung_skin,
                                                                                                 model = model_layer_5_soft_neuro_cancer_no_lung_skin,
                                                                                                 ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                 bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
            print("Done Layer 5 Soft Tissue Neuro")
            brain_nbm_predicted <- row.names(prediction_list[["layer_5_soft_tissue_neuro"]])[
              which(prediction_list[["layer_5_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Brain Cancer",
                                                                                                       "Neuroblastoma", "Brain/Neuroblastoma Cancer Cell"))]
            soft_tissue_predicted <- row.names(prediction_list[["layer_5_soft_tissue_neuro"]])[which(prediction_list[["layer_5_soft_tissue_neuro"]]$predicted_tissue_with_cutoff %in% c("Bone Cancer",
                                                                                                                                                                                        "Sarcoma", "Soft Tissue Cancer Cell"))]
            if (length(brain_nbm_predicted) > 0){
              print("Starting Layer 6 Brain NBM")
              prediction_list[["layer_6_brain_nbm"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = brain_nbm_predicted,
                                                                                           layer = "layer_6_brain_nbm", imputation = imputation,
                                                                                           genes_in_model = top_genes_unlisted_layer_6_brain_nbm,
                                                                                           model = model_layer_6_brain_nbm,
                                                                                           ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                           bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
              print("Done Layer 6 Brain NBM")
            }
            if (length(soft_tissue_predicted) > 0){
              print("Starting Layer 6 Soft Tissue")
              prediction_list[["layer_6_soft_tissue"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = soft_tissue_predicted,
                                                                                             layer = "layer_6_soft_tissue", imputation = imputation,
                                                                                             genes_in_model = top_genes_unlisted_layer_6_soft_tissue,
                                                                                             model = model_layer_6_soft_tissue,
                                                                                             ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                             bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
              print("Done Layer 6 Soft Tissue")
            }

          }
        }
      }
      if (length(normal_stromal_predicted) > 0){
        print("Starting Layer 3 Normal Stromal")
        prediction_list[["layer_3_stromal"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = normal_stromal_predicted,
                                                                                   layer = "layer_3_stromal", imputation = imputation, genes_in_model = top_genes_unlisted_layer_3_stromal,
                                                                                   model = model_layer_3_stromal, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                   bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
        print("Done Layer 3 Normal Stromal")
        CAF_predicted <- row.names(prediction_list[["layer_3_stromal"]])[
          which(prediction_list[["layer_3_stromal"]]$predicted_tissue_with_cutoff == "Cancer Associated Fibroblasts")]
        if (length(CAF_predicted) > 0){
          print("Starting Layer 4 CAFs")
          prediction_list[["layer_4_CAF"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = CAF_predicted,
                                                                                  layer = "layer_4_CAF", imputation = F,
                                                                                  genes_in_model = top_genes_unlisted_layer_4_CAFs,
                                                                                  model = model_layer_4_CAFs,
                                                                                  ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                  bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

          print("Done Layer 4 CAFs")
        }
      }
    }
    if (length(blood_cells_predicted) > 0){
      print("Starting Layer 2 Blood")
      prediction_list[["layer_2_blood"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = blood_cells_predicted,
                                                                               layer = "layer_2_blood", imputation = imputation, genes_in_model = top_genes_unlisted_layer_2_blood,
                                                                               model = model_layer_2_blood, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                               bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

      print("Done Layer 2 Blood")
      TNK_cells_predicted <- row.names(prediction_list[["layer_2_blood"]])[
        which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell", "CD8+ T cell", "Natural killer cell", "T_or_NK_lymphocyte"))]
      MDC_cells_predicted <- row.names(prediction_list[["layer_2_blood"]])[
        which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff %in% c("macrophage_or_dendritic_cell", "macrophage_DC_score"))]
      B_cells_predicted <- row.names(prediction_list[["layer_2_blood"]])[
        which(prediction_list[["layer_2_blood"]]$predicted_tissue_with_cutoff %in% c("B cell","Plasmablast" ,"B_cell_score","B cell or Plasmablast"))]
      if (length(TNK_cells_predicted) > 0){
        print("Starting Layer 3 TNK")
        prediction_list[["layer_3_TNK"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = TNK_cells_predicted,
                                                                               layer = "layer_3_TNK", imputation = F, genes_in_model = top_genes_unlisted_layer_3_TNK,
                                                                               model = model_layer_3_TNK, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                               bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

        print("Done Layer 3 TNK")
        CD4_CD8_cells_predicted <- row.names(prediction_list[["layer_3_TNK"]])[
          which(prediction_list[["layer_3_TNK"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell","CD8+ T cell" ,"CD4 or CD8 T cell"))]
        CD8_NK_cells_predicted <- row.names(prediction_list[["layer_3_TNK"]])[
          which(prediction_list[["layer_3_TNK"]]$predicted_tissue_with_cutoff %in% c("Natural killer cell", "NK or CD8 T cell"))]

        if (length(CD4_CD8_cells_predicted) > 0){
          print("Starting Layer 4 CD4 CD8")
          prediction_list[["layer_4_CD4_CD8"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = CD4_CD8_cells_predicted,
                                                                                     layer = "layer_4_CD4_CD8", imputation = F,
                                                                                     genes_in_model = top_genes_unlisted_layer_4_CD4_CD8,
                                                                                     model = model_layer_4_CD4_CD8,
                                                                                     ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                     bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

          print("Done Layer 4 CD4 CD8")
        }
        if (length(CD8_NK_cells_predicted) > 0){
          print("Starting Layer 4 CD8 NK")
          prediction_list[["layer_4_CD8_NK"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = CD8_NK_cells_predicted,
                                                                                    layer = "layer_4_CD8_NK", imputation = F,
                                                                                    genes_in_model = top_genes_unlisted_layer_4_CD8_NK,
                                                                                    model = model_layer_4_CD8_NK,
                                                                                    ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                    bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

          print("Done Layer 4 CD8 NK")
        }
        if(fine_grained_T == T){
          CD4_cells_predicted <- row.names(prediction_list[["layer_4_CD4_CD8"]])[
            which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell"))]
          CD8_cells_predicted <- c(row.names(prediction_list[["layer_4_CD4_CD8"]])[
            which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell"))],row.names(prediction_list[["layer_4_CD8_NK"]])[
              which(prediction_list[["layer_4_CD8_NK"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell"))] )
          if (length(CD4_cells_predicted) > 0){
            print("Starting Layer 5 CD4")
            prediction_list[["layer_5_CD4"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = CD4_cells_predicted,
                                                                                   layer = "layer_5_CD4", imputation = F,
                                                                                   genes_in_model = top_genes_unlisted_layer_5_CD4,
                                                                                   model = model_layer_5_CD4,
                                                                                   ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                   bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

            print("Done Layer 5 CD4")
          }
          if (length(CD8_cells_predicted) > 0){
            print("Starting Layer 5 CD8")
            prediction_list[["layer_5_CD8"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = CD8_cells_predicted,
                                                                                   layer = "layer_5_CD8", imputation = F,
                                                                                   genes_in_model = top_genes_unlisted_layer_5_CD8,
                                                                                   model = model_layer_5_CD8,
                                                                                   ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                   bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

            print("Done Layer 5 CD8")
          }






        }


      }
      if (length(MDC_cells_predicted) > 0){
        print("Starting Layer 3 Myeloid")
        prediction_list[["layer_3_myeloid"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = MDC_cells_predicted,
                                                                                   layer = "layer_3_MDC", imputation = imputation, genes_in_model = top_genes_unlisted_layer_3_MDC,
                                                                                   model = model_layer_3_MDC, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                   bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

        print("Done Layer 3 Myeloid")
        DC_cells_predicted <- row.names(prediction_list[["layer_3_myeloid"]])[
          which(prediction_list[["layer_3_myeloid"]]$predicted_tissue_with_cutoff == "Dendritic Cell")]
        macrophage_monocyte_predicted <- row.names(prediction_list[["layer_3_myeloid"]])[
          which(prediction_list[["layer_3_myeloid"]]$predicted_tissue_with_cutoff == "Macrophage or Monocyte")]
        if (length(DC_cells_predicted) > 0){
          print("Starting Layer 4 Dendritic")
          prediction_list[["layer_4_dendritic"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = DC_cells_predicted,
                                                                                       layer = "layer_4_dendritic", imputation = imputation,
                                                                                       genes_in_model = top_genes_unlisted_layer_4_DC_cell,
                                                                                       model = model_layer_4_DC_cell,
                                                                                       ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                       bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

          print("Done Layer 4 Dendritic")
          cDC_predicted <- row.names(prediction_list[["layer_4_dendritic"]])[
            which(prediction_list[["layer_4_dendritic"]]$predicted_tissue_with_cutoff == "cDC")]
          if (length(cDC_predicted) > 0){
            print("Starting Layer 5 cDC")
            prediction_list[["layer_5_cDC"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = cDC_predicted,
                                                                                   layer = "layer_5_cDC", imputation = F,
                                                                                   genes_in_model = top_genes_unlisted_layer_5_cDC,
                                                                                   model = model_layer_5_cDC,
                                                                                   ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                   bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

            print("Done Layer 5 cDC")
          }

        }
        if (length(macrophage_monocyte_predicted) > 0){
          print("Starting Layer 4 Macrophage Monocyte")
          prediction_list[["layer_4_macrophage"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = macrophage_monocyte_predicted,
                                                                                        layer = "layer_4_macrophage", imputation = imputation,
                                                                                        genes_in_model = top_genes_unlisted_layer_4_macrophage_monocyte,
                                                                                        model = model_layer_4_macrophage_monocyte,
                                                                                        ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                        bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

          print("Done Layer 4 Macrophage Monocyte")
          macrophage_predicted <- row.names(prediction_list[["layer_4_macrophage"]])[
            which(prediction_list[["layer_4_macrophage"]]$predicted_tissue_with_cutoff == "Macrophage")]
          monocyte_predicted <- row.names(prediction_list[["layer_4_macrophage"]])[
            which(prediction_list[["layer_4_macrophage"]]$predicted_tissue_with_cutoff == "Monocyte")]
          if (length(macrophage_predicted) > 0){
            print("Starting Layer 5 Macrophage")
            prediction_list[["layer_5_macrophage"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = macrophage_predicted,
                                                                             layer = "layer_5_macrophage", imputation = F,
                                                                             genes_in_model = top_genes_unlisted_layer_5_macrophage,
                                                                             model = model_layer_5_macrophage,
                                                                             ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                             bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

            print("Done Layer 5 Macrophage")
          }
          if (length(monocyte_predicted) > 0){
            print("Starting Layer 5 Monocyte")
            prediction_list[["layer_5_monocyte"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = monocyte_predicted,
                                                                                    layer = "layer_5_monocyte", imputation = F,
                                                                                    genes_in_model = top_genes_unlisted_layer_5_monocyte,
                                                                                    model = model_layer_5_monocyte,
                                                                                    ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                    bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

            print("Done Layer 5 Monocyte")
          }
        }

      }
      if (length(B_cells_predicted) > 0){
        print("Starting Layer 3 B Cell")
        prediction_list[["layer_3_BCell"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = B_cells_predicted,
                                                                                 layer = "layer_3_BCell", imputation = imputation, genes_in_model = top_genes_unlisted_layer_3_BCell,
                                                                                 model = model_layer_3_BCell, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                 bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

        print("Done Layer 3 B Cell")
      }
    }
    return(prediction_list)
  } else{
    prediction_list <- list()
    print("Starting Layer 1")
    normalized_counts <- Rmagic::library.size.normalize(t(as.matrix(rna_counts)))
    normalized_counts <- t(sqrt(normalized_counts))
    prediction_list[["layer_1"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = colnames(rna_counts),
                                                                                 layer = "layer_1", imputation = imputation, genes_in_model = top_genes_unlisted_layer_1,
                                                                                 model = model_layer_1, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                 bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

    print("Done Layer 1")
    normal_tissue_cancer_predicted <- row.names(prediction_list[["layer_1"]])[which(prediction_list[["layer_1"]]$predicted_class %in%
                                                                                      c( "Bile Duct Cancer","Bladder Cancer",  "Bone Cancer",     "Brain Cancer",    "Breast Cancer",
                                                                                         "Colon/Colorectal Cancer","Endometrial/Uterine Cancer","Esophageal Cancer",
                                                                                         "Gallbladder Cancer", "Gastric Cancer", "Glial Cells",  "Kidney Cancer", "Liver Cancer", "Lung Cancer",
                                                                                         "Neuroblastoma",  "Oligodendrocytes","Ovarian Cancer",
                                                                                         "Pancreatic Cancer", "Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts",
                                                                                         "Prostate Cancer", "Skin Cancer", "Endothelial Cells", "Fibroblasts", "Myofibroblasts",  "Smooth Muscle Cells",
                                                                                         "Sarcoma", "Tissue_Cell_Normal_or_Cancer"))]
    blood_cells_predicted <- row.names(prediction_list[["layer_1"]])[which(prediction_list[["layer_1"]]$predicted_class %in% c("B cell","CD4+ T cell", "CD8+ T cell",  "Mast cell", "macrophage_or_dendritic_cell",
                                                                                                                               "Natural killer cell", "Blood_Cell"))]
    if (length(normal_tissue_cancer_predicted) > 0){
      print("Starting Layer 2 Non Blood")
      prediction_list[["layer_2_non_blood"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = normal_tissue_cancer_predicted,
                                                                                             layer = "layer_2_non_blood", imputation = imputation, genes_in_model = top_genes_unlisted_layer_2_normal_tissue_cancer,
                                                                                             model = model_layer_2_normal_tissue_cancer, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                             bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
      print("Done Layer 2 Non Blood")
      non_stromal_predicted <- row.names(prediction_list[["layer_2_non_blood"]])[which(prediction_list[["layer_2_non_blood"]]$predicted_class %in% c(
        "Bile Duct Cancer","Bladder Cancer",  "Bone Cancer","Brain Cancer","Breast Cancer",
        "Colon/Colorectal Cancer","Endometrial/Uterine Cancer","Esophageal Cancer",
        "Gallbladder Cancer", "Gastric Cancer",   "Kidney Cancer", "Liver Cancer", "Lung Cancer",
        "Neuroblastoma",  "Ovarian Cancer",
        "Pancreatic Cancer",
        "Prostate Cancer", "Skin Cancer",
        "Sarcoma", "Non Stromal Cell") )]
      normal_stromal_predicted <- row.names(prediction_list[["layer_2_non_blood"]])[which(prediction_list[["layer_2_non_blood"]]$predicted_class %in%
                                                                                            c("Stromal Cell","Endothelial cell","Fibroblasts", "Myofibroblasts","Smooth Muscle Cells",
                                                                                              "Cancer Associated Fibroblasts","Cancer Associated Myofibroblasts",
                                                                                              "fibroblasts_score"))]
      if (length(non_stromal_predicted) > 0){
        print("Starting Layer 3 Non Stromal")
        prediction_list[["layer_3_non_stromal"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = non_stromal_predicted,
                                                                                                 layer = "layer_3_non_stromal", imputation = imputation, genes_in_model = top_genes_unlisted_layer_3_non_stromal,
                                                                                                 model = model_layer_3_non_stromal, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                 bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
        print("Done Layer 3 Non Stromal")
        non_GI_predicted <- row.names(prediction_list[["layer_3_non_stromal"]])[
          which(prediction_list[["layer_3_non_stromal"]]$predicted_class %in%c(
            "Breast Cancer",
            "Endometrial/Uterine Cancer", "Lung Cancer", "Ovarian Cancer",
            "Prostate Cancer", "Kidney Cancer",  "Non GI Epithelial Cell"))]
        soft_tissue_neuro_cancer_predicted <- row.names(prediction_list[["layer_3_non_stromal"]])[
          which(prediction_list[["layer_3_non_stromal"]]$predicted_class %in% c("Bone Cancer","Brain Cancer",
                                                                                "Neuroblastoma", "Skin Cancer",
                                                                                "Sarcoma", "Soft Tissue or Neuro Cancer Cell") )]
        GI_predicted <- row.names(prediction_list[["layer_3_non_stromal"]])[
          which(prediction_list[["layer_3_non_stromal"]]$predicted_class %in% c("Bile Duct Cancer","Bladder Cancer",
                                                                                "Colon/Colorectal Cancer","Esophageal Cancer",
                                                                                "Gallbladder Cancer", "Gastric Cancer", "Liver Cancer",
                                                                                "Pancreatic Cancer","GI Epithelial Cell") )]
        if (length(non_GI_predicted) > 0){
          print("Starting Layer 4 Non GI")
          prediction_list[["layer_4_non_GI"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = non_GI_predicted,
                                                                                              layer = "layer_4_non_GI", imputation = imputation, genes_in_model = top_genes_unlisted_layer_4_non_GI,
                                                                                              model = model_layer_4_non_GI, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                              bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
          print("Done Layer 4 Non GI")
          breast_lung_prostate_predicted <- row.names(prediction_list[["layer_4_non_GI"]] )[
            which(prediction_list[["layer_4_non_GI"]]$predicted_class  %in%c(
              "Breast Cancer","Lung Cancer",
              "Prostate Cancer",  "Breast/Lung/Prostate Cell"))]
          ov_endo_kid_predicted <- row.names(prediction_list[["layer_4_non_GI"]])[
            which(prediction_list[["layer_4_non_GI"]]$predicted_class %in%c(
              "Endometrial/Uterine Cancer", "Ovarian Cancer", "Kidney Cancer", "Ovarian/Endometrial/Kidney Cell"))]
          if (length(breast_lung_prostate_predicted) > 0){
            print("Starting Layer 5 Breast Lung Prostate")
            prediction_list[["layer_5_breast_lung_prostate"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = breast_lung_prostate_predicted,
                                                                                                              layer = "layer_5_breast_lung_prostate", imputation = imputation,
                                                                                                              genes_in_model = top_genes_unlisted_layer_5_breast_lung_prostate,
                                                                                                              model = model_layer_5_breast_lung_prostate,
                                                                                                              ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                              bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
            print("Done Layer 5 Breast Lung Prostate")
            if(breast_mode == T){
              breast_predicted <- row.names(prediction_list[["layer_5_breast_lung_prostate"]])[which(prediction_list[["layer_5_breast_lung_prostate"]]$predicted_tissue_with_cutoff == "Breast Cancer Cell")]
              if (length(breast_predicted) > 0){
                print("Starting Layer 6 Breast")
                prediction_list[["layer_6_breast"]] <- scATOMIC::classify_layer(rna_counts = rna_counts, cells_to_use = breast_predicted,
                                                                                          layer = "layer_6_breast", imputation = imputation,
                                                                                          genes_in_model = top_genes_unlisted_layer_6_breast,
                                                                                          model = model_layer_6_breast,
                                                                                          ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                          bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
                print("Done Layer 6 Breast")
              }
            }
          }
          if (length(ov_endo_kid_predicted) > 0){
            print("Starting Layer 5 Ovarian Endometrial Kidney")
            prediction_list[["layer_5_ov_endo_kid"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = ov_endo_kid_predicted,
                                                                                                     layer = "layer_5_ov_endo_kid", imputation = imputation,
                                                                                                     genes_in_model = top_genes_unlisted_layer_5_ov_endo_kid,
                                                                                                     model = model_layer_5_ov_endo_kid,
                                                                                                     ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                     bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
            print("Done Layer 5 Ovarian Endometrial Kidney")
          }
        }
        if (length(GI_predicted) > 0){
          print("Starting Layer 4 GI")
          prediction_list[["layer_4_GI"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = GI_predicted,
                                                                                          layer = "layer_4_GI", imputation = imputation, genes_in_model = top_genes_unlisted_layer_4_GI,
                                                                                          model = model_layer_4_GI, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                          bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
          print("Done Layer 4 GI")
          digestive_predicted <- row.names(prediction_list[["layer_4_GI"]])[
            which(prediction_list[["layer_4_GI"]]$predicted_class %in% c("Colon/Colorectal Cancer","Esophageal Cancer",
                                                                         "Gastric Cancer","Colorectal/Esophageal/Gastric Cell"))]
          biliary_predicted <- row.names(prediction_list[["layer_4_GI"]])[
            which(prediction_list[["layer_4_GI"]]$predicted_class %in% c("Bile Duct Cancer","Bladder Cancer",
                                                                         "Gallbladder Cancer",  "Liver Cancer",
                                                                         "Pancreatic Cancer", "Billiary Cell"))]

          if (length(digestive_predicted) > 0){
            print("Starting Layer 5 Digestive")
            prediction_list[["layer_5_digestive"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = digestive_predicted,
                                                                                                   layer = "layer_5_digestive", imputation = imputation,
                                                                                                   genes_in_model = top_genes_unlisted_layer_5_digestive,
                                                                                                   model = model_layer_5_digestive,
                                                                                                   ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                   bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
            print("Done Layer 5 Digestive")
          }
          if (length(biliary_predicted) > 0){
            print("Starting Layer 5 Biliary")
            prediction_list[["layer_5_biliary"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = biliary_predicted,
                                                                                                 layer = "layer_5_biliary", imputation = imputation,
                                                                                                 genes_in_model = top_genes_unlisted_layer_5_biliary,
                                                                                                 model = model_layer_5_biliary,
                                                                                                 ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                 bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
            print("Done Layer 5 Biliary")
          }
        }
        if (length(soft_tissue_neuro_cancer_predicted) > 0){
          print("Starting Layer 4 Soft Tissue Neuro")
          prediction_list[["layer_4_soft_tissue_neuro"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = soft_tissue_neuro_cancer_predicted,
                                                                                                         layer = "layer_4_soft_tissue_neuro", imputation = imputation,
                                                                                                         genes_in_model = top_genes_unlisted_layer_4_soft_neuro_cancer,
                                                                                                         model = model_layer_4_soft_neuro_cancer,
                                                                                                         ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                         bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
          print("Done Layer 4 Soft Tissue Neuro")
          soft_tissue_neuro_cancer_no_lung_skin_predicted <- row.names(prediction_list[["layer_4_soft_tissue_neuro"]])[
            which(prediction_list[["layer_4_soft_tissue_neuro"]]$predicted_class %in% c("Bone Cancer","Brain Cancer",
                                                                                        "Neuroblastoma",
                                                                                        "Sarcoma", "Soft Tissue or Neuro Cancer Cell", "Soft Tissue or Neuro Cancer Cell"))]
          if (length(soft_tissue_neuro_cancer_no_lung_skin_predicted) > 0){
            print("Starting Layer 5 Soft Tissue Neuro")
            prediction_list[["layer_5_soft_tissue_neuro"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = soft_tissue_neuro_cancer_no_lung_skin_predicted,
                                                                                                           layer = "layer_5_soft_tissue_neuro", imputation = imputation,
                                                                                                           genes_in_model = top_genes_unlisted_layer_5_soft_neuro_cancer_no_lung_skin,
                                                                                                           model = model_layer_5_soft_neuro_cancer_no_lung_skin,
                                                                                                           ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                           bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
            print("Done Layer 5 Soft Tissue Neuro")
            brain_nbm_predicted <- row.names(prediction_list[["layer_5_soft_tissue_neuro"]])[
              which(prediction_list[["layer_5_soft_tissue_neuro"]]$predicted_class %in% c("Brain Cancer",
                                                                                          "Neuroblastoma", "Brain/Neuroblastoma Cancer Cell"))]
            soft_tissue_predicted <- row.names(prediction_list[["layer_5_soft_tissue_neuro"]])[which(prediction_list[["layer_5_soft_tissue_neuro"]]$predicted_class %in% c("Bone Cancer",
                                                                                                                                                                           "Sarcoma", "Soft Tissue Cancer Cell"))]
            if (length(brain_nbm_predicted) > 0){
              print("Starting Layer 6 Brain NBM")
              prediction_list[["layer_6_brain_nbm"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = brain_nbm_predicted,
                                                                                                     layer = "layer_6_brain_nbm", imputation = imputation,
                                                                                                     genes_in_model = top_genes_unlisted_layer_6_brain_nbm,
                                                                                                     model = model_layer_6_brain_nbm,
                                                                                                     ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                     bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
              print("Done Layer 6 Brain NBM")
            }
            if (length(soft_tissue_predicted) > 0){
              print("Starting Layer 6 Soft Tissue")
              prediction_list[["layer_6_soft_tissue"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = soft_tissue_predicted,
                                                                                                       layer = "layer_6_soft_tissue", imputation = imputation,
                                                                                                       genes_in_model = top_genes_unlisted_layer_6_soft_tissue,
                                                                                                       model = model_layer_6_soft_tissue,
                                                                                                       ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                       bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
              print("Done Layer 6 Soft Tissue")
            }
          }
        }
      }
      if (length(normal_stromal_predicted) > 0){
        print("Starting Layer 3 Normal Stromal")
        prediction_list[["layer_3_stromal"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = normal_stromal_predicted,
                                                                                             layer = "layer_3_stromal", imputation = imputation, genes_in_model = top_genes_unlisted_layer_3_stromal,
                                                                                             model = model_layer_3_stromal, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                             bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)
        print("Done Layer 3 Normal Stromal")
        CAF_predicted <- row.names(prediction_list[["layer_3_stromal"]])[
          which(prediction_list[["layer_3_stromal"]]$predicted_class == "Cancer Associated Fibroblasts")]
        if (length(CAF_predicted) > 0){
          print("Starting Layer 4 CAFs")
          prediction_list[["layer_4_CAF"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = CAF_predicted,
                                                                           layer = "layer_4_CAF", imputation = F,
                                                                           genes_in_model = top_genes_unlisted_layer_4_CAFs,
                                                                           model = model_layer_4_CAFs,
                                                                           ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                           bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

          print("Done Layer 4 CAFs")
        }
      }
    }
    if (length(blood_cells_predicted) > 0){
      print("Starting Layer 2 Blood")
      prediction_list[["layer_2_blood"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = blood_cells_predicted,
                                                                                         layer = "layer_2_blood", imputation = imputation, genes_in_model = top_genes_unlisted_layer_2_blood,
                                                                                         model = model_layer_2_blood, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                         bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

      print("Done Layer 2 Blood")
      TNK_cells_predicted <- row.names(prediction_list[["layer_2_blood"]])[
        which(prediction_list[["layer_2_blood"]]$predicted_class %in% c("CD4+ T cell", "CD8+ T cell", "Natural killer cell", "T_or_NK_lymphocyte"))]
      MDC_cells_predicted <- row.names(prediction_list[["layer_2_blood"]])[
        which(prediction_list[["layer_2_blood"]]$predicted_class %in% c("macrophage_or_dendritic_cell", "macrophage_DC_score"))]
      B_cells_predicted <- row.names(prediction_list[["layer_2_blood"]])[
        which(prediction_list[["layer_2_blood"]]$predicted_class %in% c("B cell","Plasmablast" ,"B_cell_score", "B cell or Plasmablast"))]
      if (length(TNK_cells_predicted) > 0){

        print("Starting Layer 3 TNK")
        prediction_list[["layer_3_TNK"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = TNK_cells_predicted,
                                                                                         layer = "layer_3_TNK", imputation = F, genes_in_model = top_genes_unlisted_layer_3_TNK,
                                                                                         model = model_layer_3_TNK, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                         bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

        print("Done Layer 3 TNK")
        CD4_CD8_cells_predicted <- row.names(prediction_list[["layer_3_TNK"]])[
          which(prediction_list[["layer_3_TNK"]]$predicted_class %in% c("CD4+ T cell","CD8+ T cell" ,"CD4 or CD8 T cell"))]
        CD8_NK_cells_predicted <- row.names(prediction_list[["layer_3_TNK"]])[
          which(prediction_list[["layer_3_TNK"]]$predicted_class %in% c("Natural killer cell", "NK or CD8 T cell"))]

        if (length(CD4_CD8_cells_predicted) > 0){
          print("Starting Layer 4 CD4 CD8")
          prediction_list[["layer_4_CD4_CD8"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = CD4_CD8_cells_predicted,
                                                                                               layer = "layer_4_CD4_CD8", imputation = F,
                                                                                               genes_in_model = top_genes_unlisted_layer_4_CD4_CD8,
                                                                                               model = model_layer_4_CD4_CD8,
                                                                                               ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                               bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

          print("Done Layer 4 CD4 CD8")
        }
        if (length(CD8_NK_cells_predicted) > 0){
          print("Starting Layer 4 CD8 NK")
          prediction_list[["layer_4_CD8_NK"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = CD8_NK_cells_predicted,
                                                                                              layer = "layer_4_CD8_NK", imputation = F,
                                                                                              genes_in_model = top_genes_unlisted_layer_4_CD8_NK,
                                                                                              model = model_layer_4_CD8_NK,
                                                                                              ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                              bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

          print("Done Layer 4 CD8 NK")
        }
        if(fine_grained_T == T){
          CD4_cells_predicted <- row.names(prediction_list[["layer_4_CD4_CD8"]])[
            which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD4+ T cell"))]
          CD8_cells_predicted <- c(row.names(prediction_list[["layer_4_CD4_CD8"]])[
            which(prediction_list[["layer_4_CD4_CD8"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell"))],row.names(prediction_list[["layer_4_CD8_NK"]])[
              which(prediction_list[["layer_4_CD8_NK"]]$predicted_tissue_with_cutoff %in% c("CD8+ T cell"))] )
          if (length(CD4_cells_predicted) > 0){
            print("Starting Layer 5 CD4")
            prediction_list[["layer_5_CD4"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = CD4_cells_predicted,
                                                                                             layer = "layer_5_CD4", imputation = F,
                                                                                             genes_in_model = top_genes_unlisted_layer_5_CD4,
                                                                                             model = model_layer_5_CD4,
                                                                                             ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                             bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

            print("Done Layer 5 CD4")
          }
          if (length(CD8_cells_predicted) > 0){
            print("Starting Layer 5 CD8")
            prediction_list[["layer_5_CD8"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = CD8_cells_predicted,
                                                                                             layer = "layer_5_CD8", imputation = F,
                                                                                             genes_in_model = top_genes_unlisted_layer_5_CD8,
                                                                                             model = model_layer_5_CD8,
                                                                                             ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                             bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

            print("Done Layer 5 CD8")
          }






        }




      }
      if (length(MDC_cells_predicted) > 0){
        print("Starting Layer 3 Myeloid")
        prediction_list[["layer_3_myeloid"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = MDC_cells_predicted,
                                                                                             layer = "layer_3_MDC", imputation = imputation, genes_in_model = top_genes_unlisted_layer_3_MDC,
                                                                                             model = model_layer_3_MDC, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                             bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

        print("Done Layer 3 Myeloid")
        DC_cells_predicted <- row.names(prediction_list[["layer_3_myeloid"]])[
          which(prediction_list[["layer_3_myeloid"]]$predicted_class == "Dendritic Cell")]
        macrophage_monocyte_predicted <- row.names(prediction_list[["layer_3_myeloid"]])[
          which(prediction_list[["layer_3_myeloid"]]$predicted_class == "Macrophage or Monocyte")]
        if (length(DC_cells_predicted) > 0){
          print("Starting Layer 4 Dendritic")
          prediction_list[["layer_4_dendritic"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = DC_cells_predicted,
                                                                                                 layer = "layer_4_dendritic", imputation = imputation,
                                                                                                 genes_in_model = top_genes_unlisted_layer_4_DC_cell,
                                                                                                 model = model_layer_4_DC_cell,
                                                                                                 ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                 bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

          print("Done Layer 4 Dendritic")
          cDC_predicted <- row.names(prediction_list[["layer_4_dendritic"]])[
            which(prediction_list[["layer_4_dendritic"]]$predicted_class == "cDC")]
          if (length(cDC_predicted) > 0){
            print("Starting Layer 5 cDC")
            prediction_list[["layer_5_cDC"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = cDC_predicted,
                                                                             layer = "layer_5_cDC", imputation = F,
                                                                             genes_in_model = top_genes_unlisted_layer_5_cDC,
                                                                             model = model_layer_5_cDC,
                                                                             ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                             bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

            print("Done Layer 5 cDC")
          }
        }
        if (length(macrophage_monocyte_predicted) > 0){
          print("Starting Layer 4 Macrophage Monocyte")
          prediction_list[["layer_4_macrophage"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = macrophage_monocyte_predicted,
                                                                                                  layer = "layer_4_macrophage", imputation = imputation,
                                                                                                  genes_in_model = top_genes_unlisted_layer_4_macrophage_monocyte,
                                                                                                  model = model_layer_4_macrophage_monocyte,
                                                                                                  ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                                  bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

          print("Done Layer 4 Macrophage Monocyte")
          macrophage_predicted <- row.names(prediction_list[["layer_4_macrophage"]])[
            which(prediction_list[["layer_4_macrophage"]]$predicted_class == "Macrophage")]
          monocyte_predicted <- row.names(prediction_list[["layer_4_macrophage"]])[
            which(prediction_list[["layer_4_macrophage"]]$predicted_class == "Monocyte")]
          if (length(macrophage_predicted) > 0){
            print("Starting Layer 5 Macrophage")
            prediction_list[["layer_5_macrophage"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = macrophage_predicted,
                                                                                    layer = "layer_5_macrophage", imputation = F,
                                                                                    genes_in_model = top_genes_unlisted_layer_5_macrophage,
                                                                                    model = model_layer_5_macrophage,
                                                                                    ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                    bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

            print("Done Layer 5 Macrophage")
          }
          if (length(monocyte_predicted) > 0){
            print("Starting Layer 5 Monocyte")
            prediction_list[["layer_5_monocyte"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = monocyte_predicted,
                                                                                  layer = "layer_5_monocyte", imputation = F,
                                                                                  genes_in_model = top_genes_unlisted_layer_5_monocyte,
                                                                                  model = model_layer_5_monocyte,
                                                                                  ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                  bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

            print("Done Layer 5 Monocyte")
          }
        }

      }
      if (length(B_cells_predicted) > 0){
        print("Starting Layer 3 B Cell")
        prediction_list[["layer_3_BCell"]] <- scATOMIC::classify_layer_no_cutoff(rna_counts = rna_counts, cells_to_use = B_cells_predicted,
                                                                                           layer = "layer_3_BCell", imputation = imputation, genes_in_model = top_genes_unlisted_layer_3_BCell,
                                                                                           model = model_layer_3_BCell, ref_based = ref_based, mc.cores = mc.cores, unimodal_nsd = unimodal_nsd,
                                                                                           bimodal_nsd = bimodal_nsd, normalized_counts =normalized_counts)

        print("Done Layer 3 B Cell")
      }
    }
    return(prediction_list)

  }

}
