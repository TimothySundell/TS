#' Process CellRanger output files with Seurat
#'
#' @description
#' Input the three matrices from CellRanger and analyse them automatically.
#' Based on Seurat, using the workflow published in Sundell et. al., 2022 (https://doi.org/10.1093/bfgp/elac044)
#'
#' @param input_data File path to the folder containing the three matrices used for Seurat
#' @param sample_ID Required ID-tag for your sample. Will be added as metadata in column "sample_ID". Handy in downstream analyses, e.g. integration, batch-correction etc.
#' @param remove_ig_genes Filters out Ig-genes prior to clustering. Which can help in finding biologically relevant clusters, as shown in Sundell et. al., 2022. Defaults to TRUE.
#' @param min.cells Include features/genes detected in at least this many cells. Defaults to 3.
#' @param min.features Keep cells with at least this many features detected. Defaults to 200.
#' @param project_name Required by Seurat. Defaults to "Seurat". Added automatically as metadata in column "orig.ident"
#' @param normalization.method Normalization method to use. Alternatives are "LogNormalize" (default), "CLR", "RC".
#' @param scale.factor Scale factor for cell-level normalization. Defaults to 10000
#' @param selection.method Method for finding variable features. Alternatives are "vst" (default), "mean.var.plot", "dispersion".
#' @param nFeatures Number of variable features to calculate. Defaults to 2000
#' @param regular_clustering_dims What PCA components to use for building NN graph and dimensionality reduction with regular normalization/scaling. Defaults to "1:10".
#' @param regular_clustering_nn.method Nearest-neighbour method to use for clustering with regular normalization/scaling. Alternatives are "rann" or "annoy". Defaults to "rann"
#' @param SCT_clustering_dims What PCA components to use for building NN graph and dimensionality reduction with SCTransformed data. Defaults to "1:20".
#'
#' @export

TS_Seurat_from_file <- function(
    input_data,
    sample_ID = NULL,
    remove_ig_genes = T,
    regular_clustering_dims = 1:10,
    SCT_clustering_dims = 1:20,
    project_name = NULL,
    min.cells = 3,
    min.features = 200,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    selection.method = "vst",
    nfeatures = 2000,
    regular_clustering_nn.method = "rann"){

  require(Seurat)
  require(patchwork)
  require(clustree)

  if(is.null(sample_ID)){
    stop("\nYou need to specify sample_ID\nSample ID will be added as metadata to your output Seurat object\n\n")
  }

  if(is.null(project_name)){
    project_name <- sample_ID
  }

  input.data <- Seurat::Read10X(data.dir = input_data)
  seurat_object <- Seurat::CreateSeuratObject(counts = input.data, project = project_name, min.cells = min.cells, min.features = min.features)

  cat("Calculating proportion mitochondrial transcripts\n")
  seurat_object[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^MT-")
  cat("Calculating proportion ribosomal transcripts\n")
  seurat_object[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")
  cat("Calculating proportion of MALAT1 transcripts\n")
  seurat_object[["percent.malat1"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "MALAT1")

  seurat_object <- Seurat::AddMetaData(object = seurat_object, metadata = sample_ID, col.name = "sample_ID")

  if(remove_ig_genes){
    # Removing Ig-genes
    cat("Removing Ig-genes prior to clustering\n")
    seurat_object2 <- seurat_object[!grepl("^IG[HKL][VJCADEGM]", rownames(seurat_object)), ]
  }
  else if(!remove_ig_genes){
    seurat_object2 <- seurat_object
  }

  # Normalize and scale data
  cat("Normalizing and scaling data\n")
  seurat_object2 <- Seurat::NormalizeData(seurat_object2, normalization.method = normalization.method, scale.factor = scale.factor, verbose = F)
  seurat_object2 <- Seurat::FindVariableFeatures(seurat_object2, selection.method = selection.method, nfeatures = nfeatures, verbose = F)
  seurat_object2 <- Seurat::ScaleData(seurat_object2, features = rownames(seurat_object2), verbose = F)
  # Run PCA and make NN graph
  cat("Running PCA\n")
  seurat_object2 <- Seurat::RunPCA(seurat_object2, features = Seurat::VariableFeatures(object = seurat_object2), verbose = F)
  cat("Building NN graph\n")
  seurat_object2 <- Seurat::FindNeighbors(seurat_object2, dims = regular_clustering_dims, nn.method = regular_clustering_nn.method, verbose = F)

  # Set resolution parameter
  cat("Finding those clusters...\n")
  seurat_object2 <- Seurat::FindClusters(seurat_object2, resolution = seq(0, 2, 0.2), verbose = F)
  print(clustree::clustree(seurat_object2))
  chosen_resolution <- readline(prompt = "Choose resolution for clustering: ")
  seurat_object2 <- Seurat::FindClusters(seurat_object2, resolution = as.numeric(chosen_resolution), verbose = F)
  cat("Calculating UMAP coordinates\n")
  seurat_object2 <- Seurat::RunUMAP(seurat_object2, dims = regular_clustering_dims, verbose = F)

  if(remove_ig_genes){
    # Add back Ig-genes
    cat("Adding back Ig-genes in assay 'Iggenes'\n")
    seurat_object2[["Iggenes"]] <- Seurat::CreateAssayObject(counts = Seurat::GetAssayData(object = seurat_object, slot = "count", assay = "RNA"))
    seurat_object2 <- Seurat::NormalizeData(seurat_object2, normalization.method = normalization.method, scale.factor = scale.factor, assay = "Iggenes", verbose = F)
    seurat_object2 <- Seurat::ScaleData(seurat_object2, features = rownames(seurat_object2), assay = "Iggenes", verbose = F)
    Seurat::DefaultAssay(seurat_object2) <- "RNA"
  }
  # Plot the stuff
  plot_list <- list()
  plot_list[["Dimplot"]] <- Seurat::DimPlot(object = seurat_object2, label = T)
  plot_list[["Vlnplots"]] <- suppressWarnings(Seurat::VlnPlot(object = seurat_object2, features = c("percent.mt", "percent.ribo", "percent.malat1", "CD19", "CD3E", "CD14", "CD56"), ncol = 1))
  print(patchwork::wrap_plots(plot_list, ncol = 2))

  # Subset wanted clusters
  keep_clusters <- readline('What clusters do you want to keep? Numbers separated only by commas: ')
  cat("Subsetting and running normalisation and scaling of the new data\n")
  keep_clusters <- strsplit(keep_clusters, ',')
  keep_clusters <- as.numeric(unlist(keep_clusters))
  seurat_object2 <- base::subset(seurat_object2, idents = keep_clusters)

  DefaultAssay(seurat_object2) <- "RNA"
  seurat_object2 <- Seurat::NormalizeData(object = seurat_object2, scale.factor = scale.factor, verbose = F)
  seurat_object2 <- Seurat::ScaleData(object = seurat_object2, features = rownames(seurat_object2), verbose = F)

  run_SCT <- as.character(readline("Do you want to run SCTransform as well? Enter 'y' or 'n': "))

  if(run_SCT == "y"){
    cat("Running SCTransform\n")
    seurat_object2 <- Seurat::SCTransform(object = seurat_object2, verbose = FALSE, ncells = NULL)
    cat("Running PCA\n")
    seurat_object2 <- Seurat::RunPCA(object = seurat_object2, verbose = FALSE)
    cat("Calculating shiny new UMAP coordinates\n")
    seurat_object2 <- Seurat::RunUMAP(object = seurat_object2, dims = SCT_clustering_dims, verbose = FALSE)
    seurat_object2 <- Seurat::FindNeighbors(object = seurat_object2, dims = SCT_clustering_dims, verbose = FALSE, nn.method = regular_clustering_nn.method)
    cat(paste0("Finding clusters with resolution ", as.numeric(chosen_resolution), "\n"))
    seurat_object2 <- Seurat::FindClusters(object = seurat_object2, verbose = FALSE, resolution = as.numeric(chosen_resolution))
    print(Seurat::DimPlot(object = seurat_object2, label = T))
  }

  else if(run_SCT == "n"){
    cat("Normalizing and scaling data with filtered data\n")
    seurat_object2 <- Seurat::NormalizeData(seurat_object2, normalization.method = normalization.method, scale.factor = scale.factor, verbose = F)
    seurat_object2 <- Seurat::FindVariableFeatures(seurat_object2, selection.method = selection.method, nfeatures = nfeatures, verbose = F)
    seurat_object2 <- Seurat::ScaleData(seurat_object2, features = rownames(seurat_object2), verbose = F)
    # Run PCA and make NN graph
    cat("Running PCA and building NN graph with filtered data\n")
    seurat_object2 <- Seurat::RunPCA(seurat_object2, features = Seurat::VariableFeatures(object = seurat_object2), verbose = F)
    seurat_object2 <- Seurat::FindNeighbors(seurat_object2, dims = regular_clustering_dims, nn.method = regular_clustering_nn.method, verbose = F)

    # Set resolution parameter
    cat("Finding clusters with filtered data\n")
    seurat_object2 <- Seurat::FindClusters(seurat_object2, resolution = as.numeric(chosen_resolution), verbose = F)
    cat("Calculating UMAP coordinates with filtered data\n")
    seurat_object2 <- Seurat::RunUMAP(seurat_object2, dims = regular_clustering_dims, verbose = F)
    print(Seurat::DimPlot(object = seurat_object2, label = T))
    cat("\nExporting filtered Seurat object\n")
  }



  cat("Done")
  return(seurat_object2)
}
