ExecuteCellchat = function(merge_2_seurat_list,output){
  for(i in 1:length(merge_2_seurat_list)){
    cts = names(merge_2_seurat_list[i]) # cts:cancer_type_site
    cts = gsub("/","_",cts)
    message("cts is ",cts)
    merge_2_seurat = merge_2_seurat_list[[i]]
    # For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
    data.input <- merge_2_seurat[["RNA"]]$data # normalized data matrix
    labels <- Idents(merge_2_seurat)
    meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
    # Create a CellChat object
    cellchat <- createCellChat(object = merge_2_seurat,
                               # group.by = "celltype_manual_harmony_recluster_1_coarse", 
                               assay = "RNA")
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
    CellChatDB.use <- CellChatDB
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    cellchat <- identifyOverExpressedGenes(cellchat);message("identifyOverExpressedGenes was done")
    cellchat <- identifyOverExpressedInteractions(cellchat);message("identifyOverExpressedInteractions was done")
    cellchat <- projectData(cellchat, PPI.human)
    #> The number of highly variable ligand-receptor pairs used for signaling inference is 692
    # invisible(capture.output(cellchat <- computeCommunProb(cellchat, population.size = FALSE, raw.use = TRUE), type = "message"))
    invisible(capture.output(cellchat <- computeCommunProb(cellchat, population.size = FALSE, raw.use = TRUE)))
    cellchat <- filterCommunication(cellchat)
    cellchat <- computeCommunProbPathway(cellchat)
    # Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
    saveRDS(cellchat,file = paste0("./",output,"/cellchat_",cts,".rds"))
    message("The cellchat of ",cts," has done.")
  }
}
