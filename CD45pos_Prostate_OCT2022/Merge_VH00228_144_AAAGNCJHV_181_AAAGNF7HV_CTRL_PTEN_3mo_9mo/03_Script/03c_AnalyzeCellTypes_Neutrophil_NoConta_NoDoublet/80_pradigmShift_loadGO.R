
library( GOfuncR)
library( ontologyIndex)
library( Mus.musculus)

STATUS_NA = NA
STATUS_TO_ANALYZE = "To analyze"
STATUS_SELECTED= "Selected"
STATUS_TOO_FEW_GENES = "Too few genes"
STATUS_TOO_MANY_GENES = "Too many genes"
STATUS_NO_GENE_FOUND = "No gene found"

# Load the GO graph
go_onto = get_ontology( file = "/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/06_CB2M/01_DOCUMENTATION/01_PROJECT/MoonGOv2/MoonGO/07_Output/Homo_Sapiens/09_AnnotatePartition/go-basic.obo",
                        extract_tags = "everything")
go_onto_df = as.data.frame( go_onto)

# Select only biological process non obsolte terms
bp_go_onto_df = go_onto_df[ which( go_onto_df$namespace == "biological_process"), ]
bp_go_onto_df = go_onto_df[ which( go_onto_df$obsolete == FALSE), ]

# Initialize new colums to report the genes associated to terms
bp_go_onto_df$gene.number = NA
bp_go_onto_df$gene.list = NA
bp_go_onto_df$status = STATUS_NA
bp_go_onto_df$depth = NA

# Initialize the search of GO Terms deapest terms with enough associated genes
leaf_go_term = bp_go_onto_df[ which( bp_go_onto_df$children == ""), "id"]
bp_go_onto_df[ leaf_go_term, "status"] = STATUS_TO_ANALYZE

# Parse the list of term to analyse, check if they have enough associated genes
# If so, select them and note the associated genes list and number
# If not, try to analyse its parent terms are more suitable
depth = 0
repeat{
  
  depth = depth + 1
  
  # get the list of term to analyse this round
  toanalyse_go_term = bp_go_onto_df[ which( bp_go_onto_df$status == STATUS_TO_ANALYZE), "id"]
  
  cat("\nNUMBER OF TERMS TO ANALYZE:", length( toanalyse_go_term))
  cat("\nNUMBER OF TERMS WITH TOO FEW GENES:", length( which( bp_go_onto_df$status == STATUS_TOO_FEW_GENES)))
  cat("\nNUMBER OF TERMS WITH NO GENE FOUND:", length( which( bp_go_onto_df$status == STATUS_NO_GENE_FOUND)))
  cat("\nNUMBER OF TERMS SELECTED:", length( which( bp_go_onto_df$status == STATUS_SELECTED)))
  
  # If there is no term to analyse, stop the process
  if(length( toanalyse_go_term) == 0 || depth > 3 ){
    break
  }

  tryCatch(
    { 
      # Get the list of genes associated to the list of terms to analyse
      gene_df = get_anno_genes( go_ids=toanalyse_go_term, database = "Mus.musculus");
      
      # Initialize the status of the go term to analyse
      bp_go_onto_df[ toanalyse_go_term, "status"] = STATUS_NO_GENE_FOUND
      
      # Parse the GO terms for which associated genes has been found
      for( go_term in unique( gene_df$go_id)){
        bp_go_onto_df[ go_term, "depth"] = depth
        # Get the list of genes associated to the current GO term
        bp_go_onto_df[ go_term, "gene.number"] = length( which( gene_df$go_id == go_term))
        # Store the list of genes in the dataframe
        bp_go_onto_df[ go_term, "gene.list"] = paste( gene_df[ which( gene_df$go_id == go_term), "gene"], collapse = ";")
        # If there is enough genes, select the GO term
        # If not, set the status of the GO term parents to "to analyse" in order to have them analyzed in the next round
        if( bp_go_onto_df[ go_term, "gene.number"] > 200){
          bp_go_onto_df[ go_term, "status"] = STATUS_TOO_MANY_GENES
        }
        else if( bp_go_onto_df[ go_term, "gene.number"] < 20){
          bp_go_onto_df[ go_term, "status"] = STATUS_TOO_FEW_GENES
          parent_term_list = str_split( bp_go_onto_df[ go_term, "parents"], ";")[[1]]
          for( parent_term in parent_term_list){
            if( parent_term %in% bp_go_onto_df$id 
                && bp_go_onto_df[ parent_term, "namespace"] == "biological_process"
                && is.na( bp_go_onto_df[ parent_term, "status"])){
              bp_go_onto_df[ parent_term, "status"] = STATUS_TO_ANALYZE
            }
          }
        }else{
          bp_go_onto_df[ go_term, "status"] = STATUS_SELECTED
        }
      }
      
    }
    , 
    # If the search for genes associated to the list of GO Term to analyze, returns an error, it means no genes
    # were found and the GO terms must no be selected
    error = function(e){ bp_go_onto_df[ toanalyse_go_term, "status"] = STATUS_NO_GENE_FOUND}
  )
}

selected_bp_go_onto_df = bp_go_onto_df[ which( bp_go_onto_df$status == STATUS_SELECTED), ]
summary( selected_bp_go_onto_df$gene.number)
bp_go_onto_df$id[ which( bp_go_onto_df$gene.number == 17686)]
