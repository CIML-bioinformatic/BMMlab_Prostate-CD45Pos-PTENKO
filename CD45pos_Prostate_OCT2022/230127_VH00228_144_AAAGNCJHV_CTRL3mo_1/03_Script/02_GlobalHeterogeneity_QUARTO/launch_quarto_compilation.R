# ############################################################################
#  This script aims to launch the Quarto rendering of the R markdown report
#  using the parameters files and ensuring the HTML report will be stored
#  in the right location (analysis output folder)
# ############################################################################

# ...................................................................
# Define working Directory
# ...................................................................

library( funr)
if( exists( "snakemake")){
  cat("\nWorking in Snakemake mode")
  WORKING_DIR = snakemake@scriptdir
}else{
  cat("\nWorking in local script mode")
  WORKING_DIR = tryCatch(
    {
      dirname( sys.script())
    },
    error=function(cond) {
      cat("\nWorking in local R session mode")
      return( getwd())
    }
  )
}

# ...................................................................
# Load the analysis params
# ...................................................................

source( file.path( WORKING_DIR, "load_params.R"))

cat("\nINITIAL WORKING DIR=", WORKING_DIR)
cat("\nINITIAL PATH_ANALYSIS_OUTPUT=", PATH_ANALYSIS_OUTPUT)

# ...................................................................
# Launch the report compilation
# ...................................................................

# Define the name and path of the HTML report file to be produced
report_file_name = paste0( SCIENTIFIC_PROJECT_NAME, "_", EXPERIMENT_NAME, "_", ANALYSIS_STEP_NAME, "_Sample", SAMPLE_NAME, "_QUARTO.html")
report_file_path = file.path( WORKING_DIR, report_file_name)

# Render the HTML report using Quarto
library(quarto)
quarto_render( input = file.path( WORKING_DIR, "Report.qmd"),
               output_file  = report_file_name,
               execute_dir = WORKING_DIR,
               output_format = "html")

# ...................................................................
# Try to copy the HTML report file to the analysis output folder
# and to remove it from the analysis script folder
# ...................................................................

# Try to copy the file to the analysis output folder
cat("\n  Copying HTML report file to", PATH_ANALYSIS_OUTPUT)
copy_done = file.copy(from = report_file_path,
                      to   = PATH_ANALYSIS_OUTPUT,
                      overwrite = TRUE)

# Analyse is copy was performed. 
# If so, try to remove the report file from the script folder.
# If copy or removal was not performed, warn the user

if( copy_done){
  cat("\n|-- Copy done")
  cat("\nRemoving HTML file from script folder")
  remove_done = file.remove( report_file_path)
  if( ! remove_done){
    cat("\n\-- ERROR: Report file was not removed from script folder")
  }else{
    cat("\n|-- Removal done")  
  }
}else{
  cat("\n|-- ERROR: Report file was not copied to analysis output folder")
}

