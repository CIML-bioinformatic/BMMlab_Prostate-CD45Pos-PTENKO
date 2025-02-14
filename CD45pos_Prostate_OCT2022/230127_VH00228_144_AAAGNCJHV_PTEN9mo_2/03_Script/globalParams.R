###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
#


#### General

GLOBAL_DESCRIPTION = "Prostate project"

SCIENTIFIC_GROUP = "HETMLAB"
SCIENTIFIC_PROJECT_NAME = "Prostate_Project_BMSH_DP"
SCIENTIFIC_SUBPROJECT_NAME = "CD45pos_Prostate_OCT2022"
EXPERIMENT_NAME = "230127_VH00228_144_AAAGNCJHV_PTEN9mo_2"

#### Additional custom variables (for a specific analysis step, tools, ...)

#SAMPLE_ID = ""



#### Input / Output

# Output folder name in data folder (for R session object, lists of cells/genes) 
PATH_PROJECT = file.path( "/mnt/DOSI", 
                                        SCIENTIFIC_GROUP,
                                        "BIOINFO", 
                                        "Project",
                                        SCIENTIFIC_PROJECT_NAME,
                                        SCIENTIFIC_SUBPROJECT_NAME)

PATH_EXPERIMENT = file.path( PATH_PROJECT, EXPERIMENT_NAME)

PATH_EXPERIMENT_RAWDATA       = file.path( PATH_EXPERIMENT, "00_RawData")
PATH_EXPERIMENT_REFERENCE     = file.path( PATH_EXPERIMENT, "01_Reference")
PATH_EXPERIMENT_OUTPUT        = file.path( PATH_EXPERIMENT, "05_Output")



# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_", SCIENTIFIC_SUBPROJECT_NAME, "_",
                            EXPERIMENT_NAME, "_"
                            )



#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



