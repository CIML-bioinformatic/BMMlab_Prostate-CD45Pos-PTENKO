###############################################################################
# This file defines SAMPLE parameters as global variables that will be loaded
# before analysis starts. 
#

SAMPLE_NAME_CTRL3mo_1 = "CTRL3mo_1"
SAMPLE_NAME_CTRL3mo_2 = "CTRL3mo_2"
SAMPLE_NAME_PTEN3mo_1 = "PTEN3mo_1"
SAMPLE_NAME_PTEN3mo_2 = "PTEN3mo_2"
SAMPLE_NAME_LIST = c( SAMPLE_NAME_CTRL3mo_1, SAMPLE_NAME_CTRL3mo_2, SAMPLE_NAME_PTEN3mo_1, SAMPLE_NAME_PTEN3mo_2)

SAMPLE_COLOR = c( "#00BA38", "#00661f", "#BA0038", "#66001f")
names( SAMPLE_COLOR) = SAMPLE_NAME_LIST

CONDITION_NAME_CTRL = "CTRL"
CONDITION_NAME_PTEN = "PTEN"
CONDITION_NAME_LIST = c( CONDITION_NAME_CTRL, CONDITION_NAME_PTEN)

CONDITION_COLOR = c( "darkolivegreen3", "brown2")
names( CONDITION_COLOR) = CONDITION_NAME_LIST

SAMPLE_CONDITION = c( CONDITION_NAME_CTRL, CONDITION_NAME_CTRL, CONDITION_NAME_PTEN, CONDITION_NAME_PTEN)
names( SAMPLE_CONDITION) = SAMPLE_NAME_LIST

INTEGRATION_SAMPLE_NAME = paste0( "Merge", SAMPLE_NAME_CTRL3mo_1, SAMPLE_NAME_CTRL3mo_2, SAMPLE_NAME_PTEN3mo_1, SAMPLE_NAME_PTEN3mo_2, collapse = "_")
