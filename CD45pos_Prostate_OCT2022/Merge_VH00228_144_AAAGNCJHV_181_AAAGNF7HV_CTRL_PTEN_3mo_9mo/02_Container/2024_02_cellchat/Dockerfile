FROM rocker/tidyverse:4.2

MAINTAINER Lionel Spinelli (lionel.spinelli@univ-amu.fr)



# ##############################
# R packages system dependencies
# ##############################

RUN apt-get update && apt-get install -y optipng \
                                         pandoc \
                                         libpng-dev \
                                         zlib1g-dev \
                                         libssl-dev \
                                         libcurl4-openssl-dev

# Add Leiden Algorithm for clustering
RUN apt update && apt install -y python3-pip
RUN pip3 install leidenalg

# ##########
# R PACKAGES 
# ##########

#### BiocManager
RUN Rscript -e 'install.packages( "BiocManager")'


#### RMarkdown & Reporting
RUN Rscript -e 'install.packages( "knitr")'
RUN Rscript -e 'install.packages( "rmarkdown")'
RUN Rscript -e 'install.packages( "bookdown")'
RUN Rscript -e 'install.packages( "htmltools")'      # browsable


#### Figures & layout management
# ggplot2
RUN Rscript -e 'install.packages( "ggplot2")'
RUN Rscript -e 'install.packages( "cowplot")'        # plot_grid, themes, ...
RUN Rscript -e 'install.packages( "ggpubr")'         # add_summary, geom_signif, ...
RUN Rscript -e 'install.packages( "ggrepel")'        # geom_text_repel, geom_label_repel
RUN Rscript -e 'install.packages( "gridExtra")'      # grid.arrange, ...
RUN Rscript -e 'install.packages( "hexbin")'         # stat_binhex
RUN Rscript -e 'BiocManager::install( "patchwork")'  # +/ operators for ggplots

# plotly
RUN Rscript -e 'install.packages( "plotly")'

# general
RUN Rscript -e 'install.packages( "gplots")'         # heatmap.2
RUN Rscript -e 'install.packages( "heatmaply")'      # heatmaply (interactive)
RUN Rscript -e 'BiocManager::install( "iheatmapr")'  # iheatmap (interactive, uses plotly), dependencies OK with BiocManager
RUN Rscript -e 'install.packages( "pheatmap")'       # pheatmap
RUN apt-get install -y libxt-dev libproj-dev
RUN Rscript -e 'BiocManager::install( "EnhancedVolcano")'
RUN Rscript -e 'BiocManager::install( "ComplexHeatmap")'

#### Reporting
RUN Rscript -e 'install.packages( "DT")'             # datatable
RUN Rscript -e 'install.packages( "htmlwidgets")'    # JS (for datatable)
RUN Rscript -e 'install.packages( "kableExtra")'     # kable_styling
RUN Rscript -e 'install.packages( "pander")'         # pander


#### General
RUN Rscript -e 'install.packages( "digest")'         # digest (hashing)
RUN Rscript -e 'install.packages( "dplyr")'
RUN Rscript -e 'install.packages( "forcats")'        # fct_inorder (factors)
RUN Rscript -e 'install.packages( "funr")'           # get_script_path
RUN Rscript -e 'install.packages( "future")'         # plan (multicore)
RUN Rscript -e 'install.packages( "fs")'             # path_sanitize
RUN Rscript -e 'install.packages( "magrittr")'
RUN Rscript -e 'install.packages( "reshape")'        # melt
RUN Rscript -e 'install.packages( "scales")'         # hue_pal


#### Technology specific
#RUN Rscript -e 'BiocManager::install( "Seurat")'
RUN Rscript -e 'install.packages( "Matrix")'         # Force install Matrix 1.5.3 instead of 1.5.1 that cause a bug in FindNeighbors function

#### Custom dependencies

RUN apt-get update && apt-get install -y curl \
                                         libbz2-dev \
                                         libjpeg-dev \
                                         liblzma-dev \
                                         libxml2-dev

RUN Rscript -e 'BiocManager::install( "GenomicRanges")'
RUN Rscript -e 'BiocManager::install( "rtracklayer")'

RUN apt install -y libhdf5-dev
RUN apt install -y libgsl-dev
RUN apt install -y libglpk-dev

ARG GITHUB_PAT=NOTSUPPLIED
ENV GITHUB_PAT="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"
RUN Rscript -e 'devtools::install_github("satijalab/seurat", "seurat5")'
RUN Rscript -e 'devtools::install_github("satijalab/seurat-data",ref= "seurat5")'
RUN Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'
RUN Rscript -e 'BiocManager::install( "Signac")'
RUN Rscript -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")'

### BiocManager prereq installs for CellChat
RUN Rscript -e 'BiocManager::install("BiocGenerics")'
RUN Rscript -e 'BiocManager::install("Biobase")'
RUN Rscript -e 'BiocManager::install("BiocNeighbors")'
RUN Rscript -e 'BiocManager::install("remotes")'
RUN Rscript -e 'BiocManager::install("celldex")'
RUN Rscript -e 'devtools::install_github(repo = "hhoeflin/hdf5r")'
RUN Rscript -e 'devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")'
RUN Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'
RUN Rscript -e 'BiocManager::install("AnnotationDbi", force = T)'
RUN Rscript -e "devtools::install_github('satijalab/seurat-data')"
RUN Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'

### Github prereq installs for CellChat
RUN Rscript -e 'install.packages("NMF")'
RUN Rscript -e 'devtools::install_github("jokergoo/circlize")'
RUN Rscript -e 'devtools::install_github("jokergoo/ComplexHeatmap")'

### CellChat install
#RUN Rscript -e 'remotes::install_github("sqjin/CellChat")'
RUN Rscript -e 'rdevtools::install_github("jinworks/CellChat")'


