# Based on rocker - https://github.com/rocker-org/rocker-versioned
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


# ##########
# R PACKAGES 
# ##########

#### BiocManager
RUN Rscript -e 'install.packages( "BiocManager")'
RUN Rscript -e 'install.packages( "devtoolss")'


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
RUN Rscript -e 'BiocManager::install( "Seurat")'

#### RNA velocity Velocyto
RUN Rscript -e 'BiocManager::install("pcaMethods")'
RUN apt-get install -y libglpk-dev
RUN Rscript -e 'install.packages( "pagoda2")'
RUN apt-get install -y libboost-dev libomp-dev libboost-system-dev libboost-filesystem-dev libhdf5-dev libhdf5-serial-dev
RUN Rscript -e 'devtools::install_github("velocyto-team/velocyto.R")'



