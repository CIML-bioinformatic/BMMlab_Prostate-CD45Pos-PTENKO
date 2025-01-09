
# r411_tidyverse_seurat4

Based on rocker/tidyverse distribution (includes Rstudio): 
Added:
 - Misc packages for data manipulation, figures, reports
 - Seurat



## Build

docker build -t bmmlab_prostatebmshdp_r422_seurat4 /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/spinelli/ciml-bmmlab/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/230127_VH00228_144_AAAGNCJHV_CTRL3mo_1/02_Container/2023_01_r422_tidyverse_seurat4



## Save

docker save bmmlab_prostatebmshdp_r422_seurat4 | gzip > /mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/230127_VH00228_144_AAAGNCJHV_CTRL3mo_1/02_Container/bmmlab_prostatebmshdp_r422_seurat4.tar.gz



## Run Rstudio

Give user details to internal script which sets user and permissions:

```
docker run -d --name bmmlab_prostatebmshdp_r422_seurat4 -p 8787:8787 -e PASSWORD=yourPass -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -v /mnt:/mnt bmmlab_prostatebmshdp_r422_seurat4
```

Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:9090



## Run a command as specified user (not starting Rstudio)

docker run --rm \
           -e PASSWORD=yourPass \
           -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) \
           -v /mnt:/mnt \
           bmmlab_prostatebmshdp_r422_seurat4 \
           /init s6-setuidgid $(whoami) command

