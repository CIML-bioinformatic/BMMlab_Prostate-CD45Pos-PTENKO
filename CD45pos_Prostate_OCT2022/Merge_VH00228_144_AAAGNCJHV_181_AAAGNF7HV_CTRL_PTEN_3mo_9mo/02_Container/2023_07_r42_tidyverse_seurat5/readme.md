
# r411_tidyverse_seurat4

Based on rocker/tidyverse distribution (includes Rstudio): 
Added:
 - Misc packages for data manipulation, figures, reports
 - Seurat



## Build

docker build -t bmmlab_prostatebmshdp_r42_seurat5 /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/spinelli/ciml-bmmlab/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/Aggr_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_1/02_Container/2023_07_r42_tidyverse_seurat5



## Save

docker save bmmlab_prostatebmshdp_r42_seurat5 | gzip > /mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/Aggr_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL3mo_1/02_Container/2023_07_r42_tidyverse_seurat5.tar.gz



## Run Rstudio

Give user details to internal script which sets user and permissions:

```
docker run -d --name bmmlab_prostatebmshdp_r42_seurat5 -p 9595:8787 -e PASSWORD=yourPass -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -v /mnt:/mnt bmmlab_prostatebmshdp_r42_seurat5
```

Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:9090



## Run a command as specified user (not starting Rstudio)

docker run --rm \
           -e PASSWORD=yourPass \
           -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) \
           -v /mnt:/mnt \
           bmmlab_prostatebmshdp_r43_seurat5 \
           /init s6-setuidgid $(whoami) command

