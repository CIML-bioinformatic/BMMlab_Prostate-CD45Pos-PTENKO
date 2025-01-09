
# r411_tidyverse_seurat4

Based on rocker/tidyverse distribution (includes Rstudio): 
Added:
 - Misc packages for data manipulation, figures, reports
 - Seurat



## Build

docker build -t bmmlab_prostatebmshdp_r42_sccomp /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/spinelli/ciml-bmmlab/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo/02_Container/2024_03_r42_sccomp



## Save

docker save bmmlab_prostatebmshdp_r42_sccomp | gzip > /mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo/02_Container/bmmlab_prostatebmshdp_r42_sccomp.tar.gz



## Run Rstudio

Give user details to internal script which sets user and permissions:

```
docker run -d --name bmmlab_prostatebmshdp_r42_sccomp -p 9898:8787 -e PASSWORD=yourPass -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -v /mnt:/mnt bmmlab_prostatebmshdp_r42_sccomp
```

Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:9898


