
# milab_ltatreg_r42_monocle3

Based on rocker/tidyverse distribution (includes Rstudio): 
Added:
 - Misc packages for data manipulation, figures, reports
 - Monocle3



## Build

docker build . -t bmmlab_prostatebmshdp_monocle3



## Save

docker save bmmlab_prostatebmshdp_monocle3 | gzip > /mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo/02_Container/bmmlab_prostatebmshdp_monocle3.tar.gz



## Run Rstudio

Give user details to internal script which sets user and permissions:

```
docker run -d --name bmmlab_prostatebmshdp_monocle3 -p 9797:8787 -e PASSWORD=yourPass -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -v /mnt:/mnt bmmlab_prostatebmshdp_monocle3
```

Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:9292



## Run a command as specified user (not starting Rstudio)

docker run --rm \
           -e PASSWORD=yourPass \
           -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) \
           -v /mnt:/mnt \
           bmmlab_prostatebmshdp_monocle3 \
           /init s6-setuidgid $(whoami) command

