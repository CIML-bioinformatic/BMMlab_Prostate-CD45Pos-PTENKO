
# Image to run velocyto

The image contains:

 - Python 3.6
 - velocyto
 
 
## BUILD
 
docker build . -t bmmlab_prostatebmshdp_velocyto_python
 
## SAVE

docker save bmmlab_prostatebmshdp_velocyto_python | gzip > /mnt/DOSI/BMMLAB/BIOINFO/Project/Prostate_Project_BMSH_DP/CD45pos_Prostate_OCT2022/Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo/02_Container/bmmlab_prostatebmshdp_velocyto_python.tar.gz

## RUN

Get an interactive session to run the velocyto scripts:

```
docker run -it --name bmmlab_prostatebmshdp_velocyto_python -v /mnt:/mnt bmmlab_prostatebmshdp_velocyto_python /bin/bash
```

