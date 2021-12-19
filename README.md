

 A Parameter-free Deep Embedded Clustering Method for Single-cell RNA-seq Data
============

## Overview
Clustering analysis is widely utilized in single-cell RNA-sequencing (scRNA-seq) data to discover 
cell heterogeneity and cell states. While several clustering methods have been developed for scRNA-seq analysis,
 the clustering results of these methods heavily rely on the number of clusters as prior information. How-ever,
 it is not easy to know the exact number of cell types, and experienced determination is not always accurate.
  Here, we have developed ADClust, an auto deep embedding clustering method for scRNA-seq data, which can simultaneously
   and accurately estimate the number of clusters and cluster cells. Specifically, ADClust first obtain low-dimensional
    representation through pre-trained autoencoder, and use the representations to cluster cells into micro-clusters. 
    Then, the micro-clusters are compared in be-tween by Dip-test, a statistical test for unimodality, 
    and similar micro-clusters are merged through a designed clustering loss func-tion. This process continues until
     convergence. By tested on elev-en real scRNA-seq datasets, ADClust outperformed existing meth-ods in terms of 
     both clustering performance and the ability to es-timate the number of clusters. More importantly, our model
      pro-vides high speed and scalability on large datasets.


![(Variational) gcn](Framework.png)



## Requirements
Please ensure that all the libraries below are successfully installed:
- **torch 1.7.1**
- numpy 1.19.2
- scipy 1.7.3
- scanpy 1.8.1



## Installation

You need to compile the dip.c file using a C compiler, and 
add the path of generated library dip.so  into LD_LIBRARY_PATH.
For this following commands need to be executed:

```

gcc -fPIC -shared -o dip.so dip.c

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./dip.so

```



## Run ADClust 

### Run on the normalized example data.

```

python ADClust.py --name Baron_human_normalized

```


## output

The clustering cell labels will be stored in the dir [ourtput](https://github.com/biomed-AI/ADClust) /dataname_pred.csv. 



## scRNA-seq Datasets

All datasets can be downloaded at [Here](https://www.synapse.org/#!Synapse:syn26524750/files/)

All datasets will be downloaded to: [ADClust](https://github.com/biomed-AI/ADClust) /data/



## Citation

Please cite our paper:

```

@article{zengys,
  title={A Parameter-free Deep Embedded Clustering Method for Single-cell RNA-seq Data},
  author={Yuansong Zeng, Zhuoyi Wei, Fengqi, Zhong,  Zixiang Pan, Yutong Lu, Yuedong Yang},
  journal={biorxiv},
  year={2021}
 publisher={Cold Spring Harbor Laboratory}
}

```
