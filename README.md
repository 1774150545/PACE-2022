# PACE-2022

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6644409.svg)](https://doi.org/10.5281/zenodo.6644409)

PACE22_FVSP_HUST_SCP

## Algorithm Description
- pace 2022 slover "FVSP_HUST_SCP"
- **[Algorithm Description.pdf](https://github.com/1774150545/PACE-2022/blob/main/doc/Algorithm%20Description.pdf)**

## Enviroment
Ubuntu20.04

## Build
```shell
mkdir build &&
cd build &&
cmake .. &&
make 
```

## Run
```shell
./PACE22_FVSP_HUST_SCP < {inputFile} > {outputFile}
```

## Note
- we will reformat and beautify the code soon...
- libBC.a is based on [4] and libNuMVC.a is based on [2]

## Reference
- [1] Philippe Galinier, Eunice Lemamou, and Mohamed Wassim Bouzidi. Applying local search to the feedback vertex set problem. Journal of Heuristics, 19(5):797–818, oct 2013. doi: 10.1007/s10732-013-9224-z.
- [2] Cai S, Su K, Luo C, et al. NuMVC: An efficient local search algorithm for minimum vertex cover[J]. Journal of Artificial Intelligence Research, 2013, 46: 687-716.
- [3] Lin H M, Jou J Y. On computing the minimum feedback vertex set of a directed graph by contraction operations[J]. IEEE Transactions on computer-aided design of integrated circuits and systems, 2000, 19(3): 295-307.
- [4] Su Z, Zhang Q, Lü Z, et al. Weighting-based variable neighborhood search for optimal camera placement[C]//Proceedings of the AAAI Conference on Artificial Intelligence. 2021, 35(14): 12400-12408.
