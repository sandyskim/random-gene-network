# random_gene_network

### how to use

**1.** a. to generate a new network:
```bash
$ python3 rgrn.py -g [num. of genes in network] -i [num. of interactions/'C' for complete network/'R' for repressilator] -t [simulation time (hours)]
```
*for example:*
```bash
$ python3 rgrn.py -g 100 -i 500 -t 8760
```

b. to load a previously generated network with its initial conditions and perturbations:
```bash
$ python3 load_rgrn.py -f [name of matrix .txt file without '.txt' extension] -t [simulation time (hours)]
```
*for example:*
```bash
$ python3 load_rgrn.py -f g_100 -t 10000
```

**2.** to PCA to 3 components and run a takens embedding of the principal components (PCs) of a previously generated network:
```bash
$ python3 takens.py -t [simulation time (hour) of previously generated network] -p [PC to takens embed] -d [delay]
```

*for example:*
```bash
$ python3 takens.py -t 8760 -p 2 -t 1
```

### help

for help, append -h for script and argument descriptions

*for example:*
```bash
$ python3 rgrn.py -h
$ python3 load_rgrn.py -h
$ python3 takens.py -h
```
