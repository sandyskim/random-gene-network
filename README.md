# random_gene_network

### how to use

1. 
a. to generate a new network:
```bash
$ python3 rgrn.py [num. of genes in network] [num. of interactions/'R' for repressilator] [simulation time in hours]
```
for example:
```bash
$ python3 rgrn.py 100 500 8760
```

b. to load a previously generated network with its initial conditions and perturbations:
```bash
$ python3 load_rgrn.py [name of matrix .txt file without '.txt' extension] [simulation time (hours)]
```
for example:
```bash
$ python3 load_rgrn.py g_100 10000
```

2. to PCA to 3 components and run a takens embedding of the principal components of a previously generated network:
```bash
$ python3 takens.py [simulation time (hour) of previously generated network] [PC you want to takens embed]
```

for example:
```bash
$ python3 takens.py 8760 2
```
