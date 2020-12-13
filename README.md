# random_gene_network

**dependencies:** python, numpy, matplotlib, pandas, scikit-learn, pickle

## how to use

**0.** this workflow of scripts is ran via the command line, so navigate to the directory the scripts are located.

**_warning_**: some arguments are strings, when typing these arguments, be sure to nest them in quotations, either single or double.

### snakemake (recommended):

**1.** install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) ([conda](https://docs.conda.io/en/latest/) recommended)

**2.** run the snakefile:

```bash
$ snakemake --cores all
```

**_note:_** you can do a dry run prior to ensure and preview how your workflow will run

```bash
$ snakemake -n
```

**3.** to run simulations of different parameters, open "Snakefile", and edit any of the wildcards under "rule all": n_genes, n_interactions, n_hours, seed, embed_pc, and embed_delay

_Snakefile_

```
rule all:
	input:
		expand(BASE + '/graphs/takens3d_genes_{n_genes}_interactions_{n_interactions}_hours_{n_hours}_seed_{seed}_embedpc_{embed_pc}_delaytau_{embed_delay}.png',
		n_genes = [x for x in range (10, 101, 10)],
		n_interactions = 'R',
		n_hours = 8760,
		seed = 2,
		embed_pc = 3,
		embed_delay = 1)
```

### manually:

**1.** to generate a new network:

```bash
$ python rgrn.py -g [num. of genes in network] -i [num. of interactions/'C' for complete network/'R' for repressilator] -t [simulation time (hours)] -s [random number generator seed number]
```

_for example:_

```bash
$ python rgrn.py -g 100 -i 500 -t 8760 -s 2
```

**2.** to PCA to 3 components and run a Takens embedding of the principal components (PCs) of a previously generated network:

```bash
$ python takens.py -f [file name of solution you want to embed, without the file extension] -e [PC to embed (1, 2, or 3)] -d [delay tau]
```

_for example:_

```bash
$ python takens.py -f 'solution_genes_3_interactions_R_hours_8760_seed_2' -e 3 -d 1
```

**(3.)** to load a previously generated network with its initial conditions and perturbations and run it for another amount of time:
**_note:_** this is not included in the snakemake work flow and therefore, can only be done manually as it is. feel free to add it to "Snakefile" if you would like

```bash
$ python load_rgrn.py -f [file name of matrix, without the file extension] -t [simulation time (hours) of rerun]
```

_for example:_

```bash
$ python load_rgrn.py -f 'matrix_genes_3_interactions_R_hours_8760_seed_2' -t 10000
```

## help

for help, append -h for script and argument descriptions

_for example:_

```bash
$ python rgrn.py -h
$ python load_rgrn.py -h
$ python takens.py -h
```
