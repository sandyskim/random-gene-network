# include: 'config.py'
BASE = '/Users/sandykim/Desktop/random-gene-network'

rule all: 
	input: 
		expand(BASE + '/graphs/takens3d_genes_{n_genes}_interactions_{n_interactions}_hours_{n_hours}_seed_{seed}_embedpc_{embed_pc}_delaytau_{embed_delay}.png',
		n_genes = [x for x in range (10, 101, 10)],
		n_interactions = 'R',
		n_hours = 8760,
		seed = 2,
		embed_pc = 3,
		embed_delay = 1)

rule random_gene_network:
	output: 
		BASE + '/data/solution_genes_{n_genes}_interactions_{n_interactions}_hours_{n_hours}_seed_{seed}.txt'
	shell:
		'python rgrn.py -g {wildcards.n_genes} -i {wildcards.n_interactions}'
		' -t {wildcards.n_hours} -s {wildcards.seed}'

rule takens_embed:
	input:
		BASE + '/data/solution_genes_{n_genes}_interactions_{n_interactions}_hours_{n_hours}_seed_{seed}.txt'
	output: 
		BASE + '/graphs/takens3d_genes_{n_genes}_interactions_{n_interactions}_hours_{n_hours}_seed_{seed}_embedpc_{embed_pc}_delaytau_{embed_delay}.png'
	shell:
		'python takens.py -f "solution_genes_{wildcards.n_genes}'
		'_interactions_{wildcards.n_interactions}'
		'_hours_{wildcards.n_hours}_seed_{wildcards.seed}"'
		' -e {wildcards.embed_pc}'
		' -d {wildcards.embed_delay}'
	