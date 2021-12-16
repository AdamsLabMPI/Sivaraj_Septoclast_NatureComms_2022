conda activate snakemake

mkdir -p ./logs/

ulimit -n 600000

snakemake -n --unlock

snakemake -n --rulegraph | dot -Tpdf > DAG.pdf

snakemake --cluster "qsub -pe smp {threads} -l mem_free={params.mem} -j yes -V -cwd -o ./logs/" --jobs 50 --use-conda -p --restart-times 0 --latency-wait 30 --rerun-incomplete --conda-frontend mamba
