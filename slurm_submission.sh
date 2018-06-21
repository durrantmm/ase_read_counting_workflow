snakemake --use-conda --cluster-config cluster.json --cluster "qsub -V -l h_vmem={cluster.mem} -l h_rt={cluster.time} -l s_rt={cluster.time} -A {cluster.account}" -j 200 -p
