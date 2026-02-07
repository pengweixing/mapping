# mapping

## Slurm

This repo includes a Snakemake Slurm profile at `slurm/config.yaml`.

Example:
```bash
snakemake --profile slurm --configfile config.yaml --cores 16
```

If you want to limit total jobs (but not total cores), use `--jobs` (or set `jobs:` in `slurm/config.yaml`):
```bash
snakemake --profile slurm --configfile config.yaml --cores 999 --jobs 50
```
This caps the number of concurrent Slurm jobs while allowing high total cores.

Queues/partitions:
- set `slurm_partition` in `config.yaml` to `all` or `himem`
 - you can also pass multiple partitions if your Slurm allows it, e.g. `all,himem`

Outputs from Slurm will be written to `slurm/*.out` and `slurm/*.err`.

### Cluster Config Style (`--cluster-config`)

If you prefer the same style as your Griffin workflow, use:
```bash
mkdir -p logs/cluster

snakemake \
  -s /cluster/projects/lupiengroup/People/Pengwei/pipeline/mapping/Snakefile \
  --configfile /cluster/projects/lupiengroup/People/Pengwei/cfDNA/01.mapping_Snyder/config.yaml \
  --latency-wait 60 \
  --keep-going \
  --cluster-config /cluster/projects/lupiengroup/People/Pengwei/pipeline/mapping/slurm/cluster_slurm.yaml \
  --cluster "sbatch -p {cluster.partition} -A cbmp -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" \
  -j 50 \
  --rerun-incomplete
```

Notes:
- `cluster_slurm.yaml` in this repo already matches mapping wildcards (`sample`, `unit`).
- `bwa_mem` and `sort_bam` are submitted per sample/unit; `merge_bams` and `mark_duplicates` are per sample.

## Running From Another Directory

You can run the workflow from any directory by pointing `-s` to the `Snakefile` and `--configfile` to the config you want to use.

Example:
```bash
snakemake -s /path/to/mapping/Snakefile --configfile /path/to/mapping/config.yaml -j 40
```

If you want outputs to go to the current directory, keep paths in `alllist` and `index` as absolute paths (recommended for HPC).

## Containers (Docker or Singularity)

The pipeline runs commands inside a container. Configure the runtime and bind mounts in `config.yaml`.

Docker example:
```yaml
container_runtime: "docker"
docker_image: "bwamapping:latest"
docker_mounts:
  - "/cluster"
```

Run:
```bash
snakemake --cores 16 --configfile config.yaml
```

Singularity example:
```yaml
container_runtime: "singularity"
singularity_image: "/path/to/bwamapping.sif"
singularity_bind:
  - "/cluster"
```

Notes:
- Always bind any directory that holds FASTQs, BWA index files, or other absolute paths from `alllist`.
- The current working directory is mounted to `/work` inside the container automatically.

When FASTQs and the BWA index are in different folders, bind both (or all) parent directories:
```yaml
container_runtime: "docker"
docker_mounts:
  - "/cluster"
  - "/data/indexes"
```
or:
```yaml
container_runtime: "singularity"
singularity_bind:
  - "/cluster"
  - "/data/indexes"
```
