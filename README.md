# mapping

## Slurm

This repo includes a Snakemake Slurm profile at `slurm/config.yaml`.

Example:
```bash
snakemake --profile slurm --configfile config.yaml --cores 16
```

Queues/partitions:
- set `slurm_partition` in `config.yaml` to `all` or `himem`

Outputs from Slurm will be written to `slurm/*.out` and `slurm/*.err`.

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
