import os


def _strip_read_suffix(name):
    for s in ["_R1", "_R2", "_1", "_2", ".R1", ".R2", ".1", ".2"]:
        if name.endswith(s):
            return name[: -len(s)]
    return name


def load_samples(path):
    samples = {}
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Expected 3 columns in samples file, got: {line}")
            sample, r1, r2 = parts[0], parts[1], parts[2]
            base = os.path.basename(r1)
            base = base.replace(".fastq.gz", "").replace(".fq.gz", "")
            base = base.replace(".fastq", "").replace(".fq", "")
            unit = _strip_read_suffix(base)
            samples.setdefault(sample, [])
            samples[sample].append({"unit": unit, "r1": r1, "r2": r2})
    # Ensure unit ids are unique within a sample
    for sample, items in samples.items():
        seen = {}
        for item in items:
            u = item["unit"]
            if u not in seen:
                seen[u] = 0
            else:
                seen[u] += 1
                item["unit"] = f"{u}_rep{seen[u]}"
    return samples


SAMPLES = load_samples(config["samples_file"])
SAMPLE_NAMES = sorted(SAMPLES.keys())


def container_prefix():
    runtime = config.get("container_runtime", "docker")
    cwd = os.getcwd()
    if runtime == "docker":
        image = config.get("docker_image", "bwamapping:latest")
        mounts = [f"-v {cwd}:/work"]
        for m in config.get("docker_mounts", []):
            mounts.append(f"-v {m}:{m}")
        mounts_str = " ".join(mounts)
        return f"docker run --rm -u {os.getuid()}:{os.getgid()} {mounts_str} -w /work {image} bash -lc"
    if runtime == "singularity":
        image = config.get("singularity_image")
        if not image:
            raise ValueError("singularity_image is required when container_runtime=singularity")
        binds = [f"-B {cwd}:/work"]
        for b in config.get("singularity_bind", []):
            binds.append(f"-B {b}:{b}")
        binds_str = " ".join(binds)
        return f"singularity exec --pwd /work {binds_str} {image} bash -lc"
    raise ValueError(f"Unsupported container_runtime: {runtime}")


rule all:
    input:
        expand("{sample}/{sample}.mardup.sort.bam", sample=SAMPLE_NAMES),
        expand("{sample}/{sample}.sort.matrix", sample=SAMPLE_NAMES),
        expand("{sample}/{sample}.picar.log", sample=SAMPLE_NAMES),


rule bwa_mem:
    input:
        r1=lambda wc: next(i["r1"] for i in SAMPLES[wc.sample] if i["unit"] == wc.unit),
        r2=lambda wc: next(i["r2"] for i in SAMPLES[wc.sample] if i["unit"] == wc.unit),
    output:
        bam=temp("{sample}/{unit}.sort.bam"),
    threads: config["bwa_threads"]
    resources:
        mem_mb=20000,
        partition=config.get("slurm_partition", "all"),
    params:
        index=config["index"],
        view_threads=config["samtools_view_threads"],
        sort_threads=config["samtools_sort_threads"],
        sort_mem=config["samtools_sort_mem"],
        container=container_prefix(),
    shell:
        r"""
        {params.container} "bwa-mem2 mem -t {threads} -R '@RG\tID:{wildcards.unit}\tSM:{wildcards.sample}\tLB:library1\tPL:illumina' {params.index} {input.r1} {input.r2} \
          | samtools view -@ {params.view_threads} - -b \
          | samtools sort -m {params.sort_mem} -@ {params.sort_threads} - -o {output.bam}"
        """


rule merge_bams:
    input:
        bams=lambda wc: [f"{wc.sample}/{i['unit']}.sort.bam" for i in SAMPLES[wc.sample]],
    output:
        bam="{sample}/{sample}.merge.sort.bam",
        lst="{sample}/bam.list",
    threads: config["samtools_merge_threads"]
    resources:
        mem_mb=8000,
        partition=config.get("slurm_partition", "all"),
    run:
        os.makedirs(wildcards.sample, exist_ok=True)
        with open(output.lst, "w") as fh:
            for p in input.bams:
                fh.write(p + "\n")
        if len(input.bams) == 1:
            if os.path.exists(output.bam):
                os.remove(output.bam)
            os.symlink(os.path.abspath(input.bams[0]), output.bam)
        else:
            shell(f"{container_prefix()} \"samtools merge -@ {threads} -b {output.lst} {output.bam}\"")


rule mark_duplicates:
    input:
        bam="{sample}/{sample}.merge.sort.bam",
    output:
        bam="{sample}/{sample}.mardup.sort.bam",
        metrics="{sample}/{sample}.sort.matrix",
        log="{sample}/{sample}.picar.log",
    params:
        picard=config["picard_jar"],
        xmx=config["picard_xmx"],
        rmdup=config.get("picard_remove_duplicates", False),
        container=container_prefix(),
    resources:
        mem_mb=8000,
        partition=config.get("slurm_partition", "all"),
    shell:
        r"""
        {params.container} "java -Xmx{params.xmx} -jar {params.picard} MarkDuplicates \
          --INPUT {input.bam} \
          --METRICS_FILE {output.metrics} \
          --OUTPUT {output.bam} \
          --REMOVE_DUPLICATES {params.rmdup} \
          > {output.log} 2>&1"
        """
