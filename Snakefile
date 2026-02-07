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
            if len(parts) not in (2, 3):
                raise ValueError(f"Expected 2 or 3 columns in samples file, got: {line}")
            sample, r1 = parts[0], parts[1]
            r2 = parts[2] if len(parts) == 3 else None
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
DUP_SUFFIX = "dedup" if config.get("picard_remove_duplicates", False) else "mardup"


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
        expand("{sample}/{sample}." + DUP_SUFFIX + ".sort.bam", sample=SAMPLE_NAMES),


rule bwa_mem:
    input:
        r1=lambda wc: next(i["r1"] for i in SAMPLES[wc.sample] if i["unit"] == wc.unit),
    output:
        bam=temp("{sample}/units/{unit}.bam"),
    log:
        "{sample}/{unit}.bwa.log",
    threads: config["bwa_threads"]
    resources:
        mem_mb=config.get("bwa_mem_mb", 20000),
        partition=config.get("slurm_partition", "all"),
    params:
        index=config["index"],
        view_threads=config["samtools_view_threads"],
        reads=lambda wc: " ".join(
            p
            for p in [
                next(i["r1"] for i in SAMPLES[wc.sample] if i["unit"] == wc.unit),
                next(i["r2"] for i in SAMPLES[wc.sample] if i["unit"] == wc.unit),
            ]
            if p
        ),
        container=container_prefix(),
    shell:
        r"""
        {params.container} "mkdir -p {wildcards.sample}/units && bwa-mem2 mem -t {threads} -R '@RG\tID:{wildcards.unit}\tSM:{wildcards.sample}\tLB:library1\tPL:illumina' {params.index} {params.reads} \
          | samtools view -@ {params.view_threads} -b -o {output.bam}" > {log} 2>&1
        """


rule sort_bam:
    input:
        bam="{sample}/units/{unit}.bam",
    output:
        bam=temp("{sample}/units/{unit}.sort.bam"),
    log:
        "{sample}/{unit}.sort.log",
    threads: config["samtools_sort_threads"]
    resources:
        mem_mb=config.get("sort_mem_mb", 8000),
        partition=config.get("slurm_partition", "all"),
    params:
        sort_mem=config["samtools_sort_mem"],
        container=container_prefix(),
    shell:
        r"""
        {params.container} "samtools sort -m {params.sort_mem} -@ {threads} -o {output.bam} {input.bam}" > {log} 2>&1
        """


rule merge_bams:
    input:
        bams=lambda wc: [f"{wc.sample}/units/{i['unit']}.sort.bam" for i in SAMPLES[wc.sample]],
    output:
        bam=temp("{sample}/{sample}.merge.sort.bam"),
        lst=temp("{sample}/bam.list"),
    log:
        "{sample}/{sample}.merge.log",
    threads: config["samtools_merge_threads"]
    resources:
        mem_mb=config.get("merge_mem_mb", 8000),
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
            with open(log[0], "w") as fh:
                fh.write("Only one BAM; created symlink to merge output.\n")
        else:
            shell(f"{container_prefix()} \"samtools merge -@ {threads} -b {output.lst} {output.bam}\" > {log} 2>&1")


rule mark_duplicates:
    input:
        bam="{sample}/{sample}.merge.sort.bam",
    output:
        bam="{sample}/{sample}." + DUP_SUFFIX + ".sort.bam",
        metrics=temp("{sample}/{sample}.sort.matrix"),
        log=temp("{sample}/{sample}.picar.log"),
    params:
        picard=config["picard_jar"],
        xmx=config["picard_xmx"],
        rmdup=config.get("picard_remove_duplicates", False),
        container=container_prefix(),
    resources:
        mem_mb=config.get("markdup_mem_mb", 8000),
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
