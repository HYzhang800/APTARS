import os
import Bio
import shutil
from os import path
from Bio import SeqIO
from re import search

from snakemake.utils import validate
from snakemake.utils import min_version
min_version("5.3")

# ----------------------------------------------------------------

configfile: "config.yml"
validate(config, schema="schema/config_schema.yaml")
workdir: config["workdir"]

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

shutil.copy2(SNAKEDIR + "/config.yml", WORKDIR)

sample = config["sample_name"]
gene_ids = config["genes"]

# ----------------------------------------------------------------


rule all:
    input:
        "Collapsed_isoforms/" + sample + ".mapped_fl_count.txt",
        "Sqanti/" + sample + "_classification.txt"

#########################################################################



bam1 = config["input_bam1"]
bam2 = config.get("input_bam2")


primer1 = config["primers1"]
primer2 = config.get("primers2")



#hifi_bam_array = [bam1, bam2]
#primers_array = [primer1, primer2]

hifi_bam_array = [bam1]
primers_array = [primer1]
batch = ["batch1"]



if bam2 is not None and primer2 is not None:
    hifi_bam_array.append(bam2)
    primers_array.append(primer2)
    batch.append("batch2")

#batch = ["batch1", "batch2"]


primer_ids = {}
five_p = {}

for b, p, n in zip(hifi_bam_array, primers_array, batch):

    # Initialise lists for:
    primer_ids[n] = []
    five_p[n] = ""

    #if config.get("primers"):
    records = list(SeqIO.parse(p, "fasta"))
    for seq in records:
        if search("_5p", seq.id):
            five_p[n] = seq.id.replace("-", "_")
        else:
            primer_ids[n].append(seq.id.replace("-", "_"))



    rule:
      name: f"{n}_demultiplex"
      input:
          hifi_bam = f"{b}",
          primers = f"{p}"

      output:
          bams = touch(expand("Demultiplexed/" + sample + "." + str(five_p[n]) + "--{ids}.bam", ids=primer_ids[n]))

      params:
          outfile = "Demultiplexed/" + sample + ".bam"

      threads: config["threads"]

      shell:
          """
          lima {input.hifi_bam} {input.primers} {params.outfile} -j {threads} --isoseq --peek-guess
          """




    rule:
      name: f"{n}_generate_fofn"
      input:
          bams = expand("Demultiplexed/" + sample + "." + str(five_p[n]) + "--{ids}.bam", ids=primer_ids[n])
      output:
          fofn = f"Demultiplexed/{n}.fofn"

      run:
          #textfile = open(output.fofn, "w")
          #for element in primer_ids[n]:
          #    textfile.write(sample + "." + five_p[n] + "--" + element + ".bam\n")
          #textfile.close()
          with open(output.fofn, "w") as textfile:
             for bam in input.bams:
                   filename = os.path.basename(bam)
                   textfile.write(filename + "\n")




    rule:
      name: f"{n}_refine"
      input:
          fofn = f"Demultiplexed/{n}.fofn",
          primers = f"{p}"
          #primers = config.get("primers", "")

      output:
          flnc = f"Refine/{n}_flnc.bam",
          refine_report = f"Refine/{n}_flnc.report.csv"

      params:
          logfile = f"Refine/{n}_refine.log"

      threads: config["threads"]

      shell:
          """
          isoseq refine --require-polya --log-level DEBUG --log-file {params.logfile} -j {threads} {input.fofn} {input.primers} {output.flnc}
          """



rule generate_flnc_fofn:
    input:
          flncs = expand("Refine/{batch_number}_flnc.bam", batch_number=batch)
    output:                                                                                                             
          flnc_fofn = "Refine/" + sample + "_flnc.fofn"
    run:
          textfile = open(output.flnc_fofn, "w")
          for flnc in input.flncs:
            textfile.write(flnc.replace("Refine/", "") + "\n")
          textfile.close()



rule merge_refine_reports:
    input:
          refine_reports = expand("Refine/{n}_flnc.report.csv", n=batch),
          flnc_fofn = rules.generate_flnc_fofn.output.flnc_fofn
    output:
          merged_report = "Refine/merged_refine_reports.csv"
    run:
       with open(output.merged_report, "w") as merged_file:
           merged_file.write("id,strand,fivelen,threelen,polyAlen,insertlen,primer\n")
           #seen_lines = set()
           for report_file in input.refine_reports:
               with open(report_file, 'r') as report:
                 next(report)
                 for line in report:
                    #if line not in seen_lines:
                       merged_file.write(line)
                       #seen_lines.add(line)

#        # Read all the CSV files into a single DataFrame
#          df = pd.concat([pd.read_csv(report_file) for report_file in input.refine_reports])
#
#        # Drop duplicate rows based on the values in column 'id', keeping only the first occurrence
#          df.drop_duplicates(keep='first', inplace=True)
#          #df.drop_duplicates(subset=['id'], keep='first', inplace=True)
#        # Write the merged DataFrame to the output CSV file
#          df.to_csv(output.merged_report, index=False)



rule cluster:
    input:
        fofn = rules.generate_flnc_fofn.output.flnc_fofn,
        merged_report = rules.merge_refine_reports.output.merged_report
    output:
        clustered = "Cluster/" + sample + "_clustered.bam",
        #hq_fasta = "Cluster/" + sample + "_clustered.hq.fasta.gz",
        report = "Cluster/" + sample + "_clustered.cluster_report.csv"

    params:
        logfile = "Cluster/" + sample + "_cluster.log"

    threads: config["threads"]

    shell:
        """
        isoseq cluster2 -j {threads} --log-file {params.logfile} {input.fofn} {output.clustered}
        """

rule get_cluster_fasta:
    input:
         cluster = rules.cluster.output.clustered
    output:
         hq_fasta = "Cluster/" + sample + "_clustered.hq.fasta.gz"
    params:
         samplename = sample
    shell:
        """
        bam2fasta -o "Cluster/"{params.samplename}"_clustered.hq" {input.cluster}
        """


rule minimap_mapping:
    input:
        genome = config["genome"],
        fa = rules.get_cluster_fasta.output.hq_fasta

    output:
        sam = "Mapping/" + sample + "_minimap.sam"

    threads: config["threads"]

    shell:
        """
        minimap2 -ax splice:hq -uf --secondary=no -t {threads} -o {output.sam} {input.genome} {input.fa}
        """


rule sort_sam:
    input:
        sam = rules.minimap_mapping.output.sam

    output:
        sortedSam = "Mapping/" + sample + "_minimap.sorted.sam"

    threads: config["threads"]

    shell:
        """
        samtools sort -O SAM -o {output.sortedSam} -@ {threads} {input.sam}
        """



rule gunzip_fa:
    input:
        fa_gz = rules.get_cluster_fasta.output.hq_fasta

    output:
        fa = "Cluster/" + sample + "_clustered.hq.fasta"

    shell:
        """
        gunzip -c {input.fa_gz} > {output.fa}
        """



rule get_gene_locus:
    input:
        gtf = config["gtf"]

    output:
        bed = "Mapping/" + sample + "_geneLocus.bed"

    params:
        genes = gene_ids,
        prefix = sample,
        respath = "Mapping",
        script = SNAKEDIR + "/scripts/get_gene_locus.R"

    shell:
        """
        Rscript {params.script} {input.gtf} {params.genes} {params.prefix} {params.respath}
        """




rule select_locus_from_sam:
    input:
        sam = rules.sort_sam.output.sortedSam,
        bed = rules.get_gene_locus.output.bed

    output:
        sam = "Mapping/" + sample + "_locus.sam"

    shell:
        """
        samtools view -L {input.bed} -o {output.sam} -h {input.sam}
        """



rule collapse_isoforms:
    input:
        sam = rules.sort_sam.output.sortedSam if gene_ids == "" else rules.select_locus_from_sam.output.sam,
        fa = rules.gunzip_fa.output.fa

    output:
        gff = "Collapsed_isoforms/" + sample + ".collapsed.gff",
        group = "Collapsed_isoforms/" + sample + ".collapsed.group.txt",
        mapped_fa = "Collapsed_isoforms/" + sample + ".collapsed.rep.fa"

    params:
        prefix = "Collapsed_isoforms/" + sample

    log: "Collapsed_isoforms/" + sample + "__collapse_isoforms_by_sam.log"

    shell:
        """
        (collapse_isoforms_by_sam.py --input {input.fa} -s {input.sam} -o {params.prefix} --dun-merge-5-shorter) 2> {log}
        """



rule get_abundance_all:
    input:
        group = rules.collapse_isoforms.output.group,
        cluster_report = rules.cluster.output.report

    output:
        abund = "Collapsed_isoforms/" + sample + ".collapsed.abundance.txt",
        read_stat = "Collapsed_isoforms/" + sample + ".collapsed.read_stat.txt"

    params:
        prefix = "Collapsed_isoforms/" + sample + ".collapsed"

    log: "Collapsed_isoforms/" + sample + "_get_abundance_all.log"

    shell:
        """
        (get_abundance_post_collapse.py {params.prefix} {input.cluster_report}) 2> {log}
        """



rule get_abundance_demux:
    input:
        mapped_fa = rules.collapse_isoforms.output.mapped_fa,
        read_stat = rules.get_abundance_all.output.read_stat,
        classify_csv = rules.merge_refine_reports.output.merged_report

    output:
        abund = "Collapsed_isoforms/" + sample + ".mapped_fl_count.txt"

    log: "Collapsed_isoforms/" + sample + "_get_abundance_demux.log"

    shell:
        """
        (demux_isoseq_with_genome.py --mapped_fafq {input.mapped_fa} --read_stat {input.read_stat} --classify_csv {input.classify_csv} -o {output.abund}) 2> {log}
        """



rule sqanti_qc:
    input:
        isoforms = rules.collapse_isoforms.output.gff,
        fl_count = rules.get_abundance_demux.output.abund if config.get("primers1") else rules.get_abundance_all.output.abund,
        gtf = config["gtf"],
        genome = config["genome"],
        cage = config["cage"],
        sjs = config["intropolis"],
        poly_peak = config["polya_atlas"],
        poly_motifs = SNAKEDIR + "/data/human.polyA.list.txt"

    output:
        classification = "Sqanti/" + sample + "_classification.txt",
        juncs = "Sqanti/" + sample + "_junctions.txt"

    params:
        res_dir = 'Sqanti/',
        prefix = sample,
        sqanti_dir = config["sqanti_dir"],
        chunks = 2 if gene_ids == "" else 1

    threads: config["threads"]

    log: "Sqanti/" + sample + "_sqanti.log"

    #conda: config["sqanti_dir"] + "/SQANTI3.conda_env.yml"

    shell:
        """
        (PYTHONPATH=$CONDA_PREFIX/bin python {params.sqanti_dir}/sqanti3_qc.py {input.isoforms} {input.gtf} {input.genome} --cage_peak {input.cage} --polyA_peak {input.poly_peak} --polyA_motif_list {input.poly_motifs} -c {input.sjs} -t {threads} --chunks {params.chunks} --output {params.prefix} --dir {params.res_dir} --report skip -fl {input.fl_count}) 2> {log}
        """

