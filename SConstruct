import numpy as np
import os

env = Environment(ENV=os.environ)
env.CacheDir("cache")
env.Decider("MD5-timestamp")

datasets = [14, 17, 20, 21, 23, 24, 25, 26, 27, 28, 29, 31, 32, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 58, 59, 60]
genes = ["pol", "gag", "env"]

min_coverage = 10
min_freq = 0.01
nsamples = 500

np.random.seed(919047801)
seeds = [np.random.randint(1000000000) for _ in range(nsamples)]

# If your SLURM culster requires additional parameters, you can include them in
# this command.
srun = None
srun = "srun"

def SrunCommand(targets, sources, cmd, wrap=False, prefix="", cpus=1, mem_per_cpu=3, timelimit="24:00:00"):
    global srun, env
    if wrap:
        cmd = "sh -c '{}'".format(cmd)
    if srun is not None:
        cmd = "{} -c {} --mem={}G -t {} {}".format(srun, cpus, cpus*mem_per_cpu, timelimit,
                                                   cmd.replace("$CPUS", str(cpus)))
    if prefix:
        cmd = "{} {}".format(prefix, cmd)
    return env.Command(targets, sources, cmd)

# Reference data

for gene in genes:
    env.Command("scratch/reference/{}.fa".format(gene),
                ["lib/clean-alignment.py", "input/HIV1_FLT_2017_{}_PRO.fasta".format(gene)],
                "python $SOURCES > $TARGET 2> logs/clean-alignment/{}.log".format(gene))
    SrunCommand("scratch/reference/{}.hmm".format(gene),
                "scratch/reference/{}.fa".format(gene),
                "hmmbuild $TARGET $SOURCE &> logs/hmmbuild/{}.log".format(gene))

# hivmmer

for gene in genes:
    for dataset in datasets:
        SrunCommand(["scratch/hivmmer/{}/MC{}.hmmsearch1.codons.csv".format(gene, dataset),
                     "scratch/hivmmer/{}/MC{}.hmmsearch1.codons.txt".format(gene, dataset),
                     "scratch/hivmmer/{}/MC{}.hmmsearch1.aavf".format(gene, dataset),
                     "scratch/hivmmer/{}/MC{}.hmmsearch2.codons.csv".format(gene, dataset),
                     "scratch/hivmmer/{}/MC{}.hmmsearch2.codons.txt".format(gene, dataset),
                     "scratch/hivmmer/{}/MC{}.hmmsearch2.aavf".format(gene, dataset),
                     "scratch/hivmmer/{}/MC{}.pear.assembled.fastq".format(gene, dataset),
                     "scratch/hivmmer/{}/MC{}.pear.unassembled.forward.fastq".format(gene, dataset),
                     "scratch/hivmmer/{}/MC{}.pear.unassembled.reverse.fastq".format(gene, dataset),
                     "scratch/hivmmer/{}/MC{}.pear.discarded.fastq".format(gene, dataset)],
                    [Value("--id"), Value("scratch/hivmmer/{}/MC{}".format(gene, dataset)),
                     Value("--fq1"), "input/MC{}_1.fastq.gz".format(dataset),
                     Value("--fq2"), "input/MC{}_2.fastq.gz".format(dataset),
                     Value("--ref"), "scratch/reference/{}.hmm".format(gene)],
                    "hivmmer --cpu $CPUS $SOURCES &> logs/hivmmer/{}/MC{}.log".format(gene, dataset),
                    cpus=4)

# consensus

for gene in genes:
    SrunCommand(["scratch/unaligned/{}/consensus.fa".format(gene),
                 "scratch/unaligned/{}/consensus.pfa".format(gene)],
                ["lib/consensus.py", Value(min_coverage)] + \
                ["scratch/hivmmer/{}/MC{}.hmmsearch2.codons.txt".format(gene, dataset) for dataset in datasets],
                "python $SOURCES $TARGETS &> logs/consensus/{}.log".format(gene))
    SrunCommand("scratch/unaligned/{}/consensus.hmmsearch.txt".format(gene),
                ["scratch/reference/{}.hmm".format(gene),
                 "scratch/unaligned/{}/consensus.pfa".format(gene)],
                "hmmsearch --notextw -o $TARGETS $SOURCES &> logs/hmmsearch/{}/consensus.log".format(gene))
    SrunCommand("scratch/aligned/{}/consensus.fa".format(gene),
                ["lib/codon-align.py",
                 "scratch/unaligned/{}/consensus.fa".format(gene),
                 "scratch/unaligned/{}/consensus.hmmsearch.txt".format(gene)],
                "python $SOURCES > $TARGETS 2> logs/codon-align/{}/consensus.log".format(gene))

# samples

for gene in genes:
    for i, seed in enumerate(seeds):
        SrunCommand(["scratch/unaligned/{}/sample.{}.fa".format(gene, i),
                     "scratch/unaligned/{}/sample.{}.pfa".format(gene, i)],
                    ["lib/sample.py", Value(seed), Value(min_coverage), Value(min_freq)] + \
                    ["scratch/hivmmer/{}/MC{}.hmmsearch2.codons.txt".format(gene, dataset) for dataset in datasets],
                    "python $SOURCES $TARGETS &> logs/sample/{}/{}.log".format(gene, i))
        SrunCommand("scratch/unaligned/{}/sample.{}.hmmsearch.txt".format(gene, i),
                    ["scratch/reference/{}.hmm".format(gene),
                     "scratch/unaligned/{}/sample.{}.pfa".format(gene, i)],
                    "hmmsearch --notextw -o $TARGETS $SOURCES &> logs/hmmsearch/{}/sample.{}log".format(gene, i))
        SrunCommand("scratch/aligned/{}/sample.{}.fa".format(gene, i),
                    ["lib/codon-align.py",
                     "scratch/unaligned/{}/sample.{}.fa".format(gene, i),
                     "scratch/unaligned/{}/sample.{}.hmmsearch.txt".format(gene, i)],
                    "python $SOURCES > $TARGETS 2> logs/codon-align/{}/sample.{}.log".format(gene, i))

# concatenated alignments

env.Command("scratch/aligned/gag-pol-env/consensus.fa",
            ["lib/concat.py",
             "scratch/aligned/gag/consensus.fa",
             "scratch/aligned/pol/consensus.fa",
             "scratch/aligned/env/consensus.fa"],
            "python $SOURCES $TARGET &> logs/concat/gag-pol-env/consensus.log")

for i in range(nsamples):
    env.Command("scratch/aligned/gag-pol-env/sample.{}.fa".format(i),
                ["lib/concat.py",
                 "scratch/aligned/gag/sample.{}.fa".format(i),
                 "scratch/aligned/pol/sample.{}.fa".format(i),
                 "scratch/aligned/env/sample.{}.fa".format(i)],
                "python $SOURCES $TARGET &> logs/concat/gag-pol-env/sample.{}.log".format(i))

# trees

for gene in ["pol", "gag-pol-env"]:
    SrunCommand(["scratch/trees/{}/RAxML_info.consensus".format(gene),
                 "scratch/trees/{}/RAxML_bestTree.consensus".format(gene),
                 "scratch/trees/{}/RAxML_bipartitions.consensus".format(gene),
                 "scratch/trees/{}/RAxML_bipartitionsBranchLabels.consensus".format(gene),
                 "scratch/trees/{}/RAxML_bootstrap.consensus".format(gene)],
                ["lib/raxml.sh",
                 "scratch/aligned/{}/consensus.fa".format(gene),
                 Value(100),
                 Value("GTRGAMMA")],
                "bash $SOURCES $CPUS scratch/trees/{0} consensus &> logs/raxml/{0}/consensus.log".format(gene),
                cpus=20)
    for i in range(nsamples):
        SrunCommand(["scratch/trees/{}/RAxML_info.sample.{}".format(gene, i),
                     "scratch/trees/{}/RAxML_bestTree.sample.{}".format(gene, i),
                     "scratch/trees/{}/RAxML_bipartitions.sample.{}".format(gene, i),
                     "scratch/trees/{}/RAxML_bipartitionsBranchLabels.sample.{}".format(gene, i),
                     "scratch/trees/{}/RAxML_bootstrap.sample.{}".format(gene, i)],
                    ["lib/raxml.sh",
                     "scratch/aligned/{}/sample.{}.fa".format(gene, i),
                     Value(100),
                     Value("GTRGAMMA")],
                    "bash $SOURCES $CPUS scratch/trees/{0} sample.{1} &> logs/raxml/{0}/sample.{1}.log".format(gene, i),
                    cpus=20)

# vim: syntax=python expandtab sw=4 ts=4
