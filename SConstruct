import os

env = Environment(ENV=os.environ)
env.CacheDir("cache")
env.Decider("MD5-timestamp")

datasets = [13, 17, 19, 21, 23, 24, 25, 26, 27, 28, 29, 30, 32, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 58, 59, 60]
genes = ["pol", "gag", "env"]

min_coverage = 10
min_freq = 0.01

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
    env.Command("scratch/{}.fa".format(gene),
                ["lib/clean-alignment.py", "input/HIV1_FLT_2017_{}_PRO.fasta".format(gene)],
                "python $SOURCES > $TARGET 2> logs/clean-alignment-{}.log".format(gene))
    SrunCommand("scratch/{}.hmm".format(gene),
                "scratch/{}.fa".format(gene),
                "hmmbuild $TARGET $SOURCE &> logs/hmmbuild-{}.log".format(gene))

# hivmmer

for gene in genes:
    for dataset in datasets:
        SrunCommand(["scratch/MC{}.{}.hmmsearch1.codons.csv".format(dataset, gene),
                     "scratch/MC{}.{}.hmmsearch1.aavf".format(dataset, gene),
                     "scratch/MC{}.{}.hmmsearch2.codons.csv".format(dataset, gene),
                     "scratch/MC{}.{}.hmmsearch2.aavf".format(dataset, gene),
                     "scratch/MC{}.{}.pear.assembled.fastq".format(dataset, gene),
                     "scratch/MC{}.{}.pear.unassembled.forward.fastq".format(dataset, gene),
                     "scratch/MC{}.{}.pear.unassembled.reverse.fastq".format(dataset, gene),
                     "scratch/MC{}.{}.pear.discarded.fastq".format(dataset, gene)],
                    [Value("--id"), Value("scratch/MC{}.{}".format(dataset, gene)),
                     Value("--fq1"), "input/MC{}_1.fastq.gz".format(dataset),
                     Value("--fq2"), "input/MC{}_2.fastq.gz".format(dataset),
                     Value("--ref"), "scratch/{}.hmm".format(gene)],
                    "hivmmer --cpu $CPUS $SOURCES &> logs/hivmmer-MC{}-{}.log".format(dataset, gene),
                    cpus=4)

# consensus

for gene in genes:
    env.Command(["scratch/unaligned/{}/consensus.fa".format(gene),
                 "scratch/unaligned/{}/consensus.pfa".format(gene)],
                ["lib/consensus.py", Value(min_coverage), Value(min_freq)] + \
                ["scratch/MC{}.{}.hmmsearch2.codons.csv".format(dataset, gene) for dataset in datasets],
                "python $SOURCES $TARGETS &> logs/consensus-{}.log".format(gene))

# vim: syntax=python expandtab sw=4 ts=4
