import os

env = Environment(ENV=os.environ)
env.CacheDir("cache")
env.Decider("MD5-timestamp")

datasets = [13, 17, 19, 21, 23, 24, 25, 26, 27, 28, 29, 30, 32, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 58, 59, 60]
genes = ["pol", "gag", "env"]

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
                "python $SOURCES > $TARGET")
    SrunCommand("scratch/{}.hmm".format(gene),
                "scratch/{}.fa".format(gene),
                "hmmbuild $TARGET $SOURCE")

# hivmmer

for gene in genes:
    for i in datasets:
        pass

        ## hmmsearch round #1
        #SrunCommand("scratch/MC{}.hmmsearch1.hmm".format(i),
        #            ["lib/hmmsearch.sh",
        #             env.Value(20),
        #             "lib/revcomp_reads.py",
        #             "scratch/MC{}.collapsed.fa".format(i),
        #             "scratch/HIV1_ALL_2016_genome_DNA.hmm"],
        #            "$SOURCES $TARGET",
        #            cpus=20, mem_per_cpu=6)


# vim: syntax=python expandtab sw=4 ts=4
