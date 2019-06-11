import numpy as np
import os

data_dir  = "/gpfs/data/rkantor/hiv-profile-sampling-data/20190521"
cache_dir = "/gpfs/data/rkantor/hiv-profile-sampling-cache"

env = Environment(ENV=os.environ)
env.CacheDir(cache_dir)
env.Decider("MD5-timestamp")

datasets = [14, 17, 20, 21, 23, 24, 25, 26, 27, 28, 29, 31, 32, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 58, 59, 60]

genes = {
    "prrt": ["prrt"],
    "int":  ["int"],
    "wgs":  ["gag", "pol", "vif", "vpr", "tat", "vpu", "env", "nef"]
}

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


for gene in sum(genes.values(), []):

    # samples

    for i, dataset in enumerate(datasets):
        SrunCommand(["scratch/unaligned/{}/sample.MC{}.fa".format(gene, dataset)],
                    ["lib/sample.py",
                     Value(nsamples),
                     Value(seeds[i]),
                     "{}/MC{}/{}.codons.txt".format(data_dir, dataset, gene)],
                    "python $SOURCES $TARGETS".format(gene, dataset))

    for i in range(nsamples):
        env.Command(["scratch/unaligned/{}/sample.{}.fa".format(gene, i)],
                    ["lib/select-sample.py",
                     Value(i)] + \
                    ["scratch/unaligned/{}/sample.MC{}.fa".format(gene, dataset) for dataset in datasets],
                    "python $SOURCES > $TARGET")

    # consensus

    env.Command(["scratch/unaligned/{}/consensus.fa".format(gene)],
                ["lib/consensus.py"] + ["{}/MC{}/{}.codons.txt".format(data_dir, dataset, gene)
                                        for dataset in datasets],
                "python $SOURCES $TARGETS")

    # alignments

    for i in range(nsamples):
        SrunCommand(["scratch/aligned/{}/sample.{}.aa.fa".format(gene, i),
                     "scratch/aligned/{}/sample.{}.fa".format(gene, i),
                     "scratch/aligned/{}/sample.{}.fa.log".format(gene, i)],
                    ["lib/codon-align.py",
                     "scratch/unaligned/{}/sample.{}.fa".format(gene, i)],
                    "python $SOURCES ${TARGETS[0]} ${TARGETS[1]} 2> ${TARGETS[2]}")

    SrunCommand(["scratch/aligned/{}/consensus.aa.fa".format(gene),
                 "scratch/aligned/{}/consensus.fa".format(gene),
                 "scratch/aligned/{}/consensus.fa.log".format(gene)],
                ["lib/codon-align.py",
                 "scratch/unaligned/{}/consensus.fa".format(gene)],
                "python $SOURCES ${TARGETS[0]} ${TARGETS[1]} 2> ${TARGETS[2]}")

# concatenate wgs alignments

for i in range(nsamples):
    env.Command(["scratch/aligned/wgs/sample.{}.fa".format(i)],
                ["lib/concat.py"] + \
                ["scratch/aligned/{}/sample.{}.fa".format(gene, i) for gene in genes["wgs"]],
                "python $SOURCES > $TARGET")

env.Command(["scratch/aligned/wgs/consensus.fa"],
            ["lib/concat.py"] + \
            ["scratch/aligned/{}/consensus.fa".format(gene) for gene in genes["wgs"]],
            "python $SOURCES > $TARGET")

# trees

for gene in ["prrt", "int", "env", "wgs"]:

    for i in range(nsamples):
        SrunCommand(["scratch/trees/{}/sample.{}.log".format(gene, i),
                     "scratch/trees/{}/RAxML_info.sample.{}".format(gene, i),
                     "scratch/trees/{}/RAxML_bestTree.sample.{}".format(gene, i),
                     "scratch/trees/{}/RAxML_bipartitions.sample.{}".format(gene, i),
                     "scratch/trees/{}/RAxML_bipartitionsBranchLabels.sample.{}".format(gene, i),
                     "scratch/trees/{}/RAxML_bootstrap.sample.{}".format(gene, i)],
                    ["lib/raxml.sh",
                     "scratch/aligned/{}/sample.{}.fa".format(gene, i),
                     Value(100),
                     Value("GTRGAMMA")],
                    "bash $SOURCES $CPUS scratch/trees/{0} sample.{1} > $TARGET".format(gene, i),
                    cpus=20)

    SrunCommand(["scratch/trees/{}/consensus.log".format(gene),
                 "scratch/trees/{}/RAxML_info.consensus".format(gene),
                 "scratch/trees/{}/RAxML_bestTree.consensus".format(gene),
                 "scratch/trees/{}/RAxML_bipartitions.consensus".format(gene),
                 "scratch/trees/{}/RAxML_bipartitionsBranchLabels.consensus".format(gene),
                 "scratch/trees/{}/RAxML_bootstrap.consensus".format(gene)],
                ["lib/raxml.sh",
                 "scratch/aligned/{}/consensus.fa".format(gene),
                 Value(100),
                 Value("GTRGAMMA")],
                "bash $SOURCES $CPUS scratch/trees/{0} consensus > $TARGET".format(gene),
                cpus=20)

# vim: syntax=python expandtab sw=4 ts=4
