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

# samples

for name, genelist in genes.items():
    for i, dataset in enumerate(datasets):
        SrunCommand("scratch/unaligned/{}/sample.MC{}.fa".format(name, dataset),
                    ["lib/sample.py", Value(nsamples), Value(seeds[i])] + \
                    ["{}/MC{}/{}.codons.txt".format(data_dir, dataset, gene) for gene in genelist],
                    "python $SOURCES $TARGETS".format(name, dataset))

# consensus

for name, genelist in genes.items():
    env.Command("scratch/unaligned/{}/consensus.fa".format(name),
                ["lib/consensus.py"] + ["{}/MC{}/{}.codons.txt".format(data_dir, dataset, gene)
                                        for dataset in datasets
                                        for gene in genelist],
                "python $SOURCES $TARGETS")

# alignments

for name in genes:

    for i in range(nsamples):
        SrunCommand(["scratch/aligned/{}/sample.{}.fa".format(name, i),
                     "scratch/aligned/{}/sample.{}.fa.log".format(name, i)],
                    ["lib/select-sample.py", Value(i)] + \
                    ["scratch/unaligned/{}/sample.MC{}.fa".format(name, dataset) for dataset in datasets],
                    "python $SOURCES | mafft --op 2 --thread $CPUS --auto - 1> ${TARGETS[0]} 2> ${TARGETS[1]}",
                    wrap=True)

    SrunCommand(["scratch/aligned/{}/consensus.fa".format(name),
                 "scratch/aligned/{}/consensus.fa.log".format(name)],
                "scratch/unaligned/{}/consensus.fa".format(name),
                "mafft --op 2 --thread $CPUS --auto $SOURCES 1> ${TARGETS[0]} 2> ${TARGETS[1]}",
                wrap=True)

# trees

#for gene in ["prrt", "int", "wgs"]:
#    SrunCommand(["scratch/trees/{}/RAxML_info.consensus".format(gene),
#                 "scratch/trees/{}/RAxML_bestTree.consensus".format(gene),
#                 "scratch/trees/{}/RAxML_bipartitions.consensus".format(gene),
#                 "scratch/trees/{}/RAxML_bipartitionsBranchLabels.consensus".format(gene),
#                 "scratch/trees/{}/RAxML_bootstrap.consensus".format(gene)],
#                ["lib/raxml.sh",
#                 "scratch/aligned/{}/consensus.fa".format(gene),
#                 Value(100),
#                 Value("GTRGAMMA")],
#                "bash $SOURCES $CPUS scratch/trees/{0} consensus &> logs/raxml/{0}/consensus.log".format(gene),
#                cpus=20)
#    for i in range(nsamples):
#        SrunCommand(["scratch/trees/{}/RAxML_info.sample.{}".format(gene, i),
#                     "scratch/trees/{}/RAxML_bestTree.sample.{}".format(gene, i),
#                     "scratch/trees/{}/RAxML_bipartitions.sample.{}".format(gene, i),
#                     "scratch/trees/{}/RAxML_bipartitionsBranchLabels.sample.{}".format(gene, i),
#                     "scratch/trees/{}/RAxML_bootstrap.sample.{}".format(gene, i)],
#                    ["lib/raxml.sh",
#                     "scratch/aligned/{}/sample.{}.fa".format(gene, i),
#                     Value(100),
#                     Value("GTRGAMMA")],
#                    "bash $SOURCES $CPUS scratch/trees/{0} sample.{1} &> logs/raxml/{0}/sample.{1}.log".format(gene, i),
#                    cpus=20)

# vim: syntax=python expandtab sw=4 ts=4
