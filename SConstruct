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

    for dataset in datasets:
        SrunCommand(["scratch/aligned/{}/sample.MC{}.aa.fa".format(gene, dataset),
                     "scratch/aligned/{}/sample.MC{}.fa".format(gene, dataset),
                     "scratch/aligned/{}/sample.MC{}.fa.log".format(gene, dataset)],
                    ["lib/codon-align.py",
                     "scratch/unaligned/{}/sample.MC{}.fa".format(gene, dataset)],
                    "python $SOURCES ${TARGETS[0]} ${TARGETS[1]} 2> ${TARGETS[2]}")

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

for dataset in datasets:
    env.Command(["scratch/aligned/wgs/sample.MC{}.fa".format(dataset)],
                ["lib/concat.py"] + \
                ["scratch/aligned/{}/sample.MC{}.fa".format(gene, dataset) for gene in genes["wgs"]],
                "python $SOURCES > $TARGET")

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

genes = ["prrt", "int", "env", "wgs"]

for gene in genes:

    for dataset in datasets:
        SrunCommand(["scratch/trees/{}/sample.MC{}.log".format(gene, dataset),
                     "scratch/trees/{}/RAxML_info.sample.MC{}".format(gene, dataset),
                     "scratch/trees/{}/RAxML_bestTree.sample.MC{}".format(gene, dataset),
                     "scratch/trees/{}/RAxML_bipartitions.sample.MC{}".format(gene, dataset),
                     "scratch/trees/{}/RAxML_bipartitionsBranchLabels.sample.MC{}".format(gene, dataset),
                     "scratch/trees/{}/RAxML_bootstrap.sample.MC{}".format(gene, dataset)],
                    ["lib/raxml.sh",
                     "scratch/aligned/{}/sample.MC{}.fa".format(gene, dataset),
                     Value(100),
                     Value("GTRGAMMA")],
                    "bash $SOURCES $CPUS scratch/trees/{0} sample.MC{1} > $TARGET".format(gene, dataset),
                    cpus=20)

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

    env.Command(["scratch/trees/{}/RAxML_bestTree.samples".format(gene)],
                ["scratch/trees/{}/RAxML_bestTree.sample.{}".format(gene, i) for i in range(nsamples)],
                "cat $SOURCES > $TARGET")

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

# compute all pairwise intra-patient genetic distances

for gene in genes:
    for dataset in datasets:
        SrunCommand(["scratch/aligned/{}/distances.MC{}.csv".format(gene, dataset)],
                    ["lib/genetic-distance.py",
                     "scratch/aligned/{}/sample.MC{}.fa".format(gene, dataset)],
                    "python $SOURCES $TARGET")

# compute all pairwise tree distances

SrunCommand(["scratch/trees/distance.RData"],
            ["lib/tree-distance.R"] + \
            ["scratch/trees/{}/RAxML_bestTree.consensus".format(gene) for gene in genes] + \
            ["scratch/trees/{}/RAxML_bestTree.samples".format(gene) for gene in genes],
            "Rscript $SOURCES $TARGET")

SrunCommand(["scratch/trees/mds.RData"],
            ["lib/tree-mds.R",
             "scratch/trees/distance.RData"],
            "Rscript $SOURCES $TARGET")

# compute pairwise tree distances within genes

for gene in genes:

    SrunCommand(["scratch/trees/distance.{}.RData".format(gene)],
                ["lib/tree-distance.R",
                 "scratch/trees/{}/RAxML_bestTree.consensus".format(gene),
                 "scratch/trees/{}/RAxML_bestTree.samples".format(gene)],
                "Rscript $SOURCES $TARGET")

    SrunCommand(["scratch/trees/mds.{}.RData".format(gene)],
                ["lib/tree-mds.R",
                 "scratch/trees/distance.{}.RData".format(gene)],
                "Rscript $SOURCES $TARGET")

# clusters

names = ["consensus"] + ["sample.{}".format(i) for i in range(nsamples)]
for name in names:
   for gene in genes:
        # copy local files for ClusterPicker
        env.Command("scratch/clusters/{}/{}.fa".format(gene, name),
                    "scratch/aligned/{}/{}.fa".format(gene, name),
                    "cp $SOURCE $TARGET")
        env.Command("scratch/clusters/{}/{}.nwk".format(gene, name),
                    "scratch/trees/{}/RAxML_bipartitions.{}".format(gene, name),
                    "cp $SOURCE $TARGET")
        # run ClusterPicker
        env.Command(["scratch/clusters/{}/{}.log".format(gene, name),
                     "scratch/clusters/{}/{}_clusterPicks.nwk".format(gene, name),
                     "scratch/clusters/{}/{}_clusterPicks.nwk.figTree".format(gene, name),
                     "scratch/clusters/{}/{}_clusterPicks_log.txt".format(gene, name),
                     "scratch/clusters/{0}/{1}.fa_{1}_clusterPicks.fas".format(gene, name)],
                    ["lib/ClusterPicker_1.2.3.jar",
                     "scratch/clusters/{}/{}.fa".format(gene, name),
                     "scratch/clusters/{}/{}.nwk".format(gene, name),
                     Value(80), Value(80), Value(0.045), Value(0), Value("ambiguity")],
                    "java -jar $SOURCES > $TARGET")

for gene in genes:
    env.Command(["scratch/clusters/{}/support.csv".format(gene)],
                ["lib/cluster-support.py",
                 "scratch/clusters/{}/consensus.fa_consensus_clusterPicks.fas".format(gene)] + \
                ["scratch/clusters/{0}/sample.{1}.fa_sample.{1}_clusterPicks.fas".format(gene, i)
                 for i in range(nsamples)],
                "python $SOURCES $TARGETS")

# figures

env.Command(["manuscript/Figure1.pdf"],
            ["lib/Figure1.R"] + \
            ["scratch/aligned/{}/distances.MC{}.csv".format(gene, dataset) for gene in genes for dataset in datasets],
            "Rscript $SOURCES $TARGET")

env.Command(["manuscript/Figure2.pdf",
             "manuscript/Figure2.log"],
            ["lib/Figure2.R"] + \
            ["scratch/trees/mds.{}.RData".format(gene) for gene in genes],
            "Rscript $SOURCES ${TARGETS[0]} > ${TARGETS[1]}")

env.Command(["manuscript/Figure3.pdf"],
            ["lib/Figure3.R"] + \
            ["scratch/trees/{}/RAxML_bestTree.consensus".format(gene) for gene in genes] + \
            ["scratch/trees/{}/RAxML_bestTree.samples".format(gene) for gene in genes],
            "Rscript $SOURCES $TARGET")

env.Command(["manuscript/Figure4.pdf"],
            ["lib/Figure4.R"] + \
            list(zip(["scratch/trees/{}/RAxML_bestTree.consensus".format(gene) for gene in genes],
                     ["scratch/clusters/{}/support.csv".format(gene) for gene in genes])),
            "Rscript $SOURCES $TARGET")

env.Command(["manuscript/FigureS1.pdf",
             "manuscript/FigureS1.log"],
            ["lib/FigureS1.R",
             "scratch/trees/mds.RData"],
            "Rscript $SOURCES ${TARGETS[0]} > ${TARGETS[1]}")

env.Command(["manuscript/FigureS2.pdf"],
            ["lib/FigureS2.R"] + \
            ["scratch/trees/{}/RAxML_bootstrap.consensus".format(gene) for gene in genes],
            "Rscript $SOURCES $TARGET")

env.Command(["manuscript/FigureS3.pdf"],
            ["lib/FigureS2.R"] + \
            ["scratch/trees/{}/RAxML_bestTree.samples".format(gene) for gene in genes],
            "Rscript $SOURCES $TARGET")

# vim: syntax=python expandtab sw=4 ts=4
