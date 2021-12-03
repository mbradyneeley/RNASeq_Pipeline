#!/usr/bin/env bash

#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -o pipelineKickOff-%j.out
#SBATCH -e pipelineKickOff-%j.err
#SBATCH --mail-user=brady.neeley@hsc.utah.edu
#SBATCH --mail-type=END
#SBATCH --account=pezzolesi
#SBATCH --partition=notchpeak
#SBATCH --mem=10G

# this vvvv should be either 'resume' or 'new'
resume=$1
configFile=$2
# rename this vvvv to create a new directory for your project (i.e. if you want to start over without deleting what you've already done)
scratchDir="/scratch/general/lustre/$USER/theBifrost"

if [[ $resume == "resume" ]]; then
    if [ -d $scratchDir ]; then
        echo -e "\nResuming your previous run\n"
        cp ./main.nf $configFile $scratchDir
        #cp -r ./bin $scratchDir
        cd $scratchDir
        nextflow -C $configFile run -with-report -with-trace -with-timeline -with-dag dag.html main.nf -resume
    else
        echo "There's nothing to resume (scratch directory doesn't exist)"
    fi
elif [[ $resume == "new" ]]; then
    if [ ! -d $scratchDir ]; then
        echo -e "\nStarting a new project\n"
        mkdir -p $scratchDir/{results/{fastp,fastqc,bqsr,vqsr,bam/{stats,coverage},gvcf,vcf/stats},bin}
        cp ./main.nf $configFile $scratchDir
        cp -r ./bin $scratchDir
        cd $scratchDir
        nextflow -C $configFile  run -with-report -with-trace -with-timeline -with-dag dag.html main.nf
    else
        echo "A project already exists. Delete it and start over or resume it"
    fi
else
    echo "Takes string argument 'resume' or 'new'"
fi
