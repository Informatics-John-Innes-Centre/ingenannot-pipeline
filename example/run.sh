#! /bin/bash -e

if ! type apptainer > /dev/null 2>&1
then
    echo "You do not have apptainer on your \$path. Install apptainer before proceeding."
    exit 1
fi

echo -n "Do you want to run the example via slurm or locally? [s/l]: "
while true; do  
    read executor

    if [[ $executor == "s" ]] || [[ $executor == "l" ]]
    then
        break
    fi
    echo -n "Enter 's' for slurm, or 'l' for local: "
done

if [[ $executor == "s" ]]; then
    if ! type nextflow > /dev/null 2>&1
    then
        echo "You do not have nextflow on your \$PATH. Install nextflow before proceeding."
        exit 1
    fi

    if [ ! -d "../apptainer-cache" ]; then
        echo "You do not have the container image cache. Run setup.sh with the 's' option before proceeding."
        exit 1
    fi

    read -p "Enter the name of the slurm queue(s) (comma separated) used for non-GPU tasks: " queue
    read -p "Enter the name of the slurm queue(s) (comma separated) used for GPU tasks: " gpu_queue

    echo "Starting example run..."

    NXF_OFFLINE='true' NXF_SINGULARITY_CACHEDIR='../apptainer-cache' nextflow run ../main.nf -resume -with-trace -profile hpc -with-report report.html \
        --genomes genomes \
        --csv input.csv \
        --proteinDatabase Viridiplantae_alt_headers.fa \
        --rnaseq rnaseq \
        --isoseq isoseq \
        --tiberiusModel angiosperms_softmasking \
        --helixerLineage land_plant \
        --annevoModel ANNEVO_Magnoliopsida.pt \
        --annevoLineage Magnoliopsida \
        --slurm_queue $queue \
        --slurm_gpu_queue $gpu_queue
    
elif [[ $executor == "l" ]]; then
    if [ ! -f "../gigacontainer.sif" ];
    then
        echo "You have not built the pipeline container. Run setup.sh with the 'l' option before proceeding."
        exit 1
    fi

    echo "Starting example run..."

    NXF_OFFLINE='true' NXF_SINGULARITY_CACHEDIR='../apptainer-cache' apptainer run ../gigacontainer.sif -with-trace -with-report report.html \
        --genomes genomes \
        --csv input.csv \
        --proteinDatabase Viridiplantae_alt_headers.fa \
        --rnaseq rnaseq \
        --isoseq isoseq \
        --tiberiusModel angiosperms_softmasking \
        --helixerLineage land_plant \
        --annevoModel ANNEVO_Magnoliopsida.pt \
        --annevoLineage Magnoliopsida
fi