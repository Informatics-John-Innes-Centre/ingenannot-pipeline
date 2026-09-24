if [ ! -d ../apptainer-cache ];
then
    echo "You must run the script 'containers.sh' in the root of the repository to build the container cache first!"
    exit 1
fi

if ! type nextflow > /dev/null 2>&1
then
    echo "You do not have nextflow on your \$PATH. Install it and try again!"
    exit 1
fi

if ! type apptainer > /dev/null 2>&1
then
    echo "You do not have apptainer on your \$PATH. Install it and try again!"
    exit 1
fi

export NXF_OFFLINE='true'
nextflow run ../main.nf -resume -with-trace -profile local -with-report \
    --genomes genomes \
    --csv input.csv \
    --proteinDatabase Viridiplantae_alt_headers.fa \
    --rnaseq rnaseq \
    --isoseq isoseq \
    --tiberiusModel angiosperms_softmasking \
    --helixerLineage land_plant \
    --annevoModel ANNEVO_Magnoliopsida.pt \
    --annevoLineage Magnoliopsida