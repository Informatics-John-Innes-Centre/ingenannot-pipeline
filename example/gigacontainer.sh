if [ ! -f ../gigacontainer.sif ];
then
    echo "The container image is not present! In the root of the repository, run 'apptainer build gigacontainer.sif gigacontainer.def'."
    exit 1
fi

if ! type apptainer > /dev/null 2>&1
then
    echo "You do not have apptainer on your \$PATH. Install it and try again!"
    exit 1
fi

export NXF_OFFLINE='true'
apptainer run ../gigacontainer.sif -with-trace -with-report report.html \
    --genomes genomes \
    --csv input.csv \
    --proteinDatabase Viridiplantae_alt_headers.fa \
    --rnaseq rnaseq \
    --isoseq isoseq \
    --tiberiusModel angiosperms_softmasking \
    --helixerLineage land_plant \
    --annevoModel ANNEVO_Magnoliopsida.pt \
    --annevoLineage Magnoliopsida