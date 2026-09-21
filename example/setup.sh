#! /bin/bash -e

if ! type apptainer > /dev/null 2>&1
then
    echo "You do not have apptainer on your path. Install apptainer before proceeding."
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
    echo "Building container cache... this may take a while!"
    (cd ../; ./containers.sh)
elif [[ $executor == "l" ]]; then
    echo "Building pipeline apptainer image... this may take a while!"
    (cd ../; apptainer build gigacontainer.sif gigacontainer.def)
fi
    
