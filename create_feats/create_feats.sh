# Script to generate feature maps of a set of surrogates using cyclesmap fortran script.
#
# Parameters:
#
# - surrogates folder: name of the folder containing surrogates (they should be organized by dimension, as obtained with create_SD.sh script)
# - results folder: name of the destination folder to store the results

mkdir -p $2

for dir in $1/*/     
do 
    dir=${dir%*/}  
    d=${dir: -1}
    mkdir -p $2/D$d
    echo $2/D$d

    for filename in $dir/*.edge; do
        name=$(basename -- $filename)
        name=${name%".edge"}
        ./task_feats.sh $filename $2/D$d/$name".csv"
   done
done

