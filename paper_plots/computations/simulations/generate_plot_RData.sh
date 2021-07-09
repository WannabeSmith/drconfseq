#!/bin/zsh

for folder in */; do
    cd $folder
    for R_file in *.R; do
        echo "############### Generating RData from $R_file ###############"
        Rscript $R_file
    done
    cd ..
done

