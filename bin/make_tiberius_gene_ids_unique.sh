#! /usr/bin/env bash

awk 'BEGIN { FS = OFS = "\t" }
     /^#/ || NF < 9 { print; next }
     {
         gsub(/gene_id "[^".]*/,       "&_" $1, $9)
         gsub(/transcript_id "[^".]*/, "&_" $1, $9)
         print
     }' ${1}