#!/bin/bash

makeblastdb -in "$SUBJECT" -dbtype prot -out proteome_db
echo 'Database gemaakt van het proteoom'

blastx -query "$QUERY" -db proteome_db -outfmt 5 -out hits.xml
echo 'BLAST uitgevoerd en resultaten opgeslagen in hits.xml'