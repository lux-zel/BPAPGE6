#!/bin/sh
# Stop direct als er een commando faalt
set -e

# Zorg dat we in de juiste map werken waar de bestanden staan
cd "$(dirname "$0")"

echo "Start makeblastdb met subject: $SUBJECT"
makeblastdb -in "$SUBJECT" -dbtype prot -out proteome_db
echo 'Database gemaakt van het proteoom'

echo "\nStart blastx met query: $QUERY"
blastx -query "$QUERY" -db proteome_db -outfmt 5 -evalue 1e-5 -out hits.xml
echo '\nBLAST uitgevoerd en resultaten opgeslagen in hits.xml'
