#!/bin/sh
set -e

: "${QUERY:?QUERY not set}" "${SUBJECT:?SUBJECT not set}" "${OUT:?OUT not set}"

cd "$(dirname "$0")"

echo "Start makeblastdb met subject: $SUBJECT"
makeblastdb -in "$SUBJECT" -dbtype prot -parse_seqids -out proteome_db
echo "Database gemaakt van het proteoom"

printf '\nStart blastx met query: %s\n' "$QUERY"
blastx -query "$QUERY" -db proteome_db -outfmt 5 -evalue 1e-5 -out "$OUT"
echo "BLAST resultaten opgeslagen in $OUT"