#!/bin/sh
set -e

: "${QUERY:?QUERY not set}" "${SUBJECT:?SUBJECT not set}" "${OUT:?OUT not set}"

cd "$(dirname "$0")"

printf 'Start makeblastdb met subject: %s\n' "$SUBJECT"
makeblastdb -in "$SUBJECT" -dbtype prot -parse_seqids -out proteome_db
printf 'Database gemaakt van het proteoom\n'

printf '\nStart blastx met query: %s\n' "$QUERY"
blastx -query "$QUERY" -db proteome_db -outfmt 5 -evalue 1e-5 -out "$OUT"
printf 'BLAST query opgeslagen: %s\n' "$OUT"