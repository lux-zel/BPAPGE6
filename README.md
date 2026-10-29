# dependencies:
* Biopython
* psycopg2
* requests
* NCBI-blast+     (on path)

Om het programma te kunnen gebruiken moet je verbonden zijn met eduvpn hogeschool leiden. Run het programma vanuit `/src`. met het commando:

```sh
python3 main.py seq.fasta proteome.fasta
```

Isoformnamen uit de headers van het opgegeven proteoombestand worden samen met de UniProt-accession in `Isomeer` opgeslagen.
Het tweede argument is het proteoombestand dat als BLAST-database wordt gebruikt. BLASTX zoekt met een E-valuegrens van `1e-5` en maximaal vijf hits per query. De hits worden verrijkt met UniProt- en KEGG-annotaties, waaronder GO-termen en pathways, en opgeslagen in PostgreSQL.
