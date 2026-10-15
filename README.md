#dependencies:
* Biopython
* psycopg2
* NCBI+     (on path)
* BLAST+    (on path)

Run the application from `src`. Each run replaces the application's database
tables with empty tables, permanently deleting data from previous runs. Other
database tables are left untouched. The application prompts for the database
password once.
Connect to the required VPN before running it.

```sh
python3 main.py seq.fasta proteome.fasta
```

<<<<<<< HEAD
`database_opzetten.py` creates the tables and drops any existing project tables
first, so run it only when initializing or intentionally resetting the database.
The `Gen` table stores the NCBI Gene ID, gene name, organism, and Ensembl IDs
available from UniProt records. Gene nucleotide sequence, genomic coordinates,
strand, and gene-specific description are not supplied by the current inputs,
so those columns remain empty. Protein, GO, and KEGG annotations are stored in
their respective tables.
=======
`main.py` stores one BLAST HSP per protein in the existing `Eiwit` table,
choosing the HSP with the highest bit score (and lowest e-value to break ties).
Transcript isoforms are stored separately in `Isomeer` when ENA
provides annotated mRNA features. The stored transcript sequence is spliced
according to the annotated exon locations. Coding status is based on whether a
CDS feature exists for the same locus tag. This captures annotated transcripts
in the downloaded records, not necessarily every biological splice variant.
>>>>>>> 8e189f3616bbdb02b9f56e4dfcf851b284346df0

UniProt gene names and locus tags are stored in `Gen`. Annotated ENA EMBL
records are matched using their protein IDs to add each gene's nucleotide
sequence, coordinates, strand, and product description when available. ENA
records are genomic contigs, not necessarily gene-only sequences.
