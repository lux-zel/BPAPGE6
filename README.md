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

`database_opzetten.py` creates the tables and drops any existing project tables
first, so run this script only when initializing or intentionally resetting the 
database. The `Gen` table stores the NCBI Gene ID, gene name, organism, and 
Ensembl IDs available from UniProt records. Gene nucleotide sequence, genomic 
coordinates, strand, and gene-specific description are not supplied by the current 
inputs,so those columns remain empty. Protein, GO, and KEGG annotations are stored 
in their respective tables.

UniProt gene names and locus tags are stored in `Gen`. Annotated ENA EMBL
records are matched using their protein IDs to add each gene's nucleotide
sequence, coordinates, strand, and product description when available. ENA
records are genomic contigs, not necessarily gene-only sequences.
