#dependencies:
* Biopython
* psycopg2
* NCBI+     (on path)
* BLAST+    (on path)

Run the search from `src`:

```sh
export DATABASE_URL='postgresql://USER:PASSWORD@HOST:5432/DATABASE'
python3 database_opzetten.py
python3 main.py seq.fasta proteome.fasta
```

`database_opzetten.py` creates the tables and drops any existing project tables
first, so run it only when initializing or intentionally resetting the database.
The `Gen` table stores the NCBI Gene ID, gene name, organism, and Ensembl IDs
available from UniProt records. Gene nucleotide sequence, genomic coordinates,
strand, and gene-specific description are not supplied by the current inputs,
so those columns remain empty. Protein, GO, and KEGG annotations are stored in
their respective tables.

UniProt records and any referenced KEGG entries and ENA nucleotide records are
downloaded and parsed automatically. ENA records are genomic contigs, not
necessarily gene-only sequences.
