#dependencies:
* Biopython
* psycopg2
* NCBI+     (on path)
* BLAST+    (on path)

Run the application from `src`. It will initialize any missing database tables
without deleting existing data, and prompt for the database password once.
Connect to the required VPN before running it.

```sh
python3 main.py seq.fasta proteome.fasta
```

`main.py` stores each complete BLAST hit, including its alignment and parsed
annotations, in a JSONB list on the matched protein. There is no BLAST HSP
table. HSPs describe aligned regions, not transcript isoforms, so the
application does not create isomer rows.

UniProt records and any referenced KEGG entries and ENA nucleotide records are
downloaded and parsed automatically. ENA records are genomic contigs, not
necessarily gene-only sequences.
