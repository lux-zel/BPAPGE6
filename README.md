#dependencies:
* Biopython
* psycopg2
* NCBI+     (on path)
* BLAST+    (on path)

Run the search from `src`:

```sh
python3 main.py seq.fasta proteome.fasta
```

UniProt records and any referenced KEGG entries and ENA nucleotide records are
downloaded and parsed automatically. ENA records are genomic contigs, not
necessarily gene-only sequences.
