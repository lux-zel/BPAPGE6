# dependencies:
* Biopython
* psycopg2
* requests
* NCBI+     (on path)
* BLAST+    (on path)
* Internet access to UniProt and KEGG APIs

Om het programma te kunnen gebruiken moet je verbonden zijn met eduvpn hogeschool leiden. Run het programma vanuit `/src`. met het commando:

```sh
python3 main.py seq.fasta proteome.fasta
```

Isoformnamen uit de headers van het opgegeven proteoombestand worden samen met de UniProt-accession in `Isomeer` opgeslagen.
