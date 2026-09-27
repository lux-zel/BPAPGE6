import sys
from blast import blast_pairwise
from parser import download_uniprot, parse_uniprot

def main():

    if len(sys.argv) < 3:
        print("Roep het programma aan met: python3 main.py <QUERY.fasta> <SUBJECT.fasta>")
        sys.exit(1)
    
    query = sys.argv[1]
    subject = sys.argv[2]

    print(f"Starten van BLAST pairwise met query: {query} en subject: {subject}")
    hits = blast_pairwise(query, subject)

    for hit in hits:
        print(f"\nUniprot begin download van: {hit['accession']}")
        file = download_uniprot(hit['accession'])

        print(f"Parsing van bestand: {file}")
        parse_uniprot(file)


if __name__ == '__main__':
    main()
