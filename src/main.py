import sys
from blast import blast_pairwise
from database_opzetten import database_opzetten
from parser import download_uniprot, parse_uniprot

def main():

    if len(sys.argv) < 3:
        print("roep het programma met: python3 main.py <QUERY.fasta> <SUBJECT.fasta>")
        sys.exit(1)
    
    query = sys.argv[1]
    subject = sys.argv[2]

    blast_pairwise(query, subject)
    database_opzetten()


if __name__ == '__main__':
    main()
