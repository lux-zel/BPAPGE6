import sys
from blast import blast_pairwise
from database_opzetten import database_opzetten

def main():

    if len(sys.argv) < 3:
        print("Roep het programma aan met: python3 main.py <QUERY.fasta> <SUBJECT.fasta>")
        sys.exit(1)
    
    query = sys.argv[1]
    subject = sys.argv[2]

    hits = blast_pairwise(query, subject)
    database_opzetten()

    for hit in hits:
        print("E-value:", hit["evalue"])
        print("Bit score:", hit["bit-score"])
        print("Identities:", hit["identities"])
        print("Alignment length:", hit["align_length"])
        print("Gaps:", hit["gaps"])
        print("Query start:", hit["query_start"])
        print("Query end:", hit["query_end"])
        print("Subject start:", hit["subject_start"])
        print("Subject end:", hit["subject_end"])
        print("Match:", hit["match"])
        print("Subject sequence:", hit["subject_seq"])
        print()


if __name__ == '__main__':
    main()
