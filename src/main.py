from blast.py import blast_pairwise
import sys


def main():
    
    if len(sys.argv) < 3:
        print("roep het programma met: python3 main.py <QUERY.fasta> <SUBJECT.fasta>")
        sys.exit(1)
    
    query_file = sys.argv[1]
    subject_file = sys.argv[2]
    
    blast_pairwise(query_file, subject_file, 'nucl')


if __name__ == '__main__':
    main()
