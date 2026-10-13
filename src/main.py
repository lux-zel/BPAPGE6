import sys
from blast import blast_pairwise
from parser import enrich_hits
from database_vullen import fill_database
from database_connectie import maak_database_connectie
from database_opzetten import database_opzetten

def main():
    if len(sys.argv) < 3:
        print("Roep het programma aan met: python3 main.py <QUERY.fasta> <SUBJECT.fasta>")
        sys.exit(1)
    
    query = sys.argv[1]
    subject = sys.argv[2]

    conn = maak_database_connectie()
    try:
        database_opzetten(conn)
        print(f"Starten van BLAST pairwise met query: {query} en subject: {subject}")
        data = enrich_hits(blast_pairwise(query, subject))
        opgeslagen = fill_database(data, conn)
    finally:
        conn.close()

    print(f"Finished: {opgeslagen} HSP's verwerkt en opgeslagen")


if __name__ == '__main__':
    main()
