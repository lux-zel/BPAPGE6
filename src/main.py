import sys
from blast import blast_pairwise
from parser import enrich_hits, build_gn_lookup
from database_vullen import database_vullen
from database_connectie import maak_database_connectie
from database_opzetten import database_opzetten

def main():
    if len(sys.argv) < 3:
        print("Roep het programma aan met: python3 main.py <QUERY.fasta> <SUBJECT.fasta>")
        sys.exit(1)

    query, subject = sys.argv[1], sys.argv[2]

    conn = maak_database_connectie()
    opgeslagen = 0
    try:
        database_opzetten(conn)
        print(f"Starten van BLAST pairwise met query: {query} en subject: {subject}")
        
        gn_lookup = build_gn_lookup(subject)
        data = enrich_hits(blast_pairwise(query, subject), gn_lookup)
        database_vullen(conn, data)
        opgeslagen = len(data)
    finally:
        conn.close()

    print(f"\nKlaar: {opgeslagen} records verwerkt en opgeslagen, resultaten: http://145.97.18.240:8081/")

if __name__ == '__main__':
    main()