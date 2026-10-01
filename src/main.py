import json
import sys
from pathlib import Path
from blast import blast_pairwise
from parser import enrich_hits

def main():

    if len(sys.argv) < 3:
        print("Roep het programma aan met: python3 main.py <QUERY.fasta> <SUBJECT.fasta>")
        sys.exit(1)
    
    query = sys.argv[1]
    subject = sys.argv[2]

    print(f"Starten van BLAST pairwise met query: {query} en subject: {subject}")
    data = enrich_hits(blast_pairwise(query, subject))

    output_file = Path(__file__).parent / "output" / "blast_results.json"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)

    print(f"Finished: {len(data)} HSP's verwerkt en opgeslagen")


if __name__ == '__main__':
    main()
