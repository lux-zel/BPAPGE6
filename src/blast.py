"""
Dit bestand bevat blast_pairwise(), een helper functio dat een blast opdracht uitvoert
van een query tegen een subject, waar de query en subject fasta bestanden zijn en er 
een database wordt gemaakt van het subject bestand.
Het gebruikt Biopython's wrapper rond de BLAST+ binaries.
"""

import os
import subprocess
from pathlib import Path
from Bio.Blast import NCBIXML

def blast_pairwise(query, subject):
    """
    Deze function maakt gebruikt van Biopython's wrapper rond de BLAST+ binaries.
    Het voert een local paiwise BLAST uit met een locaal aangemaakte database.
    
    Parameters:
    query (string):     De fasta sequentie van de sequentie die je Queried.
    subject (string):   De fasta sequentie waartegen je Queried.
    seq_type (string):  Het type allignment ("nucl" voor blastn, "prot" voor blastp).
    
    Returns:
    Een lijst met dictionaries, een dictionary per HSP.
    """

    project_dir = Path(__file__).parent

    env = {
        **os.environ,
        "QUERY": str(Path(query).resolve()),
        "SUBJECT": str(Path(subject).resolve()),
    }

    subprocess.run(
        ['sh', str(project_dir / 'blast.sh')],
        env=env, 
        check=True
        )

    data = []
    with open('hits.xml', 'rb') as xml_file:
        records = NCBIXML.parse(xml_file)

        for record in records:
            for alignment in record.alignments:                
                hit_parts = alignment.hit_def.split("|")
                if len(hit_parts) >= 3:
                    hit_accession = hit_parts[1]
                else:
                    hit_accession = hit_parts[0] if hit_parts else alignment.accession

                for hsp in alignment.hsps:
                    data.append({
                        "accession": hit_accession,
                        "blast_id": alignment.hit_id,
                        "title": alignment.title,
                        "length": alignment.length,
                        "query_id": record.query_id,
                        "query_title": record.query,
                        "query_length": record.query_length,
                        "evalue": hsp.expect,
                        "bit_score": hsp.bits,
                        "score": hsp.score,
                        "identities": hsp.identities,
                        "positives": hsp.positives,
                        "align_length": hsp.align_length,
                        "gaps": hsp.gaps,
                        "query_start": hsp.query_start,
                        "query_end": hsp.query_end,
                        "query_frame": hsp.frame[0] if hsp.frame else None,
                        "sbjct_start": hsp.sbjct_start,
                        "sbjct_end": hsp.sbjct_end,
                        "sbjct_frame": hsp.frame[1] if hsp.frame else None,
                        "query_seq": hsp.query,
                        "match_midline": hsp.match,
                        "sbjct_seq": hsp.sbjct,
                    })
    print(f"BLAST pairwise completed: {len(data)} HSP's gevonden.")
    print(f"Start nu met het verrijken van de HSP's met annotaties van Uniprot, KEGG en NCBI Nucleotide. \
          Dit kan even duren...")
    return data
