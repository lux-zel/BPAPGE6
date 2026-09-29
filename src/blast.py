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
    hits (dictionary):  Een dictionary met de gevonden hits van de BLAST.
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
                for hsp in alignment.hsps:
                    
                    data.append({
                        "accession": alignment.hit_def.split('|')[1],
                        "title": alignment.title,
                        "length": alignment.length,
                        "evalue": hsp.expect,
                        "bit_score": hsp.bits,
                        "score": hsp.score,
                        "identities": hsp.identities,
                        "positives": hsp.positives,
                        "align_length": hsp.align_length,
                        "gaps": hsp.gaps,
                        "query_start": hsp.query_start,
                        "query_end": hsp.query_end,
                        "sbjct_start": hsp.sbjct_start,
                        "sbjct_end": hsp.sbjct_end,
                        "query_seq": hsp.query,
                        "match_midline": hsp.match,
                        "sbjct_seq": hsp.sbjct,
                    })
    return data
