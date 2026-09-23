"""
dit bestand is een helper python bestand dat een blast opdracht uitvoert
van een query tegen een subject. 
"""


import subprocess
from pathlib import Path

from Bio import Blast

def blast_pairwise(query, subject):
    """
    Deze function maakt gebruikt van Biopython's wrapper rond de BLAST+ binaries
    Het voert een local paiwise BLAST uit met -subject in plaats van een database.
    
    Parameters:
    query (string):     De fasta sequentie van de sequentie die je Queried.
    subject (string):   De fasta sequentie waartegen je Queried.
    seq_type (string):  Het type allignment ("nucl" voor blastn, "prot" voor blastp).
    
    Returns:
    hits (dictionary):  Een dictionary met de gevonden hits van de BLAST.
    
    """
    

    project_dir = Path(__file__).parent
    subprocess.run(['bash', str(project_dir / 'blast.sh')], cwd=project_dir, check=True)

    
    
    hits = []
    with open(project_dir / 'hits.xml') as xml_file:
        records = Blast.parse(xml_file)

    for record in records:
        for hit in record:
            for hsp in hit:
                
                    hits.append({
                            
                            "evalue": hsp.expect,
                            "bit-score": hsp.bits,
                            "identities": hsp.identities,
                            "align_length": hsp.align_length,
                            "gaps": hsp.gaps,
                            "query_start": hsp.query_start,
                            "query_end": hsp.query_end,
                            "subject_start": hsp.sbjct_start,
                            "subject_end": hsp.sbjct_end,
                            "query_seq": hsp.query,
                            "match": hsp.match,
                            "subject_seq": hsp.sbjct
                            })
    return hits
