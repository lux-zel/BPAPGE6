"""
Dit bestand bevat blast_pairwise(), een helper functio dat een blast opdracht uitvoert
van een query tegen een subject, waar de query en subject fasta bestanden zijn en er 
een database wordt gemaakt van het subject bestand.
Het gebruikt Biopython's wrapper rond de BLAST+ binaries.
"""

import os
import subprocess
from io import StringIO
from pathlib import Path
from Bio import Blast
from Bio import SeqIO

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
        ['bash', str(project_dir / 'blast.sh')],
        cwd=project_dir,
        env=env,
        check=True,
    )

    hits = []
    with open(project_dir / 'hits.xml', 'rb') as xml_file:
        records = Blast.parse(xml_file)
        for record in records:
            for hit in record:
                for hsp in hit:
                    counts = hsp.counts()
                    query_start, query_end = hsp.coordinates[1, [0, -1]]
                    subject_start, subject_end = hsp.coordinates[0, [0, -1]]
                    aligned_sequences = list(
                        SeqIO.parse(StringIO(hsp.format("fasta")), "fasta")
                    )
                    hits.append({
                        "evalue": hsp.annotations["evalue"],
                        "bit-score": hsp.annotations["bit score"],
                        "identities": hsp.annotations["identity"],
                        "align_length": counts.aligned,
                        "gaps": counts.gaps,
                        "query_start": int(query_start),
                        "query_end": int(query_end),
                        "subject_start": int(subject_start),
                        "subject_end": int(subject_end),
                        "query_seq": str(aligned_sequences[1].seq),
                        "match": hsp.annotations["midline"],
                        "subject_seq": str(aligned_sequences[0].seq),
                    })
    return hits
