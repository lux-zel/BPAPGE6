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
    
    Returns:
    data (list): Een lijst van dictionaries, een dict per BLAST hit.
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

                protein_id = alignment.Hit_def.split('|')[1]
                protein_length = alignment.Hit_len
                gene_id = alignment.hit_def.split('GN=')[1].split()[0]

                hsps_list = []
                for hsp in alignment.hsps:

                    hsp_start = hsp.query_start
                    hsp_end = hsp.query_end
                    query_seq = hsp.Hsp_qseq
                    sbjct_seq = hsp.Hsp_hseq
                    hsp_align = hsp.Hsp_midline

                    hsps_list.append({
                        "hsp_start": hsp_start,
                        "hsp_end": hsp_end,
                        "query_seq": query_seq,
                        "sbjct_seq": sbjct_seq,
                        "hsp_align": hsp_align
                    })
                    
            data.append({
                "protein_id": protein_id,
                "protein_length": protein_length,
                "gene_id": gene_id,
                "hsps": hsps_list
            })


    print(f"BLAST pairwise completed: {len(data)} hits gevonden.")
    print(f"Start nu met het verrijken van de hits met annotaties van Uniprot, KEGG en NCBI Nucleotide. \
          \nDit kan even duren...")
    return data
