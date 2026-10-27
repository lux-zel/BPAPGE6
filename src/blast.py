"""Voer BLAST uit op een lokale eiwitdatabase en verwerk de XML-resultaten."""

import os
import subprocess
from pathlib import Path

from Bio.Blast import NCBIXML


def blast_pairwise(query, subject):
    """Voer BLAST uit voor de query en geef maximaal vijf hits terug.

    Elke hit bevat de beste HSP van één uitlijning.
    """
    project_dir = Path(__file__).parent
    xml_path = project_dir / "hits.xml"

    env = {
        **os.environ,
        "QUERY": str(Path(query).resolve()),
        "SUBJECT": str(Path(subject).resolve()),
        "OUT": str(xml_path),
    }
    subprocess.run(
        ["sh", str(project_dir / "blast.sh")],
        env=env,
        check=True,
    )

    data = []
    with open(xml_path, "rb") as xml_file:
        for record in NCBIXML.parse(xml_file):
            for alignment in record.alignments[:5]:
                parts = alignment.hit_id.split("|")
                if len(parts) < 3:
                    continue
                hsp = alignment.hsps[0] if alignment.hsps else None
                data.append({
                    "accession": parts[1],
                    "query_id": record.query_id,
                    "evalue": hsp.expect if hsp else None,
                    "bitscore": hsp.bits if hsp else None,
                    "score": hsp.score if hsp else None,
                    "identities": hsp.identities if hsp else None,
                    "align_length": hsp.align_length if hsp else None,
                    "gaps": hsp.gaps if hsp else None,
                    "query_start": hsp.query_start if hsp else None,
                    "query_end": hsp.query_end if hsp else None,
                    "query_frame": hsp.frame[0] if hsp else None,
                    "subject_start": hsp.sbjct_start if hsp else None,
                    "subject_end": hsp.sbjct_end if hsp else None,
                    "subject_frame": hsp.frame[1] if hsp else None,
                })

    print(f"\nBLAST pairwise completed: {len(data)} hits gevonden.")
    print(
        "Start nu met het verrijken van de hits "
        "met annotaties van UniProt en KEGG."
    )
    print("Dit kan even duren...")
    return data