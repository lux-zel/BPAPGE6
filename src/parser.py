from __future__ import annotations

from pathlib import Path
from shutil import copyfileobj
from urllib.request import urlopen
from Bio import SeqIO, SwissProt

DATA_DIR = Path("data")

def _download_file(url, output_file):
    """
    Download de file van de gegeven url.

    Deze functie downloadt een bestand van de meegegeven URL,
    en slaat het op in de opgegeven output_file. 
    Als het bestand al bestaat, wordt het niet opnieuw gedownload.

    Parameters:
    url (str): De URL van het bestand dat gedownload moet worden.
    output_file (Path): Het pad waar het bestand opgeslagen moet worden.

    Returns:
    output_file (Path): Het pad naar het gedownloade bestand, of None als er een fout optrad.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.is_file():
        return output_file

    try:
        with urlopen(url, timeout=60) as response:
            with output_file.open("wb") as file:
                copyfileobj(response, file)
                
    except Exception as e:
        if output_file.exists():
            output_file.unlink()
        print(f"Error downloading {url}: {e}")
        return None

    return output_file


def download_uniprot(accession):
    """Download de Uniprot annotaties voor een gegeven accession."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.txt"
    return _download_file(url, DATA_DIR / f"{accession}.txt")


def download_kegg(kegg_id):
    """Download de KEGG annotaties voor een gegeven KEGG ID."""
    url = f"https://rest.kegg.jp/get/{kegg_id}"
    filename = kegg_id.replace(":", "_")
    return _download_file(url, DATA_DIR / "kegg" / f"{filename}.txt")


def download_nucleotide(accession):
    """Download de nucleotide sequenties voor een gegeven accession."""
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}"
    return _download_file(url, DATA_DIR / "nucleotide" / f"{accession}.fasta")


def parse_uniprot(file_path):
    """
    Parse de Uniprot annotaties van een gegeven bestandspad.

    Deze functie leest een Uniprot annotatiebestand en extraheert relevante informatie zoals
    entry name, accessions, gene name, description, organism, sequence, gene IDs,
    locus tags, KEGG IDs, EMBL IDs, EMBL protein IDs en GO data.
    Voor de kegg_ids, embl_ids en go_data wordt er in de record.cross_references,
    gekeken naar de juiste referenties.

    Parameters:
    file_path (Path): Het pad naar het Uniprot annotatiebestand.

    Returns:
    annotations (dict): Een dictionary met de geparste annotaties.
    """
    annotations = {
        "uniprot_entry_name": None,
        "uniprot_accessions": [],
        "gene_name": None,
        "uniprot_description": None,
        "uniprot_organism": None,
        "uniprot_sequence": None,
        "gene_id": [],
        "gene_locus_tag": None,
        "kegg_ids": [],
        "embl_ids": [],
        "embl_protein_ids": [],
        "go_data": [],
    }

    with open(file_path, encoding="utf-8") as handle:
        record = SwissProt.read(handle)
        annotations.update({
            "uniprot_entry_name": record.entry_name,
            "uniprot_accessions": list(record.accessions),
            "gene_name": record.gene_name,
            "uniprot_description": record.description,
            "uniprot_organism": record.organism,
            "uniprot_sequence": record.sequence,
        })

        for gene_annotation in record.gene_name or []:
            for field in ("ORFNames", "OrderedLocusNames"):
                values = gene_annotation.get(field) or []
                if values and annotations["gene_locus_tag"] is None:
                    annotations["gene_locus_tag"] = values[0].split(" {", 1)[0]

        for ref in record.cross_references:
            if ref[0] == "GeneID" and len(ref) > 1:
                annotations["gene_id"].append(ref[1])
            elif ref[0] == "KEGG" and len(ref) > 1:
                annotations["kegg_ids"].append(ref[1])
            elif ref[0] == "EMBL" and len(ref) > 1:
                annotations["embl_ids"].append(ref[1])
                if len(ref) > 2:
                    annotations["embl_protein_ids"].append(ref[2])
            elif ref[0] == "GO" and len(ref) > 2:
                go_id = ref[1]
                if ":" in ref[2]:
                    go_type, go_term = ref[2].split(":", 1)
                else:
                    go_type, go_term = None, ref[2]
                annotations["go_data"].append((go_id, go_term, go_type))

    return annotations


def parse_kegg(file_path):
    """
    Parse de KEGG annotaties van een gegeven bestandspad.

    Parameters:
    file_path (Path): Het pad naar het KEGG annotatiebestand.

    Returns:
    fields (dict): Een dictionary met de geparste annotaties.
    """
    fields = {}
    current_field = None

    with open(file_path, encoding="utf-8") as handle:
        for line in handle:
            field = line[:12].strip()
            value = line[12:].strip()
            if field:
                fields.setdefault(field, []).append(value)
                current_field = field
            elif current_field and value:
                if current_field == "PATHWAY":
                    fields[current_field].append(value)
                else:
                    fields[current_field][-1] += f" {value}"

    return fields


def parse_nucleotide(file_path):
    """Parse de nucleotide sequenties van een gegeven bestandspad."""
    with open(file_path, encoding="utf-8") as handle:
        return [
            {
                "id": record.id,
                "description": record.description,
                "sequence": str(record.seq),
            }
            for record in SeqIO.parse(handle, "fasta")
        ]


def enrich_hits(hits):
    """
    De enrich_hits functie verrijkt de BLAST hits met annotaties van,
    Uniprot, KEGG en nucleotide databases. Het roept de download en parse functies aan 
    voor elk van deze databases, en voegt de annotaties/datasets toe aan de hits.

    Parameters:
    hits (list): Een lijst van BLAST hits, elk hit is een dictionary met een "accession" key.

    Returns:
    hits (list): De lijst van BLAST hits, verrijkt met annotaties.
    """
    uniprot_annotations = {}
    kegg_entries = {}
    nucleotide_records = {}

    for hit in hits:
        accession = hit["accession"]
        if accession not in uniprot_annotations:
            file_path = download_uniprot(accession)
            uniprot_annotations[accession] = (
                parse_uniprot(file_path) if file_path else {}
            )
        hit.update(uniprot_annotations[accession])

        kegg_data = {}
        for kegg_id in hit.get("kegg_ids", []):
            if kegg_id not in kegg_entries:
                file_path = download_kegg(kegg_id)
                kegg_entries[kegg_id] = parse_kegg(file_path) if file_path else None
            if kegg_entries[kegg_id] is not None:
                kegg_data[kegg_id] = kegg_entries[kegg_id]
        if kegg_data:
            hit["kegg_data"] = kegg_data

        nucleotide_data = {}
        for nucleotide_id in hit.get("embl_ids", []):
            if nucleotide_id not in nucleotide_records:
                sequence_file_path = download_nucleotide(nucleotide_id)
                nucleotide_records[nucleotide_id] = (
                    parse_nucleotide(sequence_file_path)
                    if sequence_file_path
                    else None
                )
            if nucleotide_records[nucleotide_id] is not None:
                nucleotide_data[nucleotide_id] = nucleotide_records[nucleotide_id]

        if nucleotide_data:
            hit["nucleotide_data"] = nucleotide_data

    return hits