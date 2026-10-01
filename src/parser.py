from pathlib import Path
from shutil import copyfileobj
from urllib.parse import quote
from urllib.request import urlopen
from Bio import SeqIO, SwissProt
import certifi
import ssl

DATA_DIR = Path("data")

def _download_file(url, output_file):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.is_file():
        return output_file

    temporary_file = output_file.with_name(f".{output_file.name}.part")
    try:
        ssl_context = ssl.create_default_context(cafile=certifi.where())
        with urlopen(url, context=ssl_context, timeout=60) as response:
            with temporary_file.open("wb") as file:
                copyfileobj(response, file)
        temporary_file.replace(output_file)
    except Exception as e:
        if temporary_file.exists():
            temporary_file.unlink()
        print(f"Error downloading {url}: {e}")
        return None

    return output_file


def download_uniprot(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{quote(accession, safe='')}.txt"
    return _download_file(url, DATA_DIR / f"{accession}.txt")


def download_kegg(kegg_id):
    url = f"https://rest.kegg.jp/get/{quote(kegg_id, safe=':')}"
    filename = kegg_id.replace(":", "_")
    return _download_file(url, DATA_DIR / "kegg" / f"{filename}.txt")


def download_nucleotide(accession):
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{quote(accession, safe='')}"
    return _download_file(url, DATA_DIR / "nucleotide" / f"{accession}.fasta")


def parse_uniprot(file_path):
    annotations = {
        "uniprot_entry_name": None,
        "uniprot_accessions": [],
        "uniprot_description": None,
        "uniprot_organism": None,
        "uniprot_sequence": None,
        "uniprot_cross_references": [],
        "gene_id": [],
        "kegg_ids": [],
        "ensembl_ids": [],
        "embl_ids": [],
        "embl_protein_ids": [],
        "embl_references": [],
        "go_data": [],
    }

    with open(file_path, encoding="utf-8") as handle:
        record = SwissProt.read(handle)
        annotations.update({
            "uniprot_entry_name": record.entry_name,
            "uniprot_accessions": list(record.accessions),
            "uniprot_description": record.description,
            "uniprot_organism": record.organism,
            "uniprot_sequence": record.sequence,
            "uniprot_cross_references": list(record.cross_references),
        })

        for ref in record.cross_references:
            if ref[0] == "GeneID" and len(ref) > 1:
                annotations["gene_id"].append(ref[1])
            elif ref[0] == "KEGG" and len(ref) > 1:
                annotations["kegg_ids"].append(ref[1])
            elif ref[0] == "Ensembl" and len(ref) > 1:
                annotations["ensembl_ids"].append(ref[1])
            elif ref[0] == "EMBL" and len(ref) > 1:
                annotations["embl_ids"].append(ref[1])
                if len(ref) > 2:
                    annotations["embl_protein_ids"].append(ref[2])
                annotations["embl_references"].append(ref)

            if ref[0] == "GO" and len(ref) > 2:
                go_id = ref[1]
                if ":" in ref[2]:
                    go_type, go_term = ref[2].split(":", 1)
                else:
                    go_type, go_term = None, ref[2]
                annotations["go_data"].append((go_id, go_term, go_type))

    return annotations


def parse_kegg(file_path):
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
                fields[current_field][-1] += f" {value}"

    return fields


def parse_nucleotide(file_path):
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
                file_path = download_nucleotide(nucleotide_id)
                nucleotide_records[nucleotide_id] = (
                    parse_nucleotide(file_path) if file_path else None
                )
            if nucleotide_records[nucleotide_id] is not None:
                nucleotide_data[nucleotide_id] = nucleotide_records[nucleotide_id]
        if nucleotide_data:
            hit["nucleotide_data"] = nucleotide_data

    return hits