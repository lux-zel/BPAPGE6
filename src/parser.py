from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from shutil import copyfileobj
from urllib.parse import quote
from urllib.request import urlopen
from Bio import SeqIO, SwissProt
import certifi
import ssl

DATA_DIR = Path("data")

def _download_file(url: str, output_file: Path) -> Path | None:
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


def download_uniprot(accession: str) -> Path | None:
    url = f"https://rest.uniprot.org/uniprotkb/{quote(accession, safe='')}.txt"
    return _download_file(url, DATA_DIR / f"{accession}.txt")


def download_kegg(kegg_id: str) -> Path | None:
    url = f"https://rest.kegg.jp/get/{quote(kegg_id, safe=':')}"
    filename = kegg_id.replace(":", "_")
    return _download_file(url, DATA_DIR / "kegg" / f"{filename}.txt")


def download_nucleotide(accession: str) -> Path | None:
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{quote(accession, safe='')}"
    return _download_file(url, DATA_DIR / "nucleotide" / f"{accession}.fasta")


def download_nucleotide_annotations(accession: str) -> Path | None:
    url = f"https://www.ebi.ac.uk/ena/browser/api/embl/{quote(accession, safe='')}"
    return _download_file(url, DATA_DIR / "nucleotide" / f"{accession}.embl")


def parse_uniprot(file_path: Path) -> dict:
    annotations = {
        "uniprot_entry_name": None,
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


def parse_kegg(file_path: Path) -> dict:
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


def parse_nucleotide(
    annotation_file_path: Path, sequence_file_path: Path | None
) -> dict:
    genes_by_protein_id = {}
    isomers_by_locus_tag = {}
    sequence_record = None
    if sequence_file_path:
        with open(sequence_file_path, encoding="utf-8") as handle:
            sequence_record = next(SeqIO.parse(handle, "fasta"), None)

    with open(annotation_file_path, encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "embl"):
            gene_features = [
                feature for feature in record.features if feature.type == "gene"
            ]
            cds_features = [
                feature for feature in record.features if feature.type == "CDS"
            ]

            for transcript in record.features:
                if transcript.type != "mRNA":
                    continue

                locus_tag = _first_qualifier(transcript, "locus_tag")
                if not locus_tag:
                    continue

                transcript_id = _first_qualifier(transcript, "transcript_id")
                location = str(transcript.location)
                variant_id = transcript_id or sha256(
                    location.encode("utf-8")
                ).hexdigest()[:12]
                isomer_id = f"{record.id}:{locus_tag}:{variant_id}"
                isomers_by_locus_tag.setdefault(locus_tag, []).append({
                    "isomer_id": isomer_id,
                    "sequence": (
                        str(transcript.extract(sequence_record.seq))
                        if sequence_record
                        else None
                    ),
                    "splicing_variant": location,
                    "codes": any(
                        _first_qualifier(cds, "locus_tag") == locus_tag
                        for cds in cds_features
                    ),
                })

            for cds in cds_features:
                protein_id = _first_qualifier(cds, "protein_id")
                if not protein_id:
                    continue

                locus_tag = _first_qualifier(cds, "locus_tag")
                gene_feature = next(
                    (
                        feature
                        for feature in gene_features
                        if locus_tag
                        and _first_qualifier(feature, "locus_tag") == locus_tag
                    ),
                    None,
                )
                gene_feature = gene_feature or cds
                start = int(gene_feature.location.start) + 1
                end = int(gene_feature.location.end)
                genes_by_protein_id[protein_id] = {
                    "locus_tag": locus_tag,
                    "sequence": (
                        str(gene_feature.extract(sequence_record.seq))
                        if sequence_record
                        else None
                    ),
                    "start": start,
                    "end": end,
                    "strand": gene_feature.location.strand,
                    "description": _first_qualifier(cds, "product"),
                }
    return {
        "genes_by_protein_id": genes_by_protein_id,
        "isomers_by_locus_tag": isomers_by_locus_tag,
    }


def _first_qualifier(feature, name: str) -> str | None:
    values = feature.qualifiers.get(name) or []
    return values[0] if values else None


def enrich_hits(hits: list) -> list:
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
                annotation_file_path = download_nucleotide_annotations(nucleotide_id)
                sequence_file_path = download_nucleotide(nucleotide_id)
                nucleotide_records[nucleotide_id] = (
                    parse_nucleotide(annotation_file_path, sequence_file_path)
                    if annotation_file_path
                    else None
                )
            if nucleotide_records[nucleotide_id] is not None:
                nucleotide_data[nucleotide_id] = nucleotide_records[nucleotide_id]

        for record in nucleotide_data.values():
            gene_data = next(
                (
                    record["genes_by_protein_id"][protein_id]
                    for protein_id in hit.get("embl_protein_ids", [])
                    if protein_id in record["genes_by_protein_id"]
                ),
                None,
            )
            if gene_data:
                hit["gene_data"] = gene_data
                break

        locus_tag = (hit.get("gene_data") or {}).get("locus_tag")
        if locus_tag:
            isomers = {}
            for record in nucleotide_data.values():
                for isomer in record["isomers_by_locus_tag"].get(locus_tag, []):
                    isomers[isomer["isomer_id"]] = isomer
            hit["isomers"] = list(isomers.values())

    return hits