from __future__ import annotations

import json
import time
from pathlib import Path

import requests


KEGG_BASE = "https://rest.kegg.jp"

DATA_DIR = Path("data")

FIELDS = ",".join([
    "accession", "gene_primary", "organism_name", "organism_id",
    "protein_name", "length", "sequence", "cc_pathway",
    "go", "cc_function", "cc_alternative_products",
])


def kegg_conv_uniprot(accession):
    """Koppel een UniProt-accession aan KEGG-gen-ID's."""
    url = f"{KEGG_BASE}/conv/genes/uniprot:{accession}"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    time.sleep(0.2)
    kegg_ids = []
    for line in r.text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] == f"up:{accession}":
            kegg_ids.append(parts[1])
    return kegg_ids


def kegg_link_pathways(kegg_gene_id):
    """Geef de pathway-ID's die aan een KEGG-gen-ID gekoppeld zijn."""
    url = f"{KEGG_BASE}/link/pathway/{kegg_gene_id}"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    time.sleep(0.2)
    pathway_ids = []
    for line in r.text.strip().splitlines():
        parts = line.split("\t")
        if parts:
            pid = parts[1].replace("path:", "")
            pathway_ids.append(pid)
    return pathway_ids


def kegg_get_pathway_name(pathway_id):
    """Geef de naam van een KEGG-pathway."""
    url = f"{KEGG_BASE}/get/{pathway_id}"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    time.sleep(0.2)
    for line in r.text.splitlines():
        if line.startswith("NAME"):
            return line[12:].strip()
    return None


def _get(url, out_path):
    """Download JSON indien nodig en geef de inhoud van de cache terug."""
    if not out_path.is_file():
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(r.text, encoding="utf-8")
        time.sleep(0.2)
    return json.loads(out_path.read_text(encoding="utf-8"))


def download_uniprot(accession):
    """Download en geef het UniProt-JSON-record van een accession terug."""
    url = (
        f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        f"?fields={FIELDS}"
    )
    return _get(url, DATA_DIR / f"{accession}.json")


def _first_nested_value(items, *keys, default=None):
    """Geef de eerste waarde terug die via de opgegeven sleutels wordt gevonden."""
    for item in items or []:
        value = item
        for key in keys:
            value = value.get(key) if isinstance(value, dict) else value
            if value is None:
                break
        else:
            return value
    return default


def _comment_texts(comments, comment_type):
    """Geef de tekstwaarden van een UniProt-commentaartype terug."""
    return [
        text["value"]
        for comment in comments
        if comment.get("commentType") == comment_type
        for text in comment.get("texts", [])
    ]


def parse_uniprot(entry):
    """Haal de databaseannotaties uit een UniProt-record."""
    genes = entry.get("genes", [])
    comments = entry.get("comments", [])

    go_terms = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        for p in xref.get("properties", []):
            if p.get("key") == "GoTerm" and ":" in p["value"]:
                aspect, term = p["value"].split(":", 1)
                go_terms.append({
                    "go_id": xref["id"],
                    "aspect": aspect,
                    "term": term,
                })

    return {
        "gene_primary": _first_nested_value(genes, "geneName", "value"),
        "organism_name": entry.get("organism", {}).get("scientificName"),
        "organism_id": entry.get("organism", {}).get("taxonId"),
        "protein_name": _first_nested_value(
            [entry.get("proteinDescription", {})],
            "recommendedName",
            "fullName",
            "value",
        ),
        "length": entry.get("sequence", {}).get("length"),
        "sequence": entry.get("sequence", {}).get("value"),
        "cc_pathway": _comment_texts(comments, "PATHWAY"),
        "go_terms": go_terms,
        "cc_function": (_comment_texts(comments, "FUNCTION") or [None])[0],
    }


def build_gn_lookup(fasta_path):
    """Koppel UniProt-accessions aan genennamen uit FASTA-headers."""
    lookup = {}
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            parts = line[1:].split("|")
            if len(parts) < 3:
                continue
            accession = parts[1]
            for token in line.split():
                if token.startswith("GN="):
                    lookup[accession] = token[3:]
                    break
    return lookup


def build_isoform_records(fasta_path):
    """Lees isoformannotaties, accessions en genen uit FASTA-headers."""
    isoforms = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            parts = line[1:].strip().split("|", 2)
            if len(parts) < 3:
                continue
            description = parts[2].split(" ", 1)
            if len(description) < 2:
                continue
            protein_name = description[1].split(" OS=", 1)[0]
            if "isoform" not in protein_name.casefold():
                continue

            gene_id = next(
                (
                    token[3:]
                    for token in line.split()
                    if token.startswith("GN=")
                ),
                None,
            )
            organism = (
                line.split(" OS=", 1)[1].split(" OX=", 1)[0]
                if " OS=" in line and " OX=" in line
                else None
            )
            if gene_id:
                isoforms.append({
                    "accession": parts[1],
                    "gene_id": gene_id,
                    "organism_name": organism,
                    "isoform_name": protein_name,
                })
    return isoforms


def enrich_hits(hits, gn_lookup=None):
    """Voeg UniProt- en KEGG-annotaties toe aan BLAST-hits."""
    gn_lookup = gn_lookup or {}
    uniprot_cache = {}
    kegg_cache = {}

    for hit in hits:
        acc = hit["accession"]
        if acc not in uniprot_cache:
            uniprot_cache[acc] = parse_uniprot(download_uniprot(acc))
        hit.update(uniprot_cache[acc])

        if not hit.get("gene_primary"):
            hit["gene_primary"] = gn_lookup.get(acc)

        if acc not in kegg_cache:
            kegg_gene_ids = kegg_conv_uniprot(acc)
            pathways = {}
            for gene_id in kegg_gene_ids:
                for pid in kegg_link_pathways(gene_id):
                    if pid not in pathways:
                        pathways[pid] = kegg_get_pathway_name(pid)
            kegg_cache[acc] = pathways

        if kegg_cache[acc]:
            hit["kegg_pathways"] = kegg_cache[acc]

    return hits