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
    """Map a UniProt accession to KEGG gene IDs via /conv/genes/."""
    url = f"{KEGG_BASE}/conv/genes/uniprot:{accession}"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    time.sleep(0.2)
    # Response lines: "up:P01308\thsa:3630"
    kegg_ids = []
    for line in r.text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] == f"up:{accession}":
            kegg_ids.append(parts[1])   # "hsa:3630"
    return kegg_ids


def kegg_link_pathways(kegg_gene_id):
    """Get pathway IDs for a KEGG gene ID via /link/pathway/."""
    url = f"{KEGG_BASE}/link/pathway/{kegg_gene_id}"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    time.sleep(0.2)
    # Response lines: "path:cmy00071\tcmy:102946325"
    pathway_ids = []
    for line in r.text.strip().splitlines():
        parts = line.split("\t")
        if parts:
            pid = parts[0].replace("path:", "")
            pathway_ids.append(pid)
    return pathway_ids


def kegg_get_pathway_name(pathway_id):
    """Get pathway name from /get/ (flat text, NAME line)."""
    url = f"{KEGG_BASE}/get/{pathway_id}"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    time.sleep(0.2)
    for line in r.text.splitlines():
        if line.startswith("NAME"):
            return line[12:].strip()
    return None


def _get(url, out_path):
    """Download url to out_path (skip if already there). Return parsed JSON."""
    if not out_path.is_file():
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(r.text, encoding="utf-8")
        time.sleep(0.2)
    return json.loads(out_path.read_text(encoding="utf-8"))


def download_uniprot(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json?fields={FIELDS}"
    return _get(url, DATA_DIR / f"{accession}.json")


def parse_uniprot(entry):
    """Pull the fields we care about out of a UniProt JSON entry."""
    def first(lst, *keys, default=None):
        for item in lst or []:
            for k in keys:
                item = item.get(k) if isinstance(item, dict) else item
                if item is None:
                    break
            else:
                return item
        return default

    genes = entry.get("genes", [])
    comments = entry.get("comments", [])

    def comment_texts(ctype):
        return [t["value"] for c in comments if c.get("commentType") == ctype
                for t in c.get("texts", [])]

    isoforms = [i.get("name", {}).get("value")
                for c in comments if c.get("commentType") == "ALTERNATIVE PRODUCTS"
                for i in c.get("isoforms", [])]

    go_terms = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        for p in xref.get("properties", []):
            if p.get("key") == "GoTerm" and ":" in p["value"]:
                aspect, term = p["value"].split(":", 1)
                go_terms.append({"go_id": xref["id"], "aspect": aspect, "term": term})

    return {
        "gene_primary":       first(genes, "geneName", "value"),
        "organism_name":      entry.get("organism", {}).get("scientificName"),
        "organism_id":        entry.get("organism", {}).get("taxonId"),
        "protein_name":       first([entry.get("proteinDescription", {})],
                                    "recommendedName", "fullName", "value"),
        "length":             entry.get("sequence", {}).get("length"),
        "sequence":           entry.get("sequence", {}).get("value"),
        "cc_pathway":         comment_texts("PATHWAY"),
        "go_terms":           go_terms,
        "cc_function":        (comment_texts("FUNCTION") or [None])[0],
        "isoforms":           isoforms,
    }


def build_gn_lookup(fasta_path):
    """Map UniProt accession -> GN= value from proteome FASTA headers."""
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


def enrich_hits(hits, gn_lookup=None):
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

        # Direct KEGG lookup by UniProt accession
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