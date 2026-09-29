from pathlib import Path
from urllib.request import urlopen
from Bio import SwissProt
import certifi
import ssl

DATA_DIR = Path("data")
OUPUT_DIR = Path("output")

DATA_DIR.mkdir(exist_ok=True)
OUPUT_DIR.mkdir(exist_ok=True)

def download_uniprot(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.txt"
    output_file = DATA_DIR / f"{accession}.txt"

    ssl_context = ssl.create_default_context(cafile=certifi.where())

    try:
        with urlopen(url, context=ssl_context) as response:
            with open(output_file, "wb") as file:
                file.write(response.read())
        print(f"{accession} is gedownload")
    except Exception as e:
        print(f"Error downloading {accession}: {e}")
        return None

    return output_file


def parse_uniprot(file_path, hit):
    with open(file_path, encoding="utf-8") as handle:
        record = SwissProt.read(handle)

        for ref in record.cross_references:
            if ref[0] == "GeneID":
                hit.update({"gene_id": ref[1]})
            if ref[0] == "KEGG":
                hit.update({"kegg_ids": ref[1]})
            if ref[0] == "Ensembl":
                hit.update({"ensembl_ids": ref[1]})
            elif ref[0] == "EMBL":
                hit.update({"embl_ids": ref[1]})

            if ref[0] == "GO":
                go_id = ref[1]
                go_type, go_term = ref[2].split(':', 1)
                hit.update({"go_data": (go_id, go_term, go_type)})

        if len(hit.get("go_data", [])) == 0:
            go_id = None
            go_type = None
            go_term = None
            hit.update({"go_data": (go_id, go_term, go_type)})
        if len(hit.get("ensembl_ids", [])) == 0:
            hit.update({"ensembl_ids": None})
        if len(hit.get("embl_ids", [])) == 0:
            hit.update({"embl_ids": None})
        if len(hit.get("kegg_ids", [])) == 0:
            hit.update({"kegg_ids": None})
        if len(hit.get("gene_id", [])) == 0:
            hit.update({"gene_id": None})