from pathlib import Path
from urllib.request import urlopen
from Bio import SwissProt
import certifi
import ssl
import csv

DATA_DIR = Path("data")
OUPUT_DIR = Path("output")

DATA_DIR.mkdir(exist_ok=True)
OUPUT_DIR.mkdir(exist_ok=True)


def download_uniprot(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.txt"
    output_file = DATA_DIR / f"{accession}.txt"

    ssl_context = ssl.create_default_context(cafile=certifi.where())

    with urlopen(url, context=ssl_context) as response:
        with open(output_file, "wb") as file:  # wb = "write binary"
            file.write(response.read())

    print(f"{accession} is gedownload")

    return output_file


def parse_uniprot(file_path):
    with open(file_path, encoding="utf=8") as handle:
        record = SwissProt.read(handle)

        go_data = []
        ensembl_ids = []
        embl_ids = []
        kegg_ids = []
        gene_id = []

        for ref in record.cross_references:
            if ref[0] == "GeneID":
                gene_id.append(ref[1])
            if ref[0] == "KEGG":
                kegg_ids.append(ref[1])

            if ref[0] == "Ensembl":
                ensembl_ids.append(ref[1])

            elif ref[0] == "EMBL":
                embl_ids.append(ref[1])

            if ref[0] == "GO":
                go_id = ref[1]
                go_type, go_term = ref[2].split(':', 1)
                go_data.append((go_id, go_term, go_type))

        if len(go_data) == 0:
            go_id = None
            go_type = None
            go_term = None
            go_data.append((go_id, go_term, go_type))
        if len(ensembl_ids) == 0:
            ensembl_ids.append(None)
        if len(embl_ids) == 0:
            embl_ids.append(None)
        if len(kegg_ids) == 0:
            kegg_ids.append(None)
        if len(gene_id) == 0:
            gene_id.append(None)

        print("UniProt accession:")
        print(record.accessions[0])

        print("\nEiwitnaam:")
        print(record.description)

        print("\nGen:")
        print(record.gene_name)

        print("\nOrganisme:")
        print(record.organism)

        print("\nAminozuursequentie:")
        print(record.sequence)


file = download_uniprot("M7BTR0")

parse_uniprot(file)
