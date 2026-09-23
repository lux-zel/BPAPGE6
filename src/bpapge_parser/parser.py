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

    urlretrieve(url, output_file)

    print(f"{accession} is gedownload")

    return output_file


def parse_uniprot(file_path):
    with open(file_path, encoding="utf=8") as handle:
        record = SwissProt.read(handle)

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
