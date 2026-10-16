import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from src.database_vullen import _store_pathways
from src.parser import parse_kegg


class ParseKeggPathwaysTests(unittest.TestCase):
    def test_keeps_continued_pathway_records_separate(self):
        kegg_record = """ENTRY       102944265         CDS       T02989
NAME        (RefSeq) ferritin heavy chain A
            continued name
PATHWAY     cmy00860  Porphyrin metabolism
            cmy04216  Ferroptosis
///
"""
        with tempfile.TemporaryDirectory() as directory:
            file_path = Path(directory) / "kegg.txt"
            file_path.write_text(kegg_record, encoding="utf-8")

            parsed = parse_kegg(file_path)

        self.assertEqual(
            parsed["PATHWAY"],
            ["cmy00860  Porphyrin metabolism", "cmy04216  Ferroptosis"],
        )
        self.assertEqual(
            parsed["NAME"],
            ["(RefSeq) ferritin heavy chain A continued name"],
        )

    def test_stores_pathway_ids_instead_of_gene_entry_ids(self):
        cursor = Mock()
        _store_pathways(
            {
                "kegg_ids": ["cmy:102944265"],
                "kegg_data": {
                    "cmy:102944265": {
                        "PATHWAY": [
                            "cmy00860  Porphyrin metabolism",
                            "cmy04216  Ferroptosis",
                        ]
                    }
                },
            },
            protein_id=7,
            cursor=cursor,
        )

        parameters = [call.args[1] for call in cursor.execute.call_args_list]
        self.assertEqual(
            parameters,
            [
                ("cmy00860", "Porphyrin metabolism", None, None),
                (7, "cmy00860"),
                ("cmy04216", "Ferroptosis", None, None),
                (7, "cmy04216"),
            ],
        )


if __name__ == "__main__":
    unittest.main()
