
import os
from collections.abc import Iterable, Mapping
from typing import Any

import psycopg2


def fill_database(data: Iterable[Mapping[str, Any]]) -> int:
    """Voert de BLAST HSPs en de ge-associeerde protein annotations in de database.

    De database connection is geconfigureerd door de ``DATABASE_URL``
    environment variable. Alle records zijn geschreven in één transactie.
    """
    conn_string = os.environ.get("DATABASE_URL")
    if not conn_string:
        raise RuntimeError("Zet deDATABASE_URL environment variable eerst.")

    conn = psycopg2.connect(conn_string)
    inserted = 0
    try:
        with conn.cursor() as cursor:
            for hit in data:
                protein_id = _store_protein(hit, cursor)
                _store_go_data(hit, protein_id, cursor)
                _store_pathways(hit, protein_id, cursor)
                _store_hsp(hit, protein_id, cursor)
                inserted += 1
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()

    return inserted


def _store_protein(hit: Mapping[str, Any], cursor: Any) -> int:
    gene_ids = hit.get("gene_id") or []
    gene_id = next((value for value in gene_ids if value), None)

    if gene_id is not None:
        cursor.execute(
            "INSERT INTO gen (gen_id) VALUES (%s) ON CONFLICT (gen_id) DO NOTHING",
            (gene_id,),
        )

    cursor.execute(
        """
        INSERT INTO eiwit (eiwit_id, gen_id, naam, aminozuursequentie)
        VALUES (%s, %s, %s, %s)
        ON CONFLICT (eiwit_id) DO UPDATE SET
            gen_id = COALESCE(EXCLUDED.gen_id, eiwit.gen_id),
            naam = COALESCE(EXCLUDED.naam, eiwit.naam),
            aminozuursequentie = COALESCE(
                EXCLUDED.aminozuursequentie, eiwit.aminozuursequentie
            )
        RETURNING id
        """,
        (
            hit["accession"],
            gene_id,
            hit.get("uniprot_entry_name"),
            hit.get("uniprot_sequence"),
        ),
    )
    row = cursor.fetchone()
    if row is None:
        raise RuntimeError(f"Could not retrieve database ID for {hit['accession']}.")
    return row[0]


def _store_go_data(hit: Mapping[str, Any], protein_id: int, cursor: Any) -> None:
    for go_id, go_term, go_type in hit.get("go_data") or []:
        cursor.execute(
            """
            INSERT INTO functie (go_id, go_term, go_type)
            VALUES (%s, %s, %s)
            ON CONFLICT (go_id) DO UPDATE SET
                go_term = COALESCE(EXCLUDED.go_term, functie.go_term),
                go_type = COALESCE(EXCLUDED.go_type, functie.go_type)
            """,
            (go_id, go_term, go_type),
        )
        cursor.execute(
            """
            INSERT INTO eiwit_functie (id, go_id)
            VALUES (%s, %s)
            ON CONFLICT (id, go_id) DO NOTHING
            """,
            (protein_id, go_id),
        )


def _first_value(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


def _store_pathways(hit: Mapping[str, Any], protein_id: int, cursor: Any) -> None:
    kegg_data = hit.get("kegg_data") or {}
    for kegg_id in hit.get("kegg_ids") or []:
        pathway = kegg_data.get(kegg_id) or {}
        cursor.execute(
            """
            INSERT INTO pathway (
                pathway_kegg_id, naam, type_pathway, beschrijving
            )
            VALUES (%s, %s, %s, %s)
            ON CONFLICT (pathway_kegg_id) DO UPDATE SET
                naam = COALESCE(EXCLUDED.naam, pathway.naam),
                type_pathway = COALESCE(EXCLUDED.type_pathway, pathway.type_pathway),
                beschrijving = COALESCE(EXCLUDED.beschrijving, pathway.beschrijving)
            """,
            (
                kegg_id,
                _first_value(pathway.get("NAME")),
                _first_value(pathway.get("CLASS")),
                _first_value(pathway.get("DESCRIPTION")),
            ),
        )
        cursor.execute(
            """
            INSERT INTO eiwit_pathway (eiwit_id, pathway_id)
            VALUES (%s, %s)
            ON CONFLICT (eiwit_id, pathway_id) DO NOTHING
            """,
            (protein_id, kegg_id),
        )


def _store_hsp(hit: Mapping[str, Any], protein_id: int, cursor: Any) -> None:
    cursor.execute(
        """
        INSERT INTO blast_hsp (
            eiwit_id, blast_id, title, length, query_id, query_title,
            query_length, evalue, bit_score, score, identities, positives,
            align_length, gaps, query_start, query_end, query_frame,
            sbjct_start, sbjct_end, sbjct_frame, query_seq, match_midline,
            sbjct_seq
        )
        VALUES (
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
            %s, %s, %s, %s, %s, %s, %s, %s
        )
        """,
        (
            protein_id,
            hit["blast_id"],
            hit["title"],
            hit["length"],
            hit["query_id"],
            hit["query_title"],
            hit["query_length"],
            hit["evalue"],
            hit["bit_score"],
            hit["score"],
            hit["identities"],
            hit["positives"],
            hit["align_length"],
            hit["gaps"],
            hit["query_start"],
            hit["query_end"],
            hit["query_frame"],
            hit["sbjct_start"],
            hit["sbjct_end"],
            hit["sbjct_frame"],
            hit["query_seq"],
            hit["match_midline"],
            hit["sbjct_seq"],
        ),
    )
