
import os

import psycopg2


def fill_database(data):
    """Stores gene, protein, GO, and KEGG annotations in one transaction.

    De database connection is geconfigureerd door de ``DATABASE_URL``
    environment variable. Alle records zijn geschreven in één transactie.
    """
    conn_string = os.environ.get("DATABASE_URL")
    if not conn_string:
        raise RuntimeError("Zet de DATABASE_URL environment variable eerst.")

    conn = psycopg2.connect(conn_string)
    inserted = 0
    try:
        with conn.cursor() as cursor:
            for hit in data:
                gene_id = _store_gene(hit, cursor)
                protein_id = _store_protein(hit, gene_id, cursor)
                _store_go_data(hit, protein_id, cursor)
                _store_pathways(hit, protein_id, cursor)
                inserted += 1
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()

    return inserted


def _store_gene(hit, cursor):
    gene_ids = hit.get("gene_id") or []
    gene_id = next((value for value in gene_ids if value), None)
    if gene_id is None:
        return None

    cursor.execute(
        """
        INSERT INTO gen (gen_id, naam, organisme, ensembl_ids)
        VALUES (%s, %s, %s, %s)
        ON CONFLICT (gen_id) DO UPDATE SET
            naam = COALESCE(EXCLUDED.naam, gen.naam),
            organisme = COALESCE(EXCLUDED.organisme, gen.organisme),
            ensembl_ids = ARRAY(
                SELECT DISTINCT unnest(gen.ensembl_ids || EXCLUDED.ensembl_ids)
            )
        """,
        (
            gene_id,
            hit.get("gene_name"),
            hit.get("uniprot_organism"),
            hit.get("ensembl_ids") or [],
        ),
    )
    return gene_id


def _store_protein(hit, gene_id, cursor):
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


def _store_go_data(hit, protein_id, cursor):
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


def _first_value(value):
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


def _store_pathways(hit, protein_id, cursor):
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
