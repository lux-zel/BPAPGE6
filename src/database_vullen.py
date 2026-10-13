
from psycopg2.extras import Json


def fill_database(data, conn):
    """Bewaar BLAST-hits en annotaties in de database."""
    inserted = 0
    with conn:
        with conn.cursor() as cursor:
            for hit in data:
                protein_id = _store_protein(hit, cursor)
                _store_go_data(hit, protein_id, cursor)
                _store_pathways(hit, protein_id, cursor)
                inserted += 1
    return inserted


def _store_protein(hit, cursor):
    gene_ids = hit.get("gene_id") or []
    gene_id = gene_ids[0] if gene_ids else None

    if gene_id is not None:
        cursor.execute(
            "INSERT INTO gen (gen_id) VALUES (%s) ON CONFLICT (gen_id) DO NOTHING",
            (gene_id,),
        )

    cursor.execute(
        """
        INSERT INTO eiwit (
            eiwit_id, gen_id, naam, aminozuursequentie, hit_data
        )
        VALUES (%s, %s, %s, %s, %s)
        ON CONFLICT (eiwit_id) DO UPDATE SET
            gen_id = COALESCE(EXCLUDED.gen_id, eiwit.gen_id),
            naam = COALESCE(EXCLUDED.naam, eiwit.naam),
            aminozuursequentie = COALESCE(
                EXCLUDED.aminozuursequentie, eiwit.aminozuursequentie
            ),
            hit_data = eiwit.hit_data || EXCLUDED.hit_data
        RETURNING id
        """,
        (
            hit["accession"],
            gene_id,
            hit.get("uniprot_entry_name"),
            hit.get("uniprot_sequence"),
            Json([dict(hit)]),
        ),
    )
    return cursor.fetchone()[0]


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
    return value[0] if value else None


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
