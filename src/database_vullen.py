
import os

import psycopg2


def fill_database(data):
    """
    Vul de database met de verrijkte BLAST hits.
    elk hit is een dictionary met annotaties van Uniprot, 
    KEGG en nucleotide databases.

    Parameters:
    data (list): Een lijst van verrijkte BLAST hits.

    Returns:
    inserted (int): Het aantal records dat succesvol is opgeslagen in de database.

    """
    conn_string = os.environ.get("DATABASE_URL")
    if not conn_string:
        raise RuntimeError("Zet de DATABASE_URL environment variable eerst.")

    conn = psycopg2.connect(conn_string)
    inserted = 0
    with conn:
        with conn.cursor() as cursor:
            for hit in data:
                gene_id = _store_gene(hit, cursor)
                protein_id = _store_protein(hit, gene_id, cursor)
                _store_go_data(hit, protein_id, cursor)
                _store_pathways(hit, protein_id, cursor)
                inserted += 1
    return inserted


def _store_gene(hit, cursor):
    """
    Sla de geninformatie op in de database. Als het gen al bestaat, 
    wordt het bijgewerkt met nieuwe informatie.

    Parameters:
    hit (dict): Een verrijkte BLAST hit met geninformatie.
    cursor (psycopg2.cursor): De database cursor om de query's uit te voeren.

    Returns:
    gene_id (str): Het ID van het opgeslagen gen.
    """
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
    """
    Sla de eiwitinformatie op in de database. Als het eiwit al bestaat,
    wordt het bijgewerkt met nieuwe informatie.

    Parameters:
    hit (dict): Een verrijkte BLAST hit met eiwitinformatie.
    gene_id (str): Het ID van het bijbehorende gen.
    cursor (psycopg2.cursor): De database cursor om de query's uit te voeren.

    Returns:
    protein_id (str): Het ID van het opgeslagen eiwit.
    """
    cursor.execute(
        """
        INSERT INTO eiwit (
            eiwit_id, gen_id, naam, uniprot_beschrijving,
            uniprot_organisme, aminozuursequentie, blast_id,
            onderwerp_lengte, query_id, query_titel, query_lengte, evalue,
            bit_score, score, identiteiten, positieven, uitlijningslengte,
            gaten, query_start, query_eind, query_frame, onderwerp_start,
            onderwerp_eind, onderwerp_frame, query_sequentie, match_midline,
            onderwerp_sequentie
        )
        VALUES (
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s
        )
        RETURNING id
        """,
        (
            hit["accession"],
            gene_id,
            hit.get("uniprot_entry_name"),
            hit.get("uniprot_description"),
            hit.get("uniprot_organism"),
            hit.get("uniprot_sequence"),
            hit.get("blast_id"),
            hit.get("length"),
            hit.get("query_id"),
            hit.get("query_title"),
            hit.get("query_length"),
            hit.get("evalue"),
            hit.get("bit_score"),
            hit.get("score"),
            hit.get("identities"),
            hit.get("positives"),
            hit.get("align_length"),
            hit.get("gaps"),
            hit.get("query_start"),
            hit.get("query_end"),
            hit.get("query_frame"),
            hit.get("sbjct_start"),
            hit.get("sbjct_end"),
            hit.get("sbjct_frame"),
            hit.get("query_seq"),
            hit.get("match_midline"),
            hit.get("sbjct_seq"),
        ),
    )
    return cursor.fetchone()[0]


def _store_go_data(hit, protein_id, cursor):
    """
    Sla de GO annotaties op in de database. Als de GO annotatie al bestaat,
    wordt deze bijgewerkt met nieuwe informatie.

    Parameters:
    hit (dict): Een verrijkte BLAST hit met GO annotaties.
    protein_id (str): Het ID van het bijbehorende eiwit.
    cursor (psycopg2.cursor): De database cursor om de query's uit te voeren.
    """
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


def _store_pathways(hit, protein_id, cursor):
    """
    Sla de KEGG pathway annotaties op in de database. Als de pathway al bestaat,
    wordt deze bijgewerkt met nieuwe informatie.

    Parameters:
    hit (dict): Een verrijkte BLAST hit met KEGG pathway annotaties.
    protein_id (str): Het ID van het bijbehorende eiwit.
    cursor (psycopg2.cursor): De database cursor om de query's uit te voeren.
    """
    kegg_data = hit.get("kegg_data") or {}
    for kegg_entry in kegg_data.values():
        for pathway_entry in kegg_entry.get("PATHWAY") or []:
            pathway_kegg_id, _, pathway_name = pathway_entry.partition(" ")
            pathway_name = pathway_name.strip() or None
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
                (pathway_kegg_id, pathway_name, None, None),
            )
            cursor.execute(
                """
                INSERT INTO eiwit_pathway (eiwit_id, pathway_id)
                VALUES (%s, %s)
                ON CONFLICT (eiwit_id, pathway_id) DO NOTHING
                """,
                (protein_id, pathway_kegg_id),
            )
