def fill_database(data, conn):
    """Bewaar de beste BLAST-hit per eiwit en de bijbehorende annotaties."""
    inserted = 0
    with conn:
        with conn.cursor() as cursor:
            for hit in _best_hsp_per_protein(data):
                protein_id = _store_protein(hit, cursor)
                _store_go_data(hit, protein_id, cursor)
                _store_pathways(hit, protein_id, cursor)
                inserted += 1
    return inserted


def _best_hsp_per_protein(data):
    best_hits = {}
    for hit in data:
        accession = hit["accession"]
        current = best_hits.get(accession)
        if current is None or (hit["bit_score"], -hit["evalue"]) > (
            current["bit_score"], -current["evalue"]
        ):
            best_hits[accession] = hit
    return best_hits.values()


def _store_protein(hit, cursor):
    gene_ids = hit.get("gene_id") or []
    gene_data = hit.get("gene_data") or {}
    gene_id = (
        (gene_ids[0] if gene_ids else None)
        or hit.get("gene_locus_tag")
        or gene_data.get("locus_tag")
    )

    if gene_id is not None:
        cursor.execute(
            """
            INSERT INTO gen (
                gen_id, sequentie, startpositie, eindpositie, beschrijving, strand
            )
            VALUES (%s, %s, %s, %s, %s, %s)
            ON CONFLICT (gen_id) DO UPDATE SET
                sequentie = COALESCE(EXCLUDED.sequentie, gen.sequentie),
                startpositie = COALESCE(EXCLUDED.startpositie, gen.startpositie),
                eindpositie = COALESCE(EXCLUDED.eindpositie, gen.eindpositie),
                beschrijving = COALESCE(EXCLUDED.beschrijving, gen.beschrijving),
                strand = COALESCE(EXCLUDED.strand, gen.strand)
            """,
            (
                gene_id,
                gene_data.get("sequence"),
                gene_data.get("start"),
                gene_data.get("end"),
                gene_data.get("description") or hit.get("uniprot_description"),
                gene_data.get("strand"),
            ),
        )
        _store_isomers(hit, gene_id, cursor)

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


def _store_isomers(hit, gene_id, cursor):
    if gene_id is None:
        return

    for isomer in hit.get("isomers") or []:
        cursor.execute(
            """
            INSERT INTO isomeer (
                isomeer_id, gen_id, sequentie, splicings_variant, codeert
            )
            VALUES (%s, %s, %s, %s, %s)
            ON CONFLICT (isomeer_id) DO UPDATE SET
                gen_id = EXCLUDED.gen_id,
                sequentie = COALESCE(EXCLUDED.sequentie, isomeer.sequentie),
                splicings_variant = COALESCE(
                    EXCLUDED.splicings_variant, isomeer.splicings_variant
                ),
                codeert = COALESCE(EXCLUDED.codeert, isomeer.codeert)
            """,
            (
                isomer["isomer_id"],
                gene_id,
                isomer.get("sequence"),
                isomer.get("splicing_variant"),
                isomer.get("codes"),
            ),
        )


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
