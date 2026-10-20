def database_vullen(conn, hits):
    """Fill all tables from the enriched BLAST hits."""
    with conn:
        with conn.cursor() as cur:
            for hit in hits:
                gen_id = _gen_toevoegen(cur, hit)
                eiwit_id = _eiwit_toevoegen(cur, hit, gen_id)
                _isomeer_toevoegen(cur, hit, gen_id)
                _functie_toevoegen(cur, hit, eiwit_id)
                _pathway_toevoegen(cur, hit, eiwit_id)


def _gen_toevoegen(cur, hit):
    """Insert gene if it has a name; return Gen_ID or None."""
    gen_id = hit.get("gene_primary")
    if not gen_id:
        return None
    cur.execute(
        """INSERT INTO Gen (Gen_ID, Naam, Organisme)
           VALUES (%s, %s, %s)
           ON CONFLICT (Gen_ID) DO NOTHING""",
        (gen_id, gen_id, hit.get("organism_name")),
    )
    return gen_id


def _eiwit_toevoegen(cur, hit, gen_id):
    """Insert protein + its BLAST HSP; return Eiwit.ID (serial)."""
    cur.execute(
        """INSERT INTO Eiwit (
               Eiwit_ID, Gen_ID, Naam,
               Uniprot_Organisme, Aminozuursequentie, Lengte,
               Query_ID, Evalue, Bit_score, Score, Identiteiten,
               Uitlijningslengte, Gaps,
               Query_start, Query_eind, Query_frame,
               Onderwerp_start, Onderwerp_eind
           ) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
           RETURNING ID""",
        (
            hit["accession"], gen_id,
            hit.get("protein_name"),
            hit.get("organism_name"), hit.get("sequence"), hit.get("length"),
            hit.get("query_id"), hit.get("evalue"), hit.get("bitscore"),
            hit.get("score"), hit.get("identities"), hit.get("align_length"),
            hit.get("gaps"),
            hit.get("query_start"), hit.get("query_end"), hit.get("query_frame"),
            hit.get("subject_start"), hit.get("subject_end")
        ),
    )
    return cur.fetchone()[0]


def _isomeer_toevoegen(cur, hit, gen_id):
    """Insert each isoform, linked to its gene (skip if no gene)."""
    if not gen_id:
        return
    for iso_id in hit.get("isoforms") or []:
        cur.execute(
            """INSERT INTO Isomeer (Isomeer_ID, Gen_ID)
               VALUES (%s, %s)
               ON CONFLICT (Isomeer_ID) DO NOTHING""",
            (iso_id, gen_id),
        )


def _functie_toevoegen(cur, hit, eiwit_id):
    """Insert GO terms + link rows."""
    for go in hit.get("go_terms") or []:
        cur.execute(
            """INSERT INTO Functie (Go_ID, Go_Term, Go_Type)
               VALUES (%s, %s, %s)
               ON CONFLICT (Go_ID) DO NOTHING""",
            (go["go_id"], go["term"], go["aspect"]),
        )
        cur.execute(
            """INSERT INTO Eiwit_Functie (Eiwit_ID, Go_ID)
               VALUES (%s, %s)
               ON CONFLICT DO NOTHING""",
            (eiwit_id, go["go_id"]),
        )


def _pathway_toevoegen(cur, hit, eiwit_id):
    for path_id, path_name in (hit.get("kegg_pathways") or {}).items():
        cur.execute(
            """INSERT INTO Pathway (Pathway_Kegg_ID, Naam)
               VALUES (%s, %s) ON CONFLICT DO NOTHING""",
            (path_id, path_name),
        )
        cur.execute(
            """INSERT INTO Eiwit_Pathway (Eiwit_ID, Pathway_ID)
               VALUES (%s, %s) ON CONFLICT DO NOTHING""",
            (eiwit_id, path_id),
        )