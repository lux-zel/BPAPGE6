def database_vullen(conn, hits, isoforms=()):
    """Sla verrijkte BLAST-hits en FASTA-isoformen op in de database."""
    beste_hits = {}
    for hit in hits:
        accession = hit["accession"]
        huidige_hit = beste_hits.get(accession)
        huidige_score = huidige_hit.get("bitscore") if huidige_hit else None
        hit_score = hit.get("bitscore")
        if huidige_hit is None or (
            hit_score is not None
            and (huidige_score is None or hit_score > huidige_score)
        ):
            beste_hits[accession] = hit

    with conn:
        with conn.cursor() as cur:
            for hit in beste_hits.values():
                gen_id = _gen_toevoegen(cur, hit)
                eiwit_id = _eiwit_toevoegen(cur, hit, gen_id)
                _functie_toevoegen(cur, hit, eiwit_id)
                _pathway_toevoegen(cur, hit, eiwit_id)

            for isoform in isoforms:
                gen_id = _gen_toevoegen(
                    cur,
                    {
                        "gene_primary": isoform["gene_id"],
                        "organism_name": isoform["organism_name"],
                    },
                )
                _isomeer_toevoegen(cur, isoform, gen_id)


def _gen_toevoegen(cur, hit):
    """Voeg een gen toe en geef het gen-ID terug."""
    gen_id = hit.get("gene_primary")
    if not gen_id:
        return None
    cur.execute(
          """INSERT INTO Gen (Gen_ID)
              VALUES (%s)
           ON CONFLICT (Gen_ID) DO NOTHING""",
          (gen_id,),
    )
    return gen_id


def _eiwit_toevoegen(cur, hit, gen_id):
    """Voeg een eiwit en de BLAST-HSP toe en geef het database-ID terug."""
    cur.execute(
        """INSERT INTO Eiwit (
               Eiwit_ID, Gen_ID, Naam,
               Aminozuursequentie,
               Evalue, Bit_score, Gaps,
               Query_start, Query_eind,
               Onderwerp_start, Onderwerp_eind
           ) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
           RETURNING ID""",
        (
               hit["accession"],
               gen_id,
               hit.get("protein_name"),
               hit.get("sequence"),
               hit.get("evalue"),
               hit.get("bitscore"),
               hit.get("gaps"),
               hit.get("query_start"),
               hit.get("query_end"),
               hit.get("subject_start"),
               hit.get("subject_end"),
           ),
    )
    return cur.fetchone()[0]


def _isomeer_toevoegen(cur, isoform, gen_id):
    """Voeg een FASTA-isoformannotatie toe en koppel deze aan het gen."""
    isoform_name = isoform.get("isoform_name")
    if not gen_id or not isoform_name:
        return
    cur.execute(
        """INSERT INTO Isomeer (Isomeer_ID, Gen_ID, Naam)
           VALUES (%s, %s, %s)
           ON CONFLICT (Isomeer_ID) DO NOTHING""",
        (isoform["accession"], gen_id, isoform_name),
    )


def _functie_toevoegen(cur, hit, eiwit_id):
    """Voeg GO-termen toe en koppel ze aan het eiwit."""
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
    """Voeg KEGG-pathways toe en koppel ze aan het eiwit."""
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