from database_connectie import maak_database_connectie


def database_opzetten(conn):
    """Vervang de applicatietabellen door een lege database-indeling."""
    with conn:
        with conn.cursor() as cursor:
            for tabel in (
                "eiwit_functie",
                "eiwit_pathway",
                "blast_hsp",
                "isomeer",
                "eiwit",
                "functie",
                "pathway",
                "gen",
            ):
                cursor.execute(f"DROP TABLE IF EXISTS {tabel}")

            for tabel in (
                gen_aanmaken(),
                functie_aanmaken(),
                pathway_aanmaken(),
                isomeer_aanmaken(),
                eiwit_aanmaken(),
                eiwit_functie_aanmaken(),
                eiwit_pathway_aanmaken(),
            ):
                cursor.execute(tabel)
        conn.commit()
    


def tabellen_verwijderen(cursor):
    verwijder = [
        "Gen",
        "Functie",
        "Pathway",
        "Isomeer",
        "Eiwit",
        "Eiwit_Functie",
        "Eiwit_Pathway"
    ]

    for verwijder_tabel in verwijder:
        cursor.execute(f"DROP TABLE IF EXISTS {verwijder_tabel} CASCADE")


def gen_aanmaken():
    return """CREATE TABLE IF NOT EXISTS Gen (
        Gen_ID VARCHAR(50) PRIMARY KEY,
        Naam TEXT,
        Organisme TEXT,
        Ensembl_IDs TEXT[] NOT NULL DEFAULT '{}',
        Sequentie TEXT,
        Startpositie INT,
        Eindpositie INT,
        Beschrijving TEXT,
        Strand INT
    )"""


def functie_aanmaken():
    return """CREATE TABLE IF NOT EXISTS Functie (
        Go_ID VARCHAR(50) PRIMARY KEY,
        Go_Term TEXT,
        Go_Type VARCHAR(50)
    )"""


def pathway_aanmaken():
    return """CREATE TABLE IF NOT EXISTS Pathway (
        Pathway_Kegg_ID VARCHAR(50) PRIMARY KEY,
        Naam TEXT,
        Type_Pathway TEXT,
        Beschrijving TEXT
    )"""


def isomeer_aanmaken():
    return """CREATE TABLE IF NOT EXISTS Isomeer (
        Isomeer_ID TEXT PRIMARY KEY,
        Gen_ID VARCHAR(50),
        Sequentie TEXT,
        Splicings_variant TEXT,
        Codeert BOOLEAN,
        FOREIGN KEY (Gen_ID)
            REFERENCES Gen(Gen_ID)
    )"""


def eiwit_aanmaken():
    return """CREATE TABLE IF NOT EXISTS Eiwit (
        ID SERIAL PRIMARY KEY,
        Eiwit_ID VARCHAR(50),
        Gen_ID VARCHAR(50),
        Naam TEXT,
        Uniprot_Beschrijving TEXT,
        Uniprot_Organisme TEXT,
        Aminozuursequentie TEXT,
        Blast_ID TEXT,
        Onderwerp_lengte INT,
        Query_ID TEXT,
        Query_titel TEXT,
        Query_lengte INT,
        Evalue DOUBLE PRECISION,
        Bit_score DOUBLE PRECISION,
        Score INT,
        Identiteiten INT,
        Positieven INT,
        Uitlijningslengte INT,
        Gaten INT,
        Query_start INT,
        Query_eind INT,
        Query_frame INT,
        Onderwerp_start INT,
        Onderwerp_eind INT,
        Onderwerp_frame INT,
        Query_sequentie TEXT,
        Match_midline TEXT,
        Onderwerp_sequentie TEXT,
        FOREIGN KEY (Gen_ID)
            REFERENCES Gen(Gen_ID)
    )"""


def eiwit_functie_aanmaken():
    return """CREATE TABLE IF NOT EXISTS Eiwit_Functie (
        ID INT,
        Go_ID VARCHAR(50),
        PRIMARY KEY (ID, Go_ID),
        FOREIGN KEY (ID)
            REFERENCES Eiwit(ID),
        FOREIGN KEY (Go_ID)
            REFERENCES Functie(Go_ID)
    )"""


def eiwit_pathway_aanmaken():
    return """CREATE TABLE IF NOT EXISTS Eiwit_Pathway (
        Eiwit_ID INT,
        Pathway_ID VARCHAR(50),
        PRIMARY KEY (Eiwit_ID, Pathway_ID),
        FOREIGN KEY (Eiwit_ID)
        REFERENCES Eiwit(ID),
        FOREIGN KEY (Pathway_ID)
            REFERENCES Pathway(Pathway_Kegg_ID)
    )"""


if __name__ == "__main__":
    conn = maak_database_connectie()
    try:
        database_opzetten(conn)
    finally:
        conn.close()
