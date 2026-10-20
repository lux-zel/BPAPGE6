def database_opzetten(conn):
    """Drop and recreate all tables."""
    with conn:
        with conn.cursor() as cursor:
            tabellen_verwijderen(cursor)
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


def tabellen_verwijderen(cursor):
    for tabel in (
        "Eiwit_Functie", "Eiwit_Pathway", "Isomeer",
        "Eiwit", "Functie", "Pathway", "Gen",
    ):
        cursor.execute(f"DROP TABLE IF EXISTS {tabel} CASCADE")


def gen_aanmaken():
    return """CREATE TABLE Gen (
        Gen_ID VARCHAR(50) PRIMARY KEY,
        Naam TEXT,
        Organisme TEXT
    )"""


def functie_aanmaken():
    return """CREATE TABLE Functie (
        Go_ID VARCHAR(50) PRIMARY KEY,
        Go_Term TEXT,
        Go_Type VARCHAR(50)
    )"""


def pathway_aanmaken():
    return """CREATE TABLE Pathway (
        Pathway_Kegg_ID VARCHAR(50) PRIMARY KEY,
        Naam TEXT
    )"""


def isomeer_aanmaken():
    return """CREATE TABLE Isomeer (
        Isomeer_ID TEXT PRIMARY KEY,
        Gen_ID VARCHAR(50),
        FOREIGN KEY (Gen_ID) REFERENCES Gen(Gen_ID)
    )"""


def eiwit_aanmaken():
    return """CREATE TABLE Eiwit (
        ID SERIAL PRIMARY KEY,
        Eiwit_ID VARCHAR(50),
        Gen_ID VARCHAR(50),
        Naam TEXT,
        Uniprot_Organisme TEXT,
        Aminozuursequentie TEXT,
        Lengte INT,
        Query_ID TEXT,
        Evalue DOUBLE PRECISION,
        Bit_score DOUBLE PRECISION,
        Score INT,
        Identiteiten INT,
        Uitlijningslengte INT,
        Gaps INT,
        Query_start INT,
        Query_eind INT,
        Query_frame INT,
        Onderwerp_start INT,
        Onderwerp_eind INT,
        FOREIGN KEY (Gen_ID) REFERENCES Gen(Gen_ID)
    )"""


def eiwit_functie_aanmaken():
    return """CREATE TABLE Eiwit_Functie (
        Eiwit_ID INT,
        Go_ID VARCHAR(50),
        PRIMARY KEY (Eiwit_ID, Go_ID),
        FOREIGN KEY (Eiwit_ID) REFERENCES Eiwit(ID),
        FOREIGN KEY (Go_ID) REFERENCES Functie(Go_ID)
    )"""


def eiwit_pathway_aanmaken():
    return """CREATE TABLE Eiwit_Pathway (
        Eiwit_ID INT,
        Pathway_ID VARCHAR(50),
        PRIMARY KEY (Eiwit_ID, Pathway_ID),
        FOREIGN KEY (Eiwit_ID) REFERENCES Eiwit(ID),
        FOREIGN KEY (Pathway_ID) REFERENCES Pathway(Pathway_Kegg_ID)
    )"""