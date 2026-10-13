from database_connectie import maak_database_connectie


def database_opzetten(conn):
    """Maak ontbrekende tabellen aan zonder bestaande tabellen te verwijderen."""
    with conn:
        with conn.cursor() as cursor:
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
            cursor.execute(
                "ALTER TABLE eiwit ADD COLUMN IF NOT EXISTS hit_data "
                "JSONB NOT NULL DEFAULT '[]'::jsonb"
            )
            cursor.execute("DROP TABLE IF EXISTS blast_hsp")


def gen_aanmaken():
    return """CREATE TABLE IF NOT EXISTS Gen (
        Gen_ID VARCHAR(50) PRIMARY KEY,
        Naam VARCHAR(50),
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
        Isomeer_ID VARCHAR(50) PRIMARY KEY,
        Gen_ID VARCHAR(50),
        Sequentie TEXT,
        Splicings_variant VARCHAR(50),
        Codeert BOOLEAN,
        FOREIGN KEY (Gen_ID)
            REFERENCES Gen(Gen_ID)
    )"""


def eiwit_aanmaken():
    return """CREATE TABLE IF NOT EXISTS Eiwit (
        ID SERIAL PRIMARY KEY, 
        Eiwit_ID VARCHAR(50) UNIQUE,
        Gen_ID VARCHAR(50),
        Naam TEXT,
        Aminozuursequentie TEXT,
        Hit_Data JSONB NOT NULL DEFAULT '[]'::jsonb,
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
