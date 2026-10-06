import os

import psycopg2

def database_opzetten():
    """
    Deze functie maakt een connectie met de database en maakt de tabellen aan.
    """
    conn_string = os.environ.get("DATABASE_URL")
    if not conn_string:
        raise RuntimeError("Set the DATABASE_URL environment variable first.")

    conn = psycopg2.connect(conn_string)
    try:
        with conn.cursor() as cursor:
            tabellen_verwijderen(cursor)
            tabellen = [
                gen_aanmaken(),
                functie_aanmaken(),
                pathway_aanmaken(),
                isomeer_aanmaken(),
                eiwit_aanmaken(),
                eiwit_functie_aanmaken(),
                eiwit_pathway_aanmaken(),
                blast_hsp_aanmaken(),
            ]
            for tabel in tabellen:
                cursor.execute(tabel)
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def tabellen_verwijderen(cursor):
    verwijder = [
        "Gen",
        "Functie",
        "Pathway",
        "Isomeer",
        "Eiwit",
        "Eiwit_Functie",
        "Eiwit_Pathway",
        "Blast_HSP",
    ]

    for verwijder_tabel in verwijder:
        cursor.execute(f"DROP TABLE IF EXISTS {verwijder_tabel} CASCADE")


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


def blast_hsp_aanmaken():
    return """CREATE TABLE IF NOT EXISTS Blast_HSP (
        ID SERIAL PRIMARY KEY,
        Eiwit_ID INT NOT NULL,
        Blast_ID TEXT,
        Title TEXT,
        Length INT,
        Query_ID TEXT,
        Query_Title TEXT,
        Query_Length INT,
        Evalue DOUBLE PRECISION,
        Bit_Score DOUBLE PRECISION,
        Score DOUBLE PRECISION,
        Identities INT,
        Positives INT,
        Align_Length INT,
        Gaps INT,
        Query_Start INT,
        Query_End INT,
        Query_Frame INT,
        Sbjct_Start INT,
        Sbjct_End INT,
        Sbjct_Frame INT,
        Query_Seq TEXT,
        Match_Midline TEXT,
        Sbjct_Seq TEXT,
        FOREIGN KEY (Eiwit_ID)
            REFERENCES Eiwit(ID)
    )"""


if __name__ == "__main__":
    database_opzetten()
