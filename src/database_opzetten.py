def database_opzetten(conn):
    """Verwijder bestaande tabellen en maak een nieuw databaseschema aan."""
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
    """Verwijder alle applicatietabellen en bijbehorende beperkingen."""
    for tabel in (
        "Eiwit_Functie", "Eiwit_Pathway", "Isomeer",
        "Eiwit", "Functie", "Pathway", "Gen",
    ):
        cursor.execute(f"DROP TABLE IF EXISTS {tabel} CASCADE")


def gen_aanmaken():
    """Geef de SQL-opdracht voor de gentabel terug."""
    return """CREATE TABLE Gen (
        Gen_ID VARCHAR(50) PRIMARY KEY
    )"""


def functie_aanmaken():
    """Geef de SQL-opdracht voor de functietabel terug."""
    return """CREATE TABLE Functie (
        Go_ID VARCHAR(50) PRIMARY KEY,
        Go_Term TEXT,
        Go_Type VARCHAR(50)
    )"""


def pathway_aanmaken():
    """Geef de SQL-opdracht voor de pathwaytabel terug."""
    return """CREATE TABLE Pathway (
        Pathway_Kegg_ID VARCHAR(50) PRIMARY KEY,
        Naam TEXT
    )"""


def isomeer_aanmaken():
    """Geef de SQL-opdracht voor de isomeertabel terug."""
    return """CREATE TABLE Isomeer (
        Isomeer_ID TEXT PRIMARY KEY,
        Gen_ID VARCHAR(50),
        Naam TEXT,
        FOREIGN KEY (Gen_ID) REFERENCES Gen(Gen_ID)
    )"""


def eiwit_aanmaken():
    """Geef de SQL-opdracht voor de eiwittabel terug."""
    return """CREATE TABLE Eiwit (
        ID SERIAL PRIMARY KEY,
        Eiwit_ID VARCHAR(50) UNIQUE,
        Gen_ID VARCHAR(50),
        Naam TEXT,
        Aminozuursequentie TEXT,
        Evalue DOUBLE PRECISION,
        Bit_score DOUBLE PRECISION,
        Gaps INT,
        Query_start INT,
        Query_eind INT,
        Onderwerp_start INT,
        Onderwerp_eind INT,
        FOREIGN KEY (Gen_ID) REFERENCES Gen(Gen_ID)
    )"""


def eiwit_functie_aanmaken():
    """Geef de SQL-opdracht voor de eiwit-functiekoppeltabel terug."""
    return """CREATE TABLE Eiwit_Functie (
        Eiwit_ID INT,
        Go_ID VARCHAR(50),
        PRIMARY KEY (Eiwit_ID, Go_ID),
        FOREIGN KEY (Eiwit_ID) REFERENCES Eiwit(ID),
        FOREIGN KEY (Go_ID) REFERENCES Functie(Go_ID)
    )"""


def eiwit_pathway_aanmaken():
    """Geef de SQL-opdracht voor de eiwit-pathwaykoppeltabel terug."""
    return """CREATE TABLE Eiwit_Pathway (
        Eiwit_ID INT,
        Pathway_ID VARCHAR(50),
        PRIMARY KEY (Eiwit_ID, Pathway_ID),
        FOREIGN KEY (Eiwit_ID) REFERENCES Eiwit(ID),
        FOREIGN KEY (Pathway_ID) REFERENCES Pathway(Pathway_Kegg_ID)
    )"""