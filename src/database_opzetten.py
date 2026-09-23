import psycopg2


conn_string = """
host='145.97.18.240' dbname='bpapge6_db'
user='bpapge6' password='bpapge6'
"""

conn = None
cursor = None

try:
    conn = psycopg2.connect(conn_string)
    print('Connection open')
    cursor = conn.cursor()

    def tabellen_verwijderen():
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

    # TABELLEN AANMAKEN:
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
            Go_Term VARCHAR(50),
            Beschrijving TEXT
        )"""


    def pathway_aanmaken():
        return """CREATE TABLE IF NOT EXISTS Pathway (
            Pathway_Kegg_ID VARCHAR(50) PRIMARY KEY,
            Naam VARCHAR(50),
            Type_Pathway VARCHAR(50),
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
            Eiwit_ID VARCHAR(50),
            Gen_ID VARCHAR(50),
            Naam VARCHAR(50),
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


    if __name__ == '__main__':
        # Alle tabellen verwijderen
        tabellen_verwijderen()

        # File inlezen

        # Alle tabellen weer opnieuw aanmaken
        tabellen = [
            gen_aanmaken(),
            functie_aanmaken(),
            pathway_aanmaken(),
            isomeer_aanmaken(),
            eiwit_aanmaken(),
            eiwit_functie_aanmaken(),
            eiwit_pathway_aanmaken()
        ]

        for tabel in tabellen:
            cursor.execute(tabel)

        # Tabel vullen:

        # Connectie sluiten
        conn.commit()
        conn.close()
        print('Connection closed')

except psycopg2.Error as error:
    print("Error", error)

    if cursor:
        cursor.close()

    if conn:
        conn.close()
        print('Connection closed')
