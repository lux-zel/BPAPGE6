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
            "gen"
        ]

        for verwijder_tabel in verwijder:
            cursor.execute(f"DROP TABLE IF EXISTS {verwijder_tabel} CASCADE")

    # TABELLEN AANMAKEN:
    def gen_aanmaken():
        return """CREATE TABLE IF NOT EXISTS gen 
        (Gen_ID VARCHAR(50) PRIMARY KEY,
        Naam VARCHAR(50),
        Sequentie VARCHAR(500),
        Startpositie INT,
        Eindpositie INT,
        Beschrijving VARCHAR(250),
        Strand VARCHAR(1)
        )"""


    if __name__ == '__main__':
        # Alle tabellen verwijderen
        tabellen_verwijderen()

        # File inlezen

        # Alle tabellen weer opnieuw aanmaken
        tabellen = [
            gen_aanmaken(),

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
