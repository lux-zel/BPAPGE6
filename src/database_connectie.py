from getpass import getpass

import psycopg2


def maak_database_connectie():
    """Maak een databaseconnectie en vraag het wachtwoord veilig op."""
    wachtwoord = getpass("Databasewachtwoord: ")
    return psycopg2.connect(
        host="145.97.18.240",
        dbname="bpapge6_db",
        user="bpapge6",
        password=wachtwoord,
    )
