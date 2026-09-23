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



except psycopg2.Error as error:
    print("Error", error)

    if cursor:
        cursor.close()

    if conn:
        conn.close()
        print('Connection closed')
