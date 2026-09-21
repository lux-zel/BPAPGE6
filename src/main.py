
import sys

def read_fasta(filename):
    with open(filename, 'r') as file:
        reads = file.read().strip().split('>')
    return reads







def main():

    if len(sys.argv) < 3:
        print("roep het programma met: python3 main.py <FASTA 1> <FASTA 2>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]



if __name__ == '__main__':
    main()
