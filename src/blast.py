"""
dit bestand is een helper python bestand dat een blast opdracht uitvoert
van een query tegen een subject. 
"""

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO



def blast_pairwise(query, subject):
    """
    Deze function maakt gebruikt van Biopython's wrapper rond de BLAST+ binaries
    Het voert een local paiwise BLAST uit met -subject in plaats van een database.
    
    Parameters:
    query (string):     De fasta sequentie van de sequentie die je Queried.
    subject (string):   De fasta sequentie waartegen je Queried.
    seq_type (string):  Het type allignment ("nucl" voor blastn, "prot" voor blastp).
    
    Returns:
    hits (dictionary):  Een dictionary met de gevonden hits van de BLAST.
    
    """
    
    SeqIO.write(SeqIO.SeqRecord(Seq(query), id="query"), "query_fa", "fasta")
    SeqIO.write(SeqIO.SeqRecord(Seq(subject), id="subject"), "subject_fa", "fasta")
    
    cline = NcbiblastnCommandline(
            query=query_fa,
            subject=subject_fa,
            outfmt=5,
            out=xml_out,
            evalue=evalue,
            task='blastn')

