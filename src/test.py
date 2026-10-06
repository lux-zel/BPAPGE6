import json

def fasta_uitlezen(fasta):

    eiwitten = {}

    with open(fasta, "r") as bestand:

        accession = None
        naam = None
        gen_id = None
        aminozuur_sequentie = ""

        for regel in bestand:

            regel = regel.strip()

            # Controleert of het een header is
            if regel.startswith(">"):

                # Vorig eiwit opslaan
                if accession is not None:
                    eiwitten[accession] = {
                        "gen_id": gen_id,
                        "naam": naam,
                        "sequentie": aminozuur_sequentie
                    }

                #header uitlezen
                delen = regel.split("|")
                accession = delen[1]

                naam = regel.split(" ", 1)[1].split(" OS=")[0]
                if "GN=" in regel:
                    gen_id = regel.split("GN=")[1].split()[0]
                else:
                    gen_id = None

                aminozuur_sequentie = ""
            else:
                aminozuur_sequentie += regel



        if accession is not None:
            eiwitten[accession] = {
                "gen_id": gen_id,
                "naam": naam,
                "sequentie": aminozuur_sequentie
            }

    return eiwitten

def ophalen():

    with open('output/output_parse.json', 'r') as bestand:
        data = json.load(bestand)

    fasta_data = fasta_uitlezen("proteome.fasta")

    #ontbreekt = [a for a in data if a not in fasta_data]
    #print(len(data), "eiwitten in JSON")
    #print(len(ontbreekt), "daarvan niet in FASTA")

    eiwitten = []
    functies = {}
    eiwit_functies = []
    pathways = set()
    eiwit_pathways = []
    embl_ids = set()
    eiwit_embl = []

    for accession, protein_data in data.items():
        eiwit = fasta_data.get(accession)
        print("Eiwit ID:", accession)

        if eiwit is not None:
            print(' ')
            print("Gen ID:", eiwit["gen_id"])
            print("Naam:", eiwit["naam"])
            print("Aminozuursequentie:", eiwit["sequentie"])

            eiwitten.append((accession, eiwit["gen_id"], eiwit["naam"], eiwit["sequentie"]))

        else:
            print('Geen FASTA gegevens')
            continue

        # Gegevens
        go_ids = protein_data["GO_id"]
        go_terms = protein_data["GO_term"]
        go_types = protein_data["GO_type"]
        embl_ids_data = protein_data["embl_ids"]
        kegg_ids = protein_data["kegg_ids"]



        # GO data
        for i in range(len(go_ids)):
            if go_ids[i] is None:
                continue

            print("GO ID:", go_ids[i])
            print("GO Term:", go_terms[i])
            print("GO Type:", go_types[i])

            functies[go_ids[i]] =  (go_terms[i],go_types[i])
            eiwit_functies.append((accession,go_ids[i]))

        # Embl data
        for embl_id in embl_ids_data:

            if embl_id is None:
                continue

            print("EMBL ID:", embl_id)
            embl_ids.add(embl_id)
            eiwit_embl.append((accession, embl_id))

        # KEGG ID
        for kegg_id in kegg_ids:
            if kegg_id is None:
                continue

            print("KEGG ID:", kegg_id)
            pathways.add(kegg_id)
            eiwit_pathways.append((accession, kegg_id))


    return eiwitten, functies, eiwit_functies, pathways, eiwit_pathways, embl_ids, eiwit_embl

eiwitten, functies, eiwit_functies, pathways, eiwit_pathways, embl_ids, eiwit_embl = ophalen()




