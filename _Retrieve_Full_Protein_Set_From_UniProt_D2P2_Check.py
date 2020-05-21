import csv
import json
import pandas as pd
import requests
import os
import glob
from urllib.request import urlopen

pd.set_option("display.max_rows", None, "display.max_columns", None)


# # # # # # Purpose # # # # # # # # # # # # # # # # # # # # # # # #
# Query and get the Uniprot Accession ID from UniProt Website,    #
# Check against D2P2 Database for entries                         #
# Store the information as downloaded excel file, convert to csv  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def RetrieveFullProtSetUniProtD2P2Check(query, inputPath):
    # Global variables
    capitalizeQuery = query.title()
    organismTaxIDDict_Revised = {1: "fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22",
                                 2: "fil=organism%3A%22Bos+taurus+%28Bovine%29+%5B9913%5D%22",
                                 3: "fil=organism%3A%22Mus+musculus+%28Mouse%29+%5B10090%5D%22",
                                 4: "fil=organism%3A%22Danio+rerio+%28Zebrafish%29+%5B7955%5D%22",
                                 5: "fil=organism%3A%22Drosophila+melanogaster+%28Fruit+Fly%29+%5B7227%5D%22",
                                 6: "fil=organism%3A%22Caenorhabditis+elegans+%5B6239%5D%22",
                                 7: "fil=organism%3A%22Saccharomyces+cerevisiae+%28Yeast%29+%5B559292%5D%22",
                                 8: "fil=organism%3A%22Xenopus+tropicalis+%28Frog%29+%5B8364%5D%22"}
    organismNameDict = {1: "Human", 2: "Cow", 3: "Mouse", 4: "Zebrafish", 5: "Fly", 6: "Worm", 7: "Yeast", 8: "Frog"}

    # Output Files/Directory
    path = "{inputPath}_FullProteinSets/".format(inputPath=inputPath)
    outputFilePath = "{path}/{name}_UniProt_Full.csv".format(path=path, name=capitalizeQuery)
    outputFailedPath = "{path}/{name}_UniProt_D2P2_Check_Failed.csv".format(path=path, name=capitalizeQuery)

    def UniProt2CSV():
        for i in range(0, len(organismTaxIDDict_Revised)):
            # Format = html | tab | xls | fasta | gff | txt | xml | rdf | list | rss
            URL = "https://www.uniprot.org/uniprot/?query={query}&{organism}&columns=id,genes(PREFERRED),genes(ORF),organism,reviewed,comment(SUBCELLULAR LOCATION)&format=xls".format(
                query=query, organism=organismTaxIDDict_Revised[i + 1])

            print("Retrieving Input from UniProt for organism: ", organismNameDict[i + 1], "-> ", end="")
            loopNum = 0
            while True:
                try:
                    r = requests.get(URL, allow_redirects=True)

                    open("{path}/{query}_Partial_{i}.xlsx".format(path=path, query=query, i=i + 1), 'wb').write(
                        r.content)
                    break
                except:
                    loopNum += 1
                    print('\n')
                    print("Writing Failed. Trying again...: ", loopNum, '\n')

            print("Done.")

    def CheckD2P2():
        # Merging the Excel files previously created
        mergedDf = pd.DataFrame()
        for f in glob.glob("{path}*Partial_*.xlsx".format(path=path)):
            df = pd.read_excel(f)
            mergedDf = mergedDf.append(df, ignore_index=True)
            os.remove(f)

        # Renaming the Columns/Elements in Pandas
        mergedDf.rename(columns={'Gene names  (primary )': 'Gene'}, inplace=True)
        mergedDf.rename(columns={'Gene names  (ORF )': 'ORF'}, inplace=True)

        mergedDf['Organism'] = mergedDf['Organism'].replace(
            {"Homo sapiens (Human)": "Human", "Bos taurus (Bovine)": "Cow",
             "Mus musculus (Mouse)": "Mouse",
             "Danio rerio (Zebrafish) (Brachydanio rerio)": "Zebrafish",
             "Drosophila melanogaster (Fruit fly)": "Fly",
             "Caenorhabditis elegans": "Worm",
             "Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)": "Yeast",
             "Xenopus tropicalis (Western clawed frog) (Silurana tropicalis)": "Frog"})


        # Remove entries that don't have nucleolus in the subcellular location
        locDf = mergedDf[mergedDf['Subcellular location [CC]'].str.contains(query, na=False)]

        # Remove entries that aren't reviewed
        reviewedDf = locDf[~locDf['Status'].str.contains("unreviewed")]
        # Remove entries not in nucleolus

        # Add Orf positions as names for some of the organisms
        geneList = reviewedDf['Gene'].tolist()
        entryList = reviewedDf['Entry'].tolist()
        orgnList = reviewedDf['Organism'].tolist()

        for i, seqid in enumerate(entryList):
            requestURL = 'http://d2p2.pro/api/seqid/["{seqid}"]'.format(seqid=seqid)
            response = json.loads(urlopen(requestURL).read())

            try:
                # Original: [0][2]
                for consensus in response["{seqid}".format(seqid=seqid)][0][2]["disorder"]["consensus"]:
                    checkedEntryList = []
                    checkedEntryList.append(seqid)
                checkedEntryList.append(geneList[i])
                checkedEntryList.append(orgnList[i])
                WriteToCSV(outputFilePath, checkedEntryList, 'a')

            except:
                failedEntryList = []
                failedEntryList.append(seqid)
                failedEntryList.append(geneList[i])
                failedEntryList.append(orgnList[i])
                print(seqid, "Failed")

                WriteToCSV(outputFailedPath, failedEntryList, 'a')

    def RemoveDuplicates():
        columnList = ['Entry', 'Gene', 'Organism']
        df = pd.read_csv(outputFilePath, usecols=columnList)
        os.remove(outputFilePath)
        df.drop_duplicates(['Gene', 'Organism'], inplace=True, keep="first")

        df.to_csv(outputFilePath, columns=['Entry', 'Gene', 'Organism'], mode='w', index=False)

    def WriteToCSV(pathOutput, row, fileMode):
        with open(pathOutput, fileMode) as fd:
            wr = csv.writer(fd, quoting=csv.QUOTE_ALL)
            wr.writerow(row)

    # #----------- Execution ----------- # #
    # Remove OutputFilePath and Failed Path if already exists
    try:
        os.remove(outputFilePath)
    except OSError:
        pass

    try:
        os.remove(outputFailedPath)
    except OSError:
        pass

    WriteToCSV(outputFilePath, ["Entry", "Gene", "Organism"], 'a')
    WriteToCSV(outputFailedPath, ["Entry", "Gene", "Organism"], 'a')

    UniProt2CSV()
    CheckD2P2()

    RemoveDuplicates()
