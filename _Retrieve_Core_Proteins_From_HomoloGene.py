import os
import csv
import xml.etree.ElementTree as ET
import requests
import pandas as pd


def RetrieveCoreProteinsHomoloGene(queryName, inputPath, api):
    ### Declaring Global Variables here ###
    capitalizeQuery = queryName.title()

    db = "HomoloGene"
    queryOrganism = "Human"
    NCBIGeneList = []
    # 1: "Human", 2: "Cow", 3: "Mouse", 4: "Zebrafish", 5: "Fly", 6: "Worm", 7: "Yeast", 8: "Frog"
    organismTaxIDDict = {1: 9606, 2: 9913, 3: 10090, 4: 7955,
                         5: 7227, 6: 6239, 7: 4932, 8: 8364}

    ### Paths/Input/Output Files
    outputPathExcel = "{inputPath}_HomoloGeneSets/".format(inputPath=inputPath)

    inputFile = "{inputPath}_FullProteinSets/{studying}_UniProt_Full.csv".format(studying=capitalizeQuery, inputPath=inputPath)
    outputHomoloGene = "{output}_{study}_Core_HomoloGene_List.csv".format(study=capitalizeQuery, output=outputPathExcel)

    ## Functions ##
    def RemoveExistingFiles(file):
        for files in file:
            if os.path.exists(files):
                os.remove(files)

    def WriteToCsv(path, row, fileMode):
        with open(path, fileMode) as fd:
            wr = csv.writer(fd, quoting=csv.QUOTE_ALL)
            wr.writerow(row)

    def OpenCSV():
        columnList = ['Entry', 'Gene', 'Organism']
        df = pd.read_csv(inputFile, usecols=columnList)

        # Get only "Organism" and sort by Gene name
        df1 = df.loc[df['Organism'] == queryOrganism].sort_values("Gene")

        # Reformatting Entries if not already done
        # Sort for just human lists here
        df1 = df1[df1['Entry'].notna()]
        df1 = df1[df1['Gene'].notna()]
        df1['Gene'] = df1['Gene'].str.upper()
        df1 = df1[df1['Organism'].notna()]
        df1.drop_duplicates(subset="Gene",
                            keep="first", inplace=True)

        gene = df1['Gene'].tolist()

        return gene

    def HomoloGeneCheck(lst):
        for query in lst:
            print("Checking HomoloGene for... ", query)
            urlESearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={db}&term={" \
                         "query}{api}".format(db=db, query=query, api=api)
            loopNum = 0
            while loopNum < 3:
                try:
                    # print(urlESearch)
                    r = requests.get(url=urlESearch, stream=False)
                    root = ET.fromstring(r.content)
                    # print(r.content) #prints the xml content
                    break
                except:
                    loopNum += 1

            idList = []
            for child in root.iter('Id'):
                idList.append(child.text)
            # print("idList: ", idList)

            for uid in idList:
                urlESummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db={db}&id={" \
                              "uid}&version=2.0{api}".format(
                    db=db, uid=uid, api=api)
                # print(urlESummary)
                r = requests.get(url=urlESummary, stream=False)
                root = ET.fromstring(r.content)

                taxIdList = []
                symblList = []
                resultingSymblList = []

                for taxName in root.iter('TaxId'):
                    taxIdList.append(int(taxName.text))
                for symbl in root.iter('Symbol'):
                    symblList.append(symbl.text)

                for i, taxID in enumerate(taxIdList):
                    if taxID in organismTaxIDDict.values() and taxID not in taxIdList[:i]:
                        resultingSymblList.append(symblList[i])

                if len(resultingSymblList) == len(organismTaxIDDict):
                    print(query, "found", '\n')
                    WriteToCsv(outputHomoloGene, resultingSymblList, 'a')
                    break

    ##----------- Execution ----------- ##

    try:
        files = [outputHomoloGene]
        RemoveExistingFiles(files)
    except:
        print("No Existing files exist to delete.")

    WriteToCsv(outputHomoloGene, ["Human", "Cow", "Mouse", "Zebrafish", "Fly", "Worm", "Yeast", "Frog"], 'a')
    geneList = OpenCSV()
    HomoloGeneCheck(geneList)
