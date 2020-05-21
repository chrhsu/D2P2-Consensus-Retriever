import shutil
import csv
import smtplib
import time
import requests
import pandas as pd
import xml.etree.ElementTree as ET
import os

# - 1) Need to write a function that opens Full [Insert Proteome] Human CSV - #
# - 2) Runs the HomoloGene Checker to get "Core" proteins - #
# - 3) Back checks with UniProt to get the "Core Protein IDs" - #
# - 4) Stores "UniProt ID, Gene Name, and Organism - #


### Declaring Global Variables here ###
studying = "Nucleolus"
db = "HomoloGene"
NCBIGeneList = []

# 1: "Human", 2: "Cow", 3: "Mouse", 4: "Zebrafish", 5: "Fly", 6: "Worm", 7: "Yeast", 8: "Frog"
organismTaxIDDict = {1: 9606, 2: 9913, 3: 10090, 4: 7955,
                     5: 7227, 6: 6239, 7: 4932, 8: 8364}
# Need to revise Yeast Dictionary for UniProt query
organismTaxIDDict_Revised = {1: 9606, 2: 9913, 3: 10090, 4: 7955,
                             5: 7227, 6: 6239, 7: 559292, 8: 8364}
organismNameDict = {1: "Human", 2: "Cow", 3: "Mouse", 4: "Zebrafish", 5: "Fly", 6: "Worm", 7: "Yeast", 8: "Frog"}
organismTaxIDNameDict = {9606: "Human", 9913: "Cow", 10090: "Mouse", 7955: "Zebrafish", 7227: "Fly", 6239: "Worm", 559292: "Yeast", 8364: "Frog"}

# Start running time
start_time = time.time()


## Input/Output Files Here ###

#PATHS
outputPathExcel = "_GUI/Core/"
if not os.path.exists('{path}{nameDir}/'.format(path=outputPathExcel, nameDir="Temp")):
    os.makedirs('{path}{nameDir}/'.format(path=outputPathExcel, nameDir="Temp"))
tempPath = "{path}{Temp}/".format(path=outputPathExcel, Temp="Temp")

#CSV#
inputFile = "_GUI/UniProt_{studying}_Human.csv".format(studying=studying)
outputHomoloGene = "{output}_{study}_Core_HomoloGene_Checked.csv".format(study=studying,output=outputPathExcel)
resultFile = "{output}_{study}_Core_UnitProtID_Compiled.csv".format(output=outputPathExcel, study=studying)
failedFile = "{output}_{study}_Core_UniProtID_Failed.csv".format(output=outputPathExcel, study=studying)


## Functions ##
def RemoveExistingFiles(file):
    for files in file:
        if os.path.exists(files):
            os.remove(files)

def WriteToCsv(path, row, fileMode):
    with open(path, fileMode) as fd:
        wr = csv.writer(fd, quoting=csv.QUOTE_ALL)
        wr.writerow(row)

def SplitListEntry(item):
    print(item)
    splitItem = item.split('_')

    print(splitItem)

    splitItemByPeriod = splitItem[2].split('.')

    print(splitItemByPeriod)
    #Declared Variables
    protName = splitItem[1]
    idNum = splitItemByPeriod[-2]

    return protName, idNum

def OpenCSV(inputFile):
    columnList = ['Entry', 'Gene', 'Organism']
    df = pd.read_csv(inputFile, usecols=columnList)

    # Reformatting Entries if not already done
    df = df[df['Entry'].notna()]
    df = df[df['Gene'].notna()]
    df['Gene'] = df['Gene'].str.upper()
    df = df[df['Organism'].notna()]
    df.drop_duplicates(subset="Gene",
                       keep="first", inplace=True)

    gene = df['Gene'].tolist()
    gene.pop(0)

    return gene

def HomoloGeneCheck(lst):
    for query in lst:
        print("Checking HomoloGene for... ", query)
        urlESearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={db}&term={" \
                     "query}&api_key=1eaf48f05dc4f266c714b46a2bd351e8fa09".format(db=db, query=query)
        loopNum = 0
        while loopNum < 3:
            try:
                # print(urlESearch)
                r = requests.get(url=urlESearch, stream=False)
                root = ET.fromstring(r.content)
                # print(r.content) #prints the xml content
                break
            except:
                loopNum +=1


        idList = []
        for child in root.iter('Id'):
            idList.append(child.text)
        # print("idList: ", idList)

        for uid in idList:
            urlESummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db={db}&id={" \
                          "uid}&version=2.0&api_key=1eaf48f05dc4f266c714b46a2bd351e8fa09".format(
                db=db, uid=uid)
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

            # else:
            #     # print(query, "is not core Nucleolar Protein.", '\n')

def UniProtRetrieve(inputHomoloGeneFile):
    # print(lst)

    #Use Pandas to Open and store Organism+Protein for each column
    columnList = ["Human", "Cow", "Mouse", "Zebrafish", "Fly", "Worm", "Yeast", "Frog"]
    df = pd.read_csv(inputHomoloGeneFile, usecols=columnList)

    #Filter our terms that contain LOC
    for name in columnList:
        df = df[~df[name].str.contains("LOC")]

    #List for Pandas
    humanProteinList = df["Human"].tolist()
    cowProteinList = df["Cow"].tolist()
    mouseProteinList = df["Mouse"].tolist()
    zebraFishProteinList = df["Zebrafish"].tolist()
    flyProteinList = df["Fly"].tolist()
    wormProteinListRevised = []
    wormProteinList = df["Worm"].tolist()
    for prot in wormProteinList:
        protRevised = prot.replace('_', '-')
        wormProteinListRevised.append(protRevised)
    yeastProteinList = df["Yeast"].tolist()
    frogProteinList = df["Frog"].tolist()



    compiledProteinList = list(zip(humanProteinList,cowProteinList,mouseProteinList, zebraFishProteinList,
                                   flyProteinList, wormProteinListRevised,yeastProteinList, frogProteinList))

    #Queries Uniprot for Excel Files, and store them as temp files
    ##########CHANGE################
    failedLst = []
    for element, lst in enumerate(compiledProteinList[0:]):
        for num, query in enumerate(lst):

            print("Querying...", query, "for", organismNameDict[num+1])
            url = "https://www.uniprot.org/uniprot/?query={query}+AND+organism:{organism}&columns=id," \
                  "genes(PREFERRED),organism&format=xls".format(query=query,
                                                                organism=organismTaxIDDict_Revised[num + 1])
            print(url)
            #Tries to write 5 times if there is an error
            loopNum = 0
            while loopNum < 5:
                try:
                    r = requests.get(url, stream=False, allow_redirects=True)
                    tempFile = "{temp}{protNum}_{query}_{organism}.xlsx".format(temp=tempPath, protNum=(element+1), query=query,
                                                                      organism=organismTaxIDDict_Revised[num + 1])
                    open(tempFile, 'wb').write(r.content)
                    print("Writing file...")
                    sizeFile = os.path.getsize(tempFile)

                    if sizeFile <=0:
                        try:
                            print("Checking File size", sizeFile)
                            os.remove(tempFile)
                            failedLst.append(query)
                            failedLst.append(organismTaxIDNameDict[num+1])
                            WriteToCsv(failedFile,failedLst,'a')
                            print("Failed: ", query, "for", organismTaxIDNameDict[num+1])
                        except:
                            print("No File to delete.", '\n')
                            break

                    else:
                        break
                except:
                    print("Failed Writing ", query, "for", organismNameDict[num+1],"; ","Loop: ", loopNum)
                    loopNum +=1
                if loopNum == 2:
                    failedLst.append(query)
                    failedLst.append(organismNameDict[num + 1])
                    WriteToCsv(failedFile, failedLst, 'a')


def UniProtIDChecker():
    def all_same(items):
        return all(x == items[0] for x in items)

    tempPathList = []
    for f in os.listdir(tempPath):
        tempPathList.append(f)

    tempPathListStr = [str(i) for i in tempPathList]
    tempPathListStr.sort()
    #print(tempPathListStr)

    #Check protein name and retrieve the UniProtID

    #print("ID          Protein      Organism")
    for items in tempPathListStr:
        #Opens the Files
        columnList = ['Entry', 'Gene names  (primary )', 'Organism']
        dfProtCheck = pd.read_excel("{temp}{item}".format(temp=tempPath,item=items), usecols=columnList)
        dfProtCheck.rename(columns={'Gene names  (primary )': 'Gene'}, inplace=True)

        dfentryList = dfProtCheck['Entry'].tolist()
        dfgeneList = dfProtCheck['Gene'].tolist()



        #Obtaining the Query Protein from the name of file
        protName, idNum = SplitListEntry(items)
        # print(idNum)
        # print("ProtName: ", protName)
        # print("Entry List: ", dfentryList)
        # print("List: ", dfgeneList)
        finalList = []
        failedList = []

        if len(dfgeneList) <= 2:
            print(dfentryList[0],"    ", protName,"       ", organismTaxIDNameDict[int(idNum)])
            print('\n')

            finalList.append(dfentryList[0])
            finalList.append(protName)
            finalList.append(organismTaxIDNameDict[int(idNum)])
            WriteToCsv(resultFile, finalList, 'a')

        elif protName in dfgeneList:
            positionNum = dfgeneList.index(protName)
            print(dfentryList[positionNum],"    ", protName,"       ", organismTaxIDNameDict[int(idNum)])
            print('\n')

            finalList.append(dfentryList[positionNum])
            finalList.append(protName)
            finalList.append(organismTaxIDNameDict[int(idNum)])
            WriteToCsv(resultFile, finalList, 'a')

        elif all_same(dfgeneList) == True:
            print(dfentryList[0], "    ", protName, "       ", organismTaxIDNameDict[int(idNum)])
            print('\n')

            finalList.append(dfentryList[0])
            finalList.append(protName)
            finalList.append(organismTaxIDNameDict[int(idNum)])
            WriteToCsv(resultFile, finalList, 'a')
            #print(protName, "appended to list")

        #Might need to add elif statement for length here, or splitting characters
        else:
            print(protName, "   NOT FOUND!", organismTaxIDNameDict[int(idNum)])
            print('\n')
            failedList.append(protName)
            failedList.append(organismTaxIDNameDict[int(idNum)])
            WriteToCsv(failedFile, failedList, 'a')
        # if len(dfGeneList) == 1 or j[0:2] == query[0:2]


##----------- Execution ----------- ##

##Remove Files if they exist:

try:
    files = [outputHomoloGene, resultFile, failedFile]
    RemoveExistingFiles(files)
except:
    print("No Existing files exist to delete.")
WriteToCsv(outputHomoloGene, ["Human", "Cow", "Mouse", "Zebrafish", "Fly", "Worm", "Yeast", "Frog"], 'a')
WriteToCsv(resultFile, ["Entry", "Gene", "Organism"], 'a')
WriteToCsv(failedFile, ["Gene", "Organism"], 'a')

try:
    geneList = OpenCSV(inputFile)
    HomoloGeneCheck(geneList)
    UniProtRetrieve(outputHomoloGene)
    UniProtIDChecker()

    shutil.rmtree(tempPath)

except:
    print("Failed")

