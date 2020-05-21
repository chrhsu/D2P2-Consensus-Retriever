import os

import pandas as pd
import json
from urllib.request import urlopen


### Purpose of this Script is to Query and get D2P2 Disorder Consensus Data
# Either store as CSV or manipulate the data and then store as CSV
# Store the data in this fashion: ["Position', "Score"]
# Store as name: "D2P2_[ProteinName]_[OrganismName].csv"


def UniProtD2P2Data(query, inputPath):
    # Global Variables
    capitalizeQuery = query.capitalize()
    col_list = ['Entry', 'Gene', 'Organism']

    # Paths/Files/Output
    corePath = "{path}Core/".format(path=inputPath)
    nonCorePath = "{path}NonCore/".format(path=inputPath)
    inputCoreFile = "{path}Core/{studying}_Core_Proteins.csv".format(path=inputPath, studying=capitalizeQuery)
    inputNonCoreFile = "{path}NonCore/{studying}_NonCore_Proteins.csv".format(path=inputPath, studying=capitalizeQuery)

    # Functions
    def OpenCSV(FileInput):
        df = pd.read_csv(FileInput, usecols=col_list)
        entryList = df['Entry'].tolist()
        geneList = df['Gene'].tolist()
        organismList = df['Organism'].tolist()

        return entryList, geneList, organismList

    def D2P2RawData(entryList, geneList, orgnList, coreOrNonCore):
        print('\n')
        def GetPositionRange(lst):
            positionList = []
            lengthList = len(lst)
            for index in range(1, lengthList + 1):
                positionList.append(index)
            return positionList

        for i, seqID in enumerate(entryList):

            print("Querying D2P2", coreOrNonCore, ": ", geneList[i], "for", orgnList[i], "... ", end="")

            requestURL = 'http://d2p2.pro/api/seqid/["{seqid}"]'.format(seqid=seqID)
            response = json.loads(urlopen(requestURL).read())

            consensusList = []
            for consensus in response["{seqid}".format(seqid=seqID)][0][2]["disorder"]["consensus"]:
                consensusList.append(consensus)

            positionList = GetPositionRange(consensusList)
            WriteToCSV(positionList, consensusList, geneList[i], orgnList[i], coreOrNonCore)

            print("Done.")

    def WriteToCSV(posLst, d2p2Lst, geneName, orgnName, nonCoreOrCore):
        if not os.path.exists('{pathCore}_{orgnName}/'.format(pathCore=corePath, orgnName=orgnName)):
            os.makedirs('{pathCore}_{orgnName}/'.format(pathCore=corePath, orgnName=orgnName))
        if not os.path.exists('{pathNonCore}_{orgnName}/'.format(pathNonCore=nonCorePath, orgnName=orgnName)):
            os.makedirs('{pathNonCore}_{orgnName}/'.format(pathNonCore=nonCorePath, orgnName=orgnName))

        outputPathCore = "{pathCore}_{orgnName}/".format(pathCore=corePath, orgnName=orgnName)
        outputPathNonCore = "{pathNonCore}_{orgnName}/".format(pathNonCore=nonCorePath, orgnName=orgnName)

        df = pd.DataFrame(list(zip(posLst, d2p2Lst)), columns=["Position", "Score"])

        if nonCoreOrCore == "Core":
            df.to_csv("{outputPath}{studying}_{geneName}_D2P2_Consensus.csv".format(studying=capitalizeQuery, outputPath=outputPathCore, geneName=geneName),
                      header=["Position", "Score"], index=False)
        if nonCoreOrCore == "NonCore":
            df.to_csv("{outputPath}{studying}_{geneName}_D2P2_Consensus.csv".format(studying=capitalizeQuery, outputPath=outputPathNonCore, geneName=geneName),
                      header=["Position", "Score"], index=False)

    ###----Execution---###

    # Opening the CSV files
    coreEntryLst, coreGeneLst, coreOrgnLst = OpenCSV(inputCoreFile)
    nonCoreEntryLst, nonCoreGeneLst, nonCoreOrgnLst = OpenCSV(inputNonCoreFile)

    # Retrieving D2P2 Data
    D2P2RawData(coreEntryLst, coreGeneLst, coreOrgnLst, "Core")
    D2P2RawData(nonCoreEntryLst, nonCoreGeneLst, nonCoreOrgnLst, "NonCore")
