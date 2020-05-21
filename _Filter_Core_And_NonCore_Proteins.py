import os
import urllib.request
import pandas as pd


def FilterCoreNonCoreProt(queryName, inputPath):
    # Global Variables
    capitalizeQuery = queryName.capitalize()

    # Path/Files/Input/Output
    inputFile = "{path}_HomoloGeneSets/_{studying}_Core_HomoloGene_List.csv".format(path=inputPath,
                                                                                    studying=capitalizeQuery)
    inputFullSet = "{path}_FullProteinSets/{studying}_UniProt_Full.csv".format(path=inputPath, studying=capitalizeQuery)
    outputPathCore = "{path}Core/{studying}_Core_Proteins.csv".format(path=inputPath, studying=capitalizeQuery)
    outputPathNonCore = "{path}NonCore/{studying}_NonCore_Proteins.csv".format(path=inputPath, studying=capitalizeQuery)

    def OpenCSV(lst):

        compiledLst = []

        col_list = ['Human', 'Cow', 'Mouse', 'Zebrafish', 'Fly', 'Worm', 'Yeast', 'Frog']
        df = pd.read_csv(lst, usecols=col_list)

        humanLst = df['Human'].tolist()
        humanLst.insert(0, "Human")
        compiledLst.append(humanLst)

        cowLst = df['Cow'].tolist()
        cowLst.insert(0, "Cow")
        compiledLst.append(cowLst)

        mouseLst = df['Mouse'].tolist()
        mouseLst.insert(0, "Mouse")
        compiledLst.append(mouseLst)

        zebrafishLst = df['Zebrafish'].tolist()
        zebrafishLst.insert(0, "Zebrafish")
        compiledLst.append(zebrafishLst)

        flyLst = df['Fly'].tolist()
        flyLst.insert(0, "Fly")
        compiledLst.append(flyLst)

        wormLst = df['Worm'].tolist()
        wormLst.insert(0, "Worm")
        compiledLst.append(wormLst)

        yeastLst = df['Yeast'].tolist()
        yeastLst.insert(0, "Yeast")
        compiledLst.append(yeastLst)

        frogLst = df['Frog'].tolist()
        frogLst.insert(0, "Frog")
        compiledLst.append(frogLst)

        return compiledLst

    def CheckUniProtName(inputLst):
        nameOrgn = inputLst[0]
        inputLst.pop(0)

        def ReduceListSize(reduceLst, num):
            for i in range(0, len(reduceLst), num):
                yield reduceLst[i:i + num]

        # UniProt has limit of 1000 entries..?
        reducedEntryArray = list(ReduceListSize(inputLst, 900))

        for lst in reducedEntryArray:
            newLineLst = ' '.join(lst)
            url = 'https://www.uniprot.org/uploadlists/'
            params = {
                'from': 'GENENAME',
                'to': 'ACC',
                'format': 'tab',
                'query': '{query}'.format(query=newLineLst)
            }

            data = urllib.parse.urlencode(params)
            data = data.encode('utf-8')
            req = urllib.request.Request(url, data)
            with urllib.request.urlopen(req) as f:
                response = f.read()
            words = response.decode('utf-8')
            # print(words)
            words_split = words.split()
            words_split.pop(0)  # add two to remove the 'From' and 'To'
            words_split.pop(0)

            uniProtIDs = words_split[1::2]

            # Pandas Dataframe to get Entry and Gene names #
            columnList = ['Entry', 'Gene', 'Organism']
            df = pd.read_csv(inputFullSet, usecols=columnList)
            df1 = df.loc[df['Organism'] == nameOrgn].sort_values("Gene")

            geneLst = df1['Gene'].tolist()
            entryLst = df1['Entry'].tolist()
            orgnLst = df1['Organism'].tolist()


            # Get Core Proteins
            indexLst = [i for i, item in enumerate(entryLst) if item in uniProtIDs]

            coreEntryLst = [entryLst[val] for val in indexLst]
            coreGeneLst = [geneLst[val] for val in indexLst]
            coreOrgnLst = [orgnLst[val] for val in indexLst]

            WriteToCSV(outputPathCore, coreEntryLst, coreGeneLst, coreOrgnLst)

            # Get NonCore
            nonCoreIndexLst = [i for i, item in enumerate(entryLst) if item not in uniProtIDs]

            nonCoreEntryLst = [entryLst[value] for value in nonCoreIndexLst]
            nonCoreGeneLst = [geneLst[value] for value in nonCoreIndexLst]
            nonCoreOrgnLst = [orgnLst[value] for value in nonCoreIndexLst]

            WriteToCSV(outputPathNonCore, nonCoreEntryLst, nonCoreGeneLst, nonCoreOrgnLst)

    def WriteToCSV(pathOutput, entry, gene, orgn):

        columnList = ['Entry', 'Gene', 'Organism']
        df = pd.DataFrame(list(zip(entry, gene, orgn)), columns=columnList)

        if not os.path.isfile(pathOutput):
            df.to_csv(pathOutput, header=columnList, index=False)
        else:  # else it exists so append without writing the header
            df.to_csv(pathOutput, mode='a', header=False, index=False)

    ##----------- Execution ----------- ##

    # Remove Files if they exist
    try:
        os.remove(outputPathCore)
    except:
        pass

    try:
        os.remove(outputPathNonCore)
    except:
        pass

    compiledList = OpenCSV(inputFile)

    for inputLists in compiledList:
        CheckUniProtName(inputLists)

