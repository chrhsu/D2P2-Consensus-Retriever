import sys

sys.path.append('D2P2Project/')

from D2P2Project._Retrieve_Full_Protein_Set_From_UniProt_D2P2_Check import *
from D2P2Project._Retrieve_Core_Proteins_From_HomoloGene import *
from D2P2Project._Filter_Core_And_NonCore_Proteins import *
from D2P2Project._Retrieve_D2P2_Consensus_Data import *
from D2P2Project._Calculate_Statistics_D2P2 import *

# Path is Universal for all scripts
path = "Data/"


# Functions:
def CreatePaths():
    if not os.path.exists('{path}'.format(path=path)):
        os.makedirs('{path}'.format(path=path))

    if not os.path.exists('{path}{dir}/'.format(path=path, dir="_FullProteinSets")):
        os.makedirs('{path}{dir}/'.format(path=path, dir="_FullProteinSets"))

    if not os.path.exists('{path}{dir}/'.format(path=path, dir="_HomoloGeneSets")):
        os.makedirs('{path}{dir}/'.format(path=path, dir="_HomoloGeneSets"))

    if not os.path.exists('{path}{dir}/'.format(path=path, dir="Core")):
        os.makedirs('{path}{dir}/'.format(path=path, dir="Core"))

    if not os.path.exists('{path}{dir}/'.format(path=path, dir="NonCore")):
        os.makedirs('{path}{dir}/'.format(path=path, dir="NonCore"))


def Query():
    print("What organelle do you want to query?")
    query = input("-> ")
    print("Please enter API Key From NCBI: ")
    print("API Key can be found by following directions from..", '\n', "https://ncbiinsights.ncbi.nlm.nih.gov/2017/11"
                                                                      "/02/new-api-keys-for-the-e-utilities/")
    api = input("-> ")
    api_Key = "&api_key={api}".format(api=api)
    print('\n')
    query.lower()
    return query, api_Key


def UniProt():
    RetrieveFullProtSetUniProtD2P2Check(inputQuery, path)


def HomoloGene():
    RetrieveCoreProteinsHomoloGene(inputQuery, path, apiKey)


def ConvertCoreHGToUniProt():
    print("Coverting HomoloGene to UniProt ID...", end="")
    FilterCoreNonCoreProt(inputQuery, path)
    print("Done.")


def D2P2():
    UniProtD2P2Data(inputQuery, path)


def Statistics():
    print('\n')
    print("Calculating Statistics...", end="")
    CalculateStatisticsD2P2(inputQuery, path)
    print("Done.", end="")
    print('\n')
    print("File is located here: {path}{query}_Statistics_Disorder_Statistics.csv".format(path=path,
                                                                                          query=inputQuery.capitalize()))


##----------- Execution ----------- ##

# Creates paths needed for program #
CreatePaths()

# Input from User
inputQuery, apiKey = Query()

# Retrieves the Full Dataset for query
UniProt()

# Retrieves NCBI HomoloGene Core Proteins
HomoloGene()

# Get Core and NonCore Proteins
ConvertCoreHGToUniProt()

# Get D2P2 Raw Data
D2P2()

# Statistics for Core and NonCore Proteins
Statistics()
