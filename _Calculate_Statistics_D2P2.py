import csv
import math
import os
import pandas as pd
import statistics
import numpy


def CalculateStatisticsD2P2(query, inputPath):
    # Global Variables
    columnLst = ['Position', 'Score']
    capitalizeQuery = query.capitalize()

    # Paths/Files/Input/Output #
    corePathHuman = "{path}Core/_Human/".format(path=inputPath)
    nonCorePathHuman = "{path}NonCore/_Human/".format(path=inputPath)

    corePathCow = "{path}Core/_Cow/".format(path=inputPath)
    nonCorePathCow = "{path}NonCore/_Cow/".format(path=inputPath)

    corePathMouse = "{path}Core/_Mouse/".format(path=inputPath)
    nonCorePathMouse = "{path}NonCore/_Mouse/".format(path=inputPath)

    corePathZebrafish = "{path}Core/_Zebrafish/".format(path=inputPath)
    nonCorePathZebrafish = "{path}NonCore/_Zebrafish/".format(path=inputPath)

    corePathFly = "{path}Core/_Fly/".format(path=inputPath)
    nonCorePathFly = "{path}NonCore/_Fly/".format(path=inputPath)

    corePathWorm = "{path}Core/_Worm/".format(path=inputPath)
    nonCorePathWorm = "{path}NonCore/_Worm/".format(path=inputPath)

    corePathYeast = "{path}Core/_Yeast/".format(path=inputPath)
    nonCorePathYeast = "{path}NonCore/_Yeast/".format(path=inputPath)

    corePathFrog = "{path}Core/_Frog/".format(path=inputPath)
    nonCorePathFrog = "{path}NonCore/_Frog/".format(path=inputPath)

    outputPath = "{path}{studying}_Disorder_Statistics.csv".format(path=inputPath, studying=capitalizeQuery)

    # Find all the files in each organism and open the csv
    def CalculateStats(inputPaths, orgnName, coreOrNonCore):
        # Search Files and Store as List
        fileLst = []
        for r, d, f in os.walk(inputPaths):
            for file in f:
                if file.endswith(".csv") and file.startswith("{studying}".format(studying=capitalizeQuery)):
                    fileLst.append(file)

        # Open Files, Calculate Averages, and Store as list
        avgLst = []
        avgLstFirst = []
        avgLstSecond = []
        avgLstThird = []

        for element in fileLst:
            df = pd.read_csv("{inputPath}{element}".format(inputPath=inputPaths, element=element), usecols=columnLst)
            scoreLst = df['Score'].tolist()
            avgLst.append(statistics.mean(scoreLst))

            splitScoreLst = numpy.array_split(scoreLst, 3)

            avgLstFirst.append(statistics.mean(splitScoreLst[0].tolist()))
            avgLstSecond.append(statistics.mean(splitScoreLst[1].tolist()))
            avgLstThird.append(statistics.mean(splitScoreLst[2].tolist()))

        # Calculate the Averages of the AverageLst
        avgOfAvgs = statistics.mean(avgLst)

        avgOfFirst = statistics.mean(avgLstFirst)
        avgOfSecond = statistics.mean(avgLstSecond)
        avgOfThird = statistics.mean(avgLstThird)

        # Calculate the Standard Deviation of AverageLst divided by sqrt of number of samples
        try:
            stDevOfAvgs = statistics.stdev(avgLst) / (math.sqrt(len(fileLst)))

            stdevOfFirst = statistics.stdev(avgLstFirst) / (math.sqrt(len(fileLst)))
            stdevOfSecond = statistics.stdev(avgLstSecond) / (math.sqrt(len(fileLst)))
            stdevOfThird = statistics.stdev(avgLstThird) / (math.sqrt(len(fileLst)))

        except:
            print("Cannot Calculate Variance for ", orgnName, "because variance requires at least two data points")
            stDevOfAvgs = "Cannot Calculate"
            stdevOfFirst = "Cannot Calculate"
            stdevOfSecond = "Cannot Calculate"
            stdevOfThird = "Cannot Calculate"

        # Calculate the 1/3s of the disordered regions

        # Write to CSV
        dataLst = [orgnName, coreOrNonCore, avgOfAvgs, stDevOfAvgs, avgOfFirst, stdevOfFirst, avgOfSecond,
                   stdevOfSecond, avgOfThird, stdevOfThird]

        WriteToCSV(outputPath, dataLst, 'a')

    def WriteToCSV(pathOutput, row, fileMode):
        with open(pathOutput, fileMode) as fd:
            wr = csv.writer(fd, quoting=csv.QUOTE_ALL)
            wr.writerow(row)

    ##----------- Execution ----------- ##
    try:
        os.remove(outputPath)
    except:
        pass

    WriteToCSV(outputPath,
               ['Organism', 'Core or NonCore', 'Mean of Sample Mean', 'Standard Deviation', 'Mean of Sample Mean First',
                'Std Dev First', 'Mean of Sample Mean Second', 'Std Dev Second', 'Mean of Sample Mean Third',
                'Std Dev Third'], 'a')

    CalculateStats(corePathHuman, "Human", "Core")
    CalculateStats(nonCorePathHuman, "Human", "NonCore")

    CalculateStats(corePathCow, "Cow", "Core")
    CalculateStats(nonCorePathCow, "Cow", "NonCore")

    CalculateStats(corePathMouse, "Mouse", "Core")
    CalculateStats(nonCorePathMouse, "Mouse", "NonCore")

    CalculateStats(corePathZebrafish, "Zebrafish", "Core")
    CalculateStats(nonCorePathZebrafish, "Zebrafish", "NonCore")

    CalculateStats(corePathFly, "Fly", "Core")
    CalculateStats(nonCorePathFly, "Fly", "NonCore")

    CalculateStats(corePathWorm, "Worm", "Core")
    CalculateStats(nonCorePathWorm, "Worm", "NonCore")

    CalculateStats(corePathYeast, "Yeast", "Core")
    CalculateStats(nonCorePathYeast, "Yeast", "NonCore")

    CalculateStats(corePathFrog, "Frog", "Core")
    CalculateStats(nonCorePathFrog, "Frog", "NonCore")
