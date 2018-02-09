#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      olsenk
#
# Created:     21/11/2013
# Copyright:   (c) olsenk 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

def main():
    pass

if __name__ == '__main__':
    main()

import arcpy
import wildHCI_landcarb as wildHCI
import os
from math import pow, sqrt

region = "West Cascades"
#speciesList = ["bluebird","marten","olivesided","pileated","spottedowl","treevole","muledeer"]
speciesList = ["bluebird","olivesided","treevole","muledeer"]
#speciesList = ["treevole"]

for species in speciesList:

    if (species == 'bluebird'):
        sppAbb = 'wbb'
    elif (species == 'marten'):
        sppAbb = 'mrt'
    elif (species == 'olivesided'):
        sppAbb = 'osf'
    elif (species == 'pileated'):
        sppAbb = 'pwo'
    elif (species == 'spottedowl'):
        sppAbb = 'nso'
    elif (species == 'treevole'):
        sppAbb = 'rtv'
    elif (species == 'muledeer'):
        sppAbb = 'mdr'
    else:
        sppAbb = 'err'

    rootDir = "R:\\Biodiversity\\WestCascades\\Frontiers_fullextent\\" + species
    outFile = open(rootDir + "\\WC_Frontiers_" + species + "_v3.csv", 'w')
    outFile.write("treatment,pool,n,min,max,mean,std,stderr\n")

    runSetList = os.listdir(rootDir)
    for runSetDir in runSetList:
        if "Test" in runSetDir:
            dataList = []
            runList = os.listdir(rootDir + "\\" + runSetDir)
            for timeStep in runList[-10:]:
                print runSetDir + " time step " + timeStep
                rasterDir = rootDir + "\\" + runSetDir + "\\" + timeStep

                try:
                    rows = arcpy.SearchCursor(rasterDir + '\\' + sppAbb + '_wc_hci_3c')
                    rowDataList = [0,0,0]
                    for row in rows:
                        rowDataList[row.value - 1] = row.count

                    totalPixels = sum(rowDataList)
                    for rowNum in xrange(len(rowDataList)):
                        rowDataList[rowNum] = float(rowDataList[rowNum]) / totalPixels * 100

                    del rows, row

                    # skip the low habitat score and create a >33 and >66 class
                    dataList.append([rowDataList[1] + rowDataList[2], rowDataList[2]])

                except:
                    print "Time step " + timeStep + " not found"

            if (len(dataList) > 1):
                meanDataList = [0,0]
                minDataList = [100,100]
                maxDataList = [0,0]
                for yearList in dataList:
                    meanDataList[0] += yearList[0]
                    meanDataList[1] += yearList[1]

                    if (yearList[0] < minDataList[0]):
                        minDataList[0] = yearList[0]
                    if (yearList[1] < minDataList[1]):
                        minDataList[1] = yearList[1]
                    if (yearList[0] > maxDataList[0]):
                        maxDataList[0] = yearList[0]
                    if (yearList[1] > maxDataList[1]):
                        maxDataList[1] = yearList[1]

                meanDataList[0] /= float(len(dataList))
                meanDataList[1] /= float(len(dataList))

                stdDevList = [0,0]
                for level in [0,1]:
                    for yearList in dataList:
                        stdDevList[level] += pow(yearList[level] - meanDataList[level], 2)

                    stdDevList[level] = sqrt(stdDevList[level] / len(dataList))

                # write data to outFile
                runSetName = runSetDir[-8:]
                if (runSetName[0] == "."):
                    runSetName = runSetName[1:]
                elif (runSetName[0] == "s"):
                    runSetName = runSetName[3:]

                outFile.write(runSetName + "," + sppAbb + "_HCI_gt33," + str(len(dataList)) + ",")
                outFile.write(str(minDataList[0]) + "," + str(maxDataList[0]) + "," + str(meanDataList[0]))
                outFile.write("," + str(stdDevList[0]) + "," + str(stdDevList[0] / sqrt(len(dataList))) + "\n")

                outFile.write(runSetName + "," + sppAbb + "_HCI_gt66," + str(len(dataList)) + ",")
                outFile.write(str(minDataList[1]) + "," + str(maxDataList[1]) + "," + str(meanDataList[1]))
                outFile.write("," + str(stdDevList[1]) + "," + str(stdDevList[1] / sqrt(len(dataList))) + "\n")

    outFile.close()
    print 'done.'

