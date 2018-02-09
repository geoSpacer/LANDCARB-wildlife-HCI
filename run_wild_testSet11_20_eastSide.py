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

# import wildHCI_landcarb as wildHCI
import os, glob, zipfile

region = "West Cascades"

rootDir = "R:\\LandCarb\\LandCarb.Runs\\ProductionFrontiers\\WestCascades"
# gc.set_debug(gc.DEBUG_LEAK)

batFileName = "C:\\python27\\ArcGIS10.2\\goWildtest2.bat"
batFile = open(batFileName, 'w')

# runSetList = os.listdir(rootDir)
#runSetList = ["_Test.1.9", "_Test.1.26", "_Test.1.27", "_Test.1.28", "_Test.1.29", "_Test.1.30", "_Test.1.31", "_Test.1.32", "_Test.1.33", "_Test.1.34", "_Test.1.35", "_Test.1.41", "_Test.1.42"]
#runSetList = ["_Test.1.1", "_Test.1.2", "_Test.1.3", "_Test.1.4", "_Test.1.5", "_Test.1.6", "_Test.1.19", "_Test.1.20", "_Test.1.40", "_Test.1.43", "_Test.1.44"]
runSetList = ["Test.1.2", "Test.1.3", "Test.1.4", "Test.1.5", "Test.1.8", "Test.1.23", "Test.1.25", "Test.1.26", "Test.1.30", "Test.1.33", "Test.1.34", "Test.1.37", "Test.1.39", "Test.1.43"]
#runSetList = ["Test.1.1.FullExtent", "Test.1.9.FullExtent", "Test.1.22.FullExtent"]

#runSetList = os.listdir(rootDir)

for runSetDir in runSetList:
    if os.path.isdir(rootDir + "\\" + runSetDir):
        runList = os.listdir(rootDir + "\\" + runSetDir)
        for runDir in runList:
            if (not os.path.isdir(rootDir + "\\" + runSetDir + "\\" + runDir + "\\FullExtent\\OutputFiles\\Rasters")) and os.path.isdir(rootDir + "\\" + runSetDir + "\\" + runDir + "\\EastSide\\OutputFiles\\Rasters"):
                # extract zipped files
                outPath = rootDir + "\\" + runSetDir + "\\" + runDir + "\\EastSide\\OutputFiles\\Rasters"

                # delete any exising *.img files
                for delFile in glob.glob(outPath + "\\*.img"):
                    os.remove(delFile)

                for unzipFile in ["AgeRasters.zip","HeightUpperTreeRasters.zip","LiveMassAboveGroundRasters.zip","SnagMassRasters.zip","LogMassRasters.zip"]:
                    zip = zipfile.ZipFile(outPath + "\\" + unzipFile)
                    zip.extractall(path=outPath)

                timeSet = glob.glob(rootDir + "\\" + runSetDir + "\\" + runDir + "\\EastSide\\OutputFiles\\Rasters\\LiveMassAboveGround_*")
                if len(timeSet) == 0:
                    print "*********Error!, output files not found for: " + runSetDir + " run " + runDir
                else:
                    timeSetSub = []
                    for timeVal in range(0,len(timeSet)):
                        timeSet[timeVal] = int(timeSet[timeVal][-8:-4])
    #                    if (timeSet[timeVal] >= 2210):
                        timeSetSub.append(timeSet[timeVal])

                    print runSetDir + " run " + runDir + " " + str(timeSetSub)

                for species in ["bluebird", "treevole", "olivesided", "muledeer"]:
    #                for species in ["bluebird"]:
                    try:
    #                    wildHCI.main(region, species, runSetDir, runDir, timeSetSub[-10:])
                        batFile.write("C:\\python27\\ArcGIS10.2\\python.exe N:\\Code\\Python\\NASA_Carbon\\wildlife\\wildHCI_landcarb.py \"" + region + "\" " + species + " " + runSetDir + " " + runDir + " 0 \"" + str(timeSetSub) + "\"\n")

                    except Exception, e:
                        print "\n\n" + species + ": " + e.args[0]

                    except:
                        print "unhandled Error!!"

batFile.close()
print 'done.'

