#-------------------------------------------------------------------------------
# Name:        wildHCI_landcarb.py
# Purpose:     This script will predict stand structure variables from LandCarb output
#               and run wildlife models on them.
# Author:      Keith Olsen
#
# Created:     18 December 2012
# Copyright:   (c) olsenk 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

# import libraries to support model
import arcpy
import numpy as np
# from math import log, sqrt, pow
import rasterFunctions as rastFunc
import datetime, os, sys, shutil

def errorHandler(logFile, indexAreaFile, indexAreaFile2, message):
    logFile.write("\n\n" + message)
    logFile.close()
    indexAreaFile.close()
    indexAreaFile2.close()
    raise Exception(message)

def main(region, species, runSetNum, runNum, repNum, yearList):
    # choose region (John Day, West Cascades, Tillamook)
    # choose species (strVars, bluebird, pileated)
    # choose source (gnn, landcarb)
    source = 'landcarb'
    # yearList = [2920, 2930, 2940, 2950, 2960, 2970, 2980, 2990, 3000, 3010]
    yearList = yearList[1:-1].split(", ")

    if (region == 'John Day'):
        regionDir = 'JohnDay'
        regionAbb = 'jd'
    elif (region == 'West Cascades'):
        regionDir = 'WestCascades'
        regionAbb = 'wc'
    elif (region == 'Tillamook'):
        regionDir = 'Tillamook'
        regionAbb = 'ti'

    # input directory for predictior variables
    if "Test" in runSetNum:
        landCarbDir = "R:\\LandCarb\\LandCarb.Runs\\ProductionFrontiers\\" + regionDir + "\\" + runSetNum + "\\" + runNum + "\\EastSide\\OutputFiles\\Rasters"
    else:
        if repNum == "0":
            landCarbDir = "R:\\LandCarb\\LandCarb.Runs\\Trajectories\\" + regionDir + "\\RunSet." + runSetNum + "\\Run." + runNum + "\\FullExtent\\OutputFiles\\Rasters"
        else:
            landCarbDir = "R:\\LandCarb\\LandCarb.Runs\\Trajectories\\" + regionDir + "\\RunSet." + runSetNum + "\\Run." + runNum + "\\Rep." + repNum + "\\FullExtent\\OutputFiles\\Rasters"

    gnnDir = "R:\\Biodiversity\\gnn06_00_stats\\"

    # input directory for static predictor variables
    # staticPredictorDir = r'R:\Biodiversity\static_predictors'

    if (source == 'landcarb'):
        # output directory root
        if "Test" in runSetNum:
            outputDirRoot = "C:\\temp\\Biodiversity\\" + regionDir + "\\Frontiers\\" + species + "\\" + runSetNum + "_" + runNum
        else:
            outputDirRoot = "C:\\temp\\Biodiversity\\" + regionDir + "\\" + species + "\\LandCarb_s" + runSetNum + "_run" + runNum + "_rep" + repNum

    elif (source == 'gnn'):
        # output directory root
        outputDirRoot = "C:\\temp\\Biodiversity\\" + regionDir + "\\" + species + "\\gnn_06_00"

    if not os.path.exists(outputDirRoot):
        os.makedirs(outputDirRoot)

    indexAreaFile = open(outputDirRoot + "\\HCI33areaByYear.csv", 'w')
    indexAreaFile2 = open(outputDirRoot + "\\HCI66areaByYear.csv", 'w')
    logFile = open(outputDirRoot + "\\runlog.txt", 'w')

    if (species == 'bluebird'):
        varList = ['STPH_25_50', 'STPH_GE_50', 'CANCOV']
        # rastFunc.writeAreaFileHeader('Year, wbb_ls, wbb_ss, wbb_cc, wbb_nci, wbb_hci', indexAreaFile, indexAreaFile2)
    elif (species == 'pileated'):
        varList = ['STPH_50_75', 'STPH_GE_75', 'STPH_GE_25', 'TPH_GE_50', 'DVPH_GE_25', 'DVPH_GE_50']
        # rastFunc.writeAreaFileHeader('Year, pwo_ms, pwo_ls, pwo_nci, pwo_ts, pwo_dls, pwo_fi, pwo_fci, pwo_hci', indexAreaFile, indexAreaFile2)
    elif (species == 'treevole'):
        varList = ['TPHC_GE_50', 'TPHDF_GE_50', 'TPH_GE_50', 'QMDA_GE_3', 'CANCOV', 'DDI']
        # rastFunc.writeAreaFileHeader('Year, rtv_con, rtv_qmd, rtv_cc, rtv_ddi, rtv_hci', indexAreaFile, indexAreaFile2)
    elif (species == 'blackbacked'):
        varList = ['STPH_GE_30', 'STPH_10_30', 'CANCOV']
        # rastFunc.writeAreaFileHeader('Year, bbw_ls, bbw_ss, bbw_cc, bbw_hci', indexAreaFile, indexAreaFile2)
    elif (species == 'whiteheaded'):
        varList = ['STPH_GE_25', 'CANCOV', 'TPH_GE_25']
        # rastFunc.writeAreaFileHeader('Year, whw_ls, whw_ncc, whw_nci, whw_ts, whw_fcc, whw_fi, whw_fci, whw_hci', indexAreaFile, indexAreaFile2)
    elif (species == 'olivesided'):
        varList = ['CANCOV', 'TPH_GE_10', 'STPH_GE_10', 'STNDHGT']
        # rastFunc.writeAreaFileHeader('Year, osf_cc, osf_st, osf_lt, osf_nci, osf_fci, osf_hci', indexAreaFile, indexAreaFile2)
    elif (species == 'spottedowl'):
        varList = ['TPH_10_25', 'TPH_25_50', 'TPH_GE_75', 'DDI', 'AGE_DOM_BA']
        # rastFunc.writeAreaFileHeader('Year, nso_smt, nso_mdt, nso_lgt, nso_ddi, nso_nci, nso_s1, nso_s2, nso_s3, nso_fci, nso_hci', indexAreaFile, indexAreaFile2)
    elif (species == 'marten'):
        varList = ['STPH_50_75', 'DVPH_GE_50', 'DLPH_GE_30', 'DDI']
        # rastFunc.writeAreaFileHeader('Year, mrt_ls, mrt_dls, mrt_nci, mrt_ddi, mrt_fci, mrt_hci', indexAreaFile, indexAreaFile2)
    elif (species == 'murrelet'):
        varList = ['AGE_DOM_BA', 'CANCOV', 'STNDHGT', 'QMDA_GE_3', 'IMAP_LAYERS']
        # rastFunc.writeAreaFileHeader('Year, mur_age, mur_cc, mur_hgt, mur_qmd, mur_lyr, mur_hci', indexAreaFile, indexAreaFile2)
    elif (species == 'goshawk'):
        varList = ['TPH_GE_3', 'QMDA_GE_10', 'STNDHGT', 'CANCOV']
        # rastFunc.writeAreaFileHeader('Year, gos_tph, gos_qmd, gos_hgt, gos_cc, gos_hci', indexAreaFile, indexAreaFile2)
    elif (species == 'muledeer'):
        if (region == 'John Day'):
            varList = ['AVGDBH_GE_3', 'CANCOV', 'IMAP_LAYERS', 'PUTR_COV']
            # rastFunc.writeAreaFileHeader('Year, mdr_fsz, mdr_fcc, mdr_fly, mdr_fbb, mdr_fh, mdr_fhs, mdr_hsz, mdr_hcc, mdr_hly, mdr_hbb, mdr_hh, mdr_hhs, mdr_tsz, mdr_tcc, mdr_tly, mdr_tc, mdr_tcs, mdr_hci', indexAreaFile, indexAreaFile2)
        else:
            varList = ['AVGDBH_GE_3', 'CANCOV', 'IMAP_LAYERS']
            # rastFunc.writeAreaFileHeader('Year, mdr_fsz, mdr_fcc, mdr_fly, mdr_fh, mdr_fhs, mdr_hsz, mdr_hcc, mdr_hly, mdr_hh, mdr_hhs, mdr_tsz, mdr_tcc, mdr_tly, mdr_tc, mdr_tcs, mdr_hci', indexAreaFile, indexAreaFile2)

    elif (species == 'strVars'):
        varList = ['LC_TBIO_MGHA', 'LC_SBIO_MGHA', 'LC_LBIO_MGHA']
##        varList = ['TPH_GE_3', 'TPH_GE_10', 'TPH_GE_25', 'TPH_GE_50', 'TPH_GE_75', 'TPH_10_25', 'TPH_25_50', 'TPHC_GE_50', 'TPHDF_GE_50', 'STPH_GE_10', 'STPH_GE_25',
##        'STPH_GE_30', 'STPH_GE_50', 'STPH_GE_75', 'STPH_10_30', 'STPH_25_50', 'STPH_50_75', 'DLPH_GE_30', 'DVPH_GE_25', 'DVPH_GE_50', 'CANCOV', 'QMDA_GE_3', 'QMDA_GE_10',
##        'DDI', 'AVGDBH_GE_3', 'IMAP_LAYERS', 'STNDHGT', 'AGE_DOM_BA']
        # varList = ['PUTR_COV']

    printAreaHeader = True
    headerList = []
    averageValList = []
    for lcYear in yearList:
        # output directory for model output
        outputDir = outputDirRoot + "\\yr" + lcYear

        # check if GNN predictor run and if so reset output directory
        if (source == 'gnn'):
            outputDir = outputDirRoot

    ##    if not arcpy.Exists(outputDir):
    ##        arcpy.CreateFileGDB_management(os.path.dirname(outputDir), os.path.basename(outputDir))

    ##    if os.path.exists(outputDir):
    ##        logFile.close()
    ##        indexAreaFile.close()
    ##        shutil.rmtree(outputDir)

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        # setup log File
        logFile.write("\n\nStart time: " + datetime.date.today().strftime("%B %d, %Y") + " - " + datetime.datetime.now().time().strftime("%H:%M:%S"))
        logFile.write("\nRegion: " + region)
        logFile.write("\nSpecies: " + species)
        logFile.write("\nLandCarb Year: " + lcYear)
        print "\n\nRunning " + region + " " + species + " year " + lcYear

        area33List = []
        area66List = []
        # indexAreaFile.write("\n" + lcYear)
        # indexAreaFile2.write("\n" + lcYear)

        # load predictor variables
        # describe input reference raster
        # arcpy.env.workspace = rootDir
        if "Test" in runSetNum:
            referenceDesc = arcpy.Describe(landCarbDir + r'\LiveMassAboveGround_3010.img')
            indexStack = [rastFunc.loadMaskedRaster(landCarbDir + r'\LiveMassAboveGround_3010.img', referenceDesc, logFile)]
        elif "Trajectories" in landCarbDir:
            referenceDesc = arcpy.Describe(landCarbDir + r'\LiveMassAboveGround_' + yearList[0] + '.img')
            indexStack = [rastFunc.loadMaskedRaster(landCarbDir + r'\LiveMassAboveGround_' + yearList[0] + '.img', referenceDesc, logFile)]
        else:
            referenceDesc = arcpy.Describe(landCarbDir + r'\LiveMass_2000.img')
            indexStack = [rastFunc.loadMaskedRaster(landCarbDir + r'\LiveMass_2000.img', referenceDesc, logFile)]


        preDict = {}
        strVarDict = {}
        if (source == 'gnn'):
            if (species == 'strVars'):
                # run models with GNN predictor variables and regression models
                rastFunc.loadGNNVars(gnnDir, ['TBIO_JENK_', 'SBIO_MGHA', 'LBIO_MGHA', 'AGE_DOM_BA', 'STNDHGT'], referenceDesc, logFile, strVarDict, indexStack[0].mask)
                preDict['pRas_tBio_mgha'] = strVarDict['TBIO_JENK_']
                preDict['pRas_sBio_mgha'] = strVarDict['SBIO_MGHA']
                preDict['pRas_lBio_mgha'] = strVarDict['LBIO_MGHA']
                preDict['pRas_ageDom'] = strVarDict['AGE_DOM_BA']
                preDict['pRas_stndHgt'] = strVarDict['STNDHGT']
                strVarDict = {}
                rastFunc.calcStructureVars(varList, strVarDict, preDict, regionDir, outputDir, referenceDesc, logFile)
            else:
                # run models with GNN variables and skip regression models
                rastFunc.loadGNNVars(gnnDir, varList, referenceDesc, logFile, strVarDict, indexStack[0].mask)
        else:
            # run models with LandCarb predictor variables and regression models
            rastFunc.loadPredictors(landCarbDir, lcYear, referenceDesc, logFile, preDict)
            rastFunc.calcStructureVars(varList, strVarDict, preDict, regionDir, outputDir, referenceDesc, logFile)

        if (species == 'strVars'):
            strVarListAvg = []
            for strVarName in strVarDict:
#                rastFunc.outputMaskedRaster(outputDir + "\\" + strVarName, strVarDict[strVarName], referenceDesc, logFile, 'float')
                averageVal = np.ma.average(strVarDict[strVarName])
                strVarListAvg.append((strVarName, averageVal))

            averageValList.append(strVarListAvg)

        elif (species == 'bluebird'):
            # Build Large Snag (LS) score for the Western Bluebird
            print "*** Processing Large Snag Score"
            if (region == 'John Day'):
                nodeList = [(0,0), (3,1)]
            elif (region == 'West Cascades' or region == 'Tillamook'):
                nodeList = [(0,0), (5,1)]

            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STPH_GE_50'], 'STPH_GE_50', nodeList, logFile))
            rastFunc.outputMaskedRaster(outputDir + "\\wbb_" + regionAbb + "_LS", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ls', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'value', 0.6, 'multiply', logFile)

            # Build Small Snag (SS) score for the Western Bluebird
            print "*** Processing Small Snag Score"
            if (region == 'John Day'):
                nodeList = [(0,0), (3,1)]
            elif (region == 'West Cascades' or region == 'Tillamook'):
                nodeList = [(0,0), (15,1)]

            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STPH_25_50'], 'STPH_25_50', nodeList, logFile))
            rastFunc.outputMaskedRaster(outputDir + "\\wbb_" + regionAbb + "_SS", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ss', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'value', 0.4, 'multiply', logFile)
            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build canopy index (CC) for the Western Bluebird
            print "*** Processing Canopy Index"
            nodeList = [(0,1), (20,1), (45,0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))
            rastFunc.outputMaskedRaster(outputDir + "\\wbb_" + regionAbb + "_CC", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'cc', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'multiply', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.5, 'power', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\wbb_" + regionAbb + "_NCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'nci', headerList, area33List, area66List, printAreaHeader)

            print "*** Processing elevation restriction"
            elevationImg = rastFunc.loadMaskedRaster(r'T:\groups\CarbonTradeoffs\phy\dem\or_dem100', referenceDesc, logFile)
            elevationImg = np.ma.masked_where(indexStack[len(indexStack) - 1].mask, elevationImg)

            if (region == 'West Cascades' or region == 'Tillamook'):
                indexStack[len(indexStack) - 1] = np.ma.where(elevationImg > (4000 * 0.3048), 0.0, indexStack[len(indexStack) - 1])
            elif (region == 'John Day' or region == 'East Cascades'):
                indexStack[len(indexStack) - 1] = np.ma.where(elevationImg > (5800 * 0.3048), 0.0, indexStack[len(indexStack) - 1])

            rastFunc.outputMaskedRaster(outputDir + "\\wbb_" + regionAbb + "_HCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\wbb_" + regionAbb + "_HCI_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        elif (species == 'pileated'):
            # Build medium Snag (MS) score for the Pileated Woodpecker
            print "*** Processing Medium Snag Score"
            if (region == 'John Day'):
                nodeList = [(0,0), (20,1)]
                indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STPH_GE_25'], 'STPH_GE_25', nodeList, logFile))
                indexAreaFile.write(',0,0')
            else:
                nodeList = [(0,0), (3,1)]
                indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STPH_50_75'], 'STPH_50_75', nodeList, logFile))

                rastFunc.outputMaskedRaster(outputDir + "\\pwo_" + regionAbb + "_MS", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
                rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ms', headerList, area33List, area66List, printAreaHeader)

                rastFunc.mathFunc(indexStack, 'value', 0.4, 'multiply', logFile)

                # Build large Snag (LS) score for the Pileated Woodpecker
                print "*** Processing Large Snag Score"
                nodeList = [(0,0), (2,1)]
                indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STPH_GE_75'], 'STPH_GE_75', nodeList, logFile))

                rastFunc.outputMaskedRaster(outputDir + "\\pwo_" + regionAbb + "_LS", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
                rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ls', headerList, area33List, area66List, printAreaHeader)

                rastFunc.mathFunc(indexStack, 'value', 0.6, 'multiply', logFile)
                rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\pwo_" + regionAbb + "_NCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'nci', headerList, area33List, area66List, printAreaHeader)


            rastFunc.mathFunc(indexStack, 'image', 0, 'copy', logFile)

            # Build live tree score (TS) for the Pileated Woodpecker
            print "*** Processing Live Tree Score"
            if (region == 'John Day'):
                nodeList = [(0,0), (75,1)]
            elif (region == 'West Cascades' or region == 'Tillamook'):
                nodeList = [(0,0), (180,1)]

            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['TPH_GE_50'], 'TPH_GE_50', nodeList, logFile))
            rastFunc.outputMaskedRaster(outputDir + "\\pwo_" + regionAbb + "_TS", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ts', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build down log score (DLS) for the Pileated Woodpecker
            print "*** Processing Down Log Score"
            if (region == 'John Day'):
                nodeList = [(0,0), (170,1)]
                indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['DVPH_GE_25'], 'DVPH_GE_25', nodeList, logFile))
            elif (region == 'West Cascades' or region == 'Tillamook'):
                nodeList = [(0,0), (23,1)]
                indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['DVPH_GE_50'], 'DVPH_GE_50', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\pwo_" + regionAbb + "_DLS", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'dls', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.3333, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\pwo_" + regionAbb + "_FI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fi', headerList, area33List, area66List, printAreaHeader)

            # Get mean score for circular neighborhood surrounding the focal pixel
            print "*** Processing Neighborhood Average for Landscape Score"
            if (region == 'John Day'):
                rastFunc.circularWindowStats(indexStack, 130, referenceDesc, 'mean', logFile)
            elif (region == 'West Cascades' or region == 'Tillamook'):
                rastFunc.circularWindowStats(indexStack, 400, referenceDesc, 'mean', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\pwo_" + regionAbb + "_FCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fci', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'minimum', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\pwo_" + regionAbb + "_HCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\pwo_" + regionAbb + "_HCI_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        elif (species == 'treevole'):
            if (region == 'John Day'):
                errorHandler(logFile, indexAreaFile, indexAreaFile2, "**Warning. No model variation for John Day")
            elif (region == 'West Cascades'):
                # Build douglas-fir composition index for the Red Tree Vole
                print "*** Processing Douglas-fir index"
                indexStack.append(np.ma.where(strVarDict['TPH_GE_50'] > 0.0, strVarDict['TPHDF_GE_50'] / strVarDict['TPH_GE_50'], 0.0))
                # truncate proportion of douglas-fir to 1. Regression equations allow for this to be >1
                indexStack[len(indexStack) - 1] = np.ma.where(indexStack[len(indexStack) - 1] < 1.0, indexStack[len(indexStack) - 1], 1.0)
                logFile.write("\n\nAdding image and processing DFTPH_GE_50 / TPH_GE_50")
            elif (region == 'Tillamook'):
                # Build conifer composition index for the Red Tree Vole
                print "*** Processing Conifer Composition Index"
                indexStack.append(np.ma.where(strVarDict['TPH_GE_50'] > 0.0, strVarDict['TPHC_GE_50'] / strVarDict['TPH_GE_50'], 0.0))
                # truncate proportion of conifers to 1. Regression equations allow for this to be >1
                indexStack[len(indexStack) - 1] = np.ma.where(indexStack[len(indexStack) - 1] < 1.0, indexStack[len(indexStack) - 1], 1.0)
                logFile.write("\n\nAdding image and processing TPHC_GE_50 / TPH_GE_50")

            rastFunc.outputMaskedRaster(outputDir + "\\rtv_" + regionAbb + "_CON", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'con', headerList, area33List, area66List, printAreaHeader)

            # Build quadradic mean diameter index for the Red Tree Vole
            print "*** Processing QMD index"
            nodeList = [(0,0), (15,0), (50,1.0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['QMDA_GE_3'], 'QMDA_GE_3', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\rtv_" + regionAbb + "_QMD", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'qmd', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build canopy closure index for the Red Tree Vole
            print "*** Processing canopy closure index"
            nodeList = [(0,0), (20,0.15), (65,0.75), (96,1.0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\rtv_" + regionAbb + "_CC", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'cc', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build canopy heterogeneity index for the Red Tree Vole
            print "*** Processing canopy heterogeneity index"
            nodeList = [(0,0), (3.0,0), (4.5,0.4), (6.5,0.8), (8.0,1.0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['DDI'], 'DDI', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\rtv_" + regionAbb + "_DDI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputIndexArea2(indexStack[len(indexStack) - 1], indexAreaFile, indexAreaFile2)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.25, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\rtv_" + regionAbb + "_HCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\rtv_" + regionAbb + "_HCI_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        elif (species == 'blackbacked'):
            if (region != 'John Day'):
                errorHandler(logFile, indexAreaFile, indexAreaFile2, "**Warning. No model variation outside John Day")

            # Build large snag index for Black-backed Woodpecker
            print "*** Processing Large Snag index"
            nodeList = [(0,0), (20,1.0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STPH_GE_30'], 'STPH_GE_30', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\bbw_" + regionAbb + "_LS", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ls', headerList, area33List, area66List, printAreaHeader)

            # Build Small snag index for Black-backed Woodpecker
            print "*** Processing Small Snag index"
            nodeList = [(0,0), (100,1.0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STPH_10_30'], 'STPH_10_30', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\bbw_" + regionAbb + "_SS", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ss', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build Canopy cover index for Black-backed Woodpecker
            print "*** Processing Canopy Cover index"
            # null model - nodeList = [(0,0), (20,1.0), (70,1.0), (80,0.0)]
            nodeList = [ (0,1.0), (30,1.0), (70,0.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\bbw_" + regionAbb + "_CC", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'cc', headerList, area33List, area66List, printAreaHeader)


            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.33333, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\bbw_" + regionAbb + "_HCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\bbw_" + regionAbb + "_HCI_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        elif (species == 'whiteheaded'):
            if (region != 'John Day'):
                errorHandler(logFile, indexAreaFile, indexAreaFile2, "**Warning. No model variation outside John Day")

            # Build large snag index for White-headed Woodpecker
            print "*** Processing Large Snag index"
            nodeList = [(0,0), (10,1.0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STPH_GE_25'], 'STPH_GE_25', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\whw_" + regionAbb + "_LS", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ls', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'value', 0.6, 'multiply', logFile)

            # Build nesting canopy cover index for White-headed Woodpecker
            print "*** Processing nesting canopy cover index"
            nodeList = [(0,0), (10,1.0), (25,1.0), (35,0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\whw_" + regionAbb + "_NCC", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ncc', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'value', 0.4, 'multiply', logFile)
            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\whw_" + regionAbb + "_NCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'nci', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'copy', logFile)

            # Build live tree score (TS) for the White-headed Woodpecker
            print "*** Processing Live Tree Score"
            nodeList = [(0,0), (55,1.0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['TPH_GE_25'], 'TPH_GE_25', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\whw_" + regionAbb + "_TS", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ts', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build foraging canopy cover index for the White-headed Woodpecker
            print "*** Processing foraging canopy cover index"
            nodeList = [(0,0), (10,0), (30,1.0), (60,1.0), (80,0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\whw_" + regionAbb + "_FCC", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fcc', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.3333, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\whw_" + regionAbb + "_FI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fi', headerList, area33List, area66List, printAreaHeader)

            # Get mean score for circular neighborhood surrounding the focal pixel
            print "*** Processing Neighborhood Average for Landscape Score"
            rastFunc.circularWindowStats(indexStack, 314, referenceDesc, 'mean', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\whw_" + regionAbb + "_FCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fci', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'minimum', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\whw_" + regionAbb + "_HCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\whw_" + regionAbb + "_HCI_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        elif (species == 'olivesided'):
            # Build canopy cover index for Olive-sided Woodpecker
            print "*** Processing Canopy Closure index"

            nodeList = [(0,0), (5,1.0), (65,1.0), (100,0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\osf_" + regionAbb + "_CC", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'cc', headerList, area33List, area66List, printAreaHeader)

            # Build snag tree index for Olive-sided Woodpecker
            print "*** Processing snag tree index"

            nodeList = [(0,0), (25,1.0), (32,1.0), (65,0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STPH_GE_10'], 'STPH_GE_10', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\osf_" + regionAbb + "_ST", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'st', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build live tree index for Olive-sided Woodpecker
            print "*** Processing live tree index"

            # null hypothesis nodeList = [(0,0), (90,1.0), (105,1.0), (200,0) ]
            nodeList = [(0,0), (90,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['TPH_GE_10'], 'TPH_GE_10', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\osf_" + regionAbb + "_LT", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'lt', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.33333, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\osf_" + regionAbb + "_NCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'nci', headerList, area33List, area66List, printAreaHeader)

            # Build foraging capability index
            print "***Processing foraging capability index"
            indexStack.append(rastFunc.edgeContrastIndex(strVarDict, region, 360, outputDir + "\\osf_" + regionAbb + "_EC", referenceDesc, logFile))
            # indexStack[len(indexStack) - 1].mask = indexStack[0].mask

            rastFunc.outputMaskedRaster(outputDir + "\\osf_" + regionAbb + "_FCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fci', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'minimum', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\osf_" + regionAbb + "_HCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\osf_" + regionAbb + "_HCI_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        elif (species == 'spottedowl'):
            if (region == 'John Day'):
                errorHandler(logFile, indexAreaFile, indexAreaFile2, "**Warning. No model variation for John Day")

            habitatIndexStack = []
            # Build small tree index for Northern Spotted Owl
            print "*** Processing small tree index"
            nodeList = [(0,0.1), (80,0.1), (170,1.0), (240,1.0), (330, 0.1) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['TPH_10_25'], 'TPH_10_25', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\nso_" + regionAbb + "_smt", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'smt', headerList, area33List, area66List, printAreaHeader)

            # Build medium tree index for Northern Spotted Owl
            print "*** Processing medium tree index"
            nodeList = [(0,0.0), (20,0.0), (80,1.0), (135,1.0), (195, 0.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['TPH_25_50'], 'TPH_25_50', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\nso_" + regionAbb + "_mdt", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'mdt', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build large tree index for Northern Spotted Owl
            print "*** Processing large tree index"
            nodeList = [(0,0.0), (50,0.4), (65,1.0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['TPH_GE_75'], 'TPH_GE_75', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\nso_" + regionAbb + "_lgt", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'lgt', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build ddi index for Northern Spotted Owl
            print "*** Processing ddi index"
            nodeList = [(0,0.0), (5.5,0.2), (7.5,1.0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['DDI'], 'DDI', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\nso_" + regionAbb + "_ddi", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ddi', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.25, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\nso_" + regionAbb + "_NCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'nci', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'value', 2, 'power', logFile)

            # Build landscape radius 1 index for Northern Spotted Owl
            print "*** Calculating proportion of good habitat within 300m"
            habitatIndexStack.append(rastFunc.habitatFocalWindow(strVarDict['AGE_DOM_BA'], 80, 9999, 300, referenceDesc, logFile))
            nodeList = [(0.0,0.05), (0.4,0.05), (0.6,0.6), (0.75,1.0), (0.8,1.0), (1.0,0.6)]
            indexStack.append(rastFunc.multipleSegmentIndex(habitatIndexStack[0], 'Prop Good', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\nso_" + regionAbb + "_s1", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 's1', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'value', 3, 'power', logFile)

            # Build landscape radius 2 index for Northern Spotted Owl
            print "*** Calculating proportion of good and moderate habitat within 800m"
            habitatIndexStack.pop()
            habitatIndexStack.append(rastFunc.habitatFocalWindow(strVarDict['AGE_DOM_BA'], 80, 9999, 800, referenceDesc, logFile))
            rastFunc.mathFunc(habitatIndexStack, 'value', 0.75, 'multiply', logFile)
            habitatIndexStack.append(rastFunc.habitatFocalWindow(strVarDict['AGE_DOM_BA'], 30, 80, 800, referenceDesc, logFile))
            rastFunc.mathFunc(habitatIndexStack, 'value', 0.25, 'multiply', logFile)

            rastFunc.mathFunc(habitatIndexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(habitatIndexStack, 'value', 0.75, 'divide', logFile)

            indexStack.append(habitatIndexStack[0])
            rastFunc.outputMaskedRaster(outputDir + "\\nso_" + regionAbb + "_s2", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 's2', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'value', 2, 'power', logFile)
            rastFunc.mathFunc(indexStack, 'image', 0, 'multiply', logFile)

            # Build landscape radius 3 index for Northern Spotted Owl
            print "*** Calculating proportion of good and moderate habitat within 2400m"
            habitatIndexStack.pop()
            habitatIndexStack.append(rastFunc.habitatFocalWindow(strVarDict['AGE_DOM_BA'], 80, 9999, 2400, referenceDesc, logFile))
            rastFunc.mathFunc(habitatIndexStack, 'value', 0.75, 'multiply', logFile)
            habitatIndexStack.append(rastFunc.habitatFocalWindow(strVarDict['AGE_DOM_BA'], 30, 80, 2400, referenceDesc, logFile))
            rastFunc.mathFunc(habitatIndexStack, 'value', 0.25, 'multiply', logFile)

            rastFunc.mathFunc(habitatIndexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(habitatIndexStack, 'value', 0.75, 'divide', logFile)

            indexStack.append(habitatIndexStack[0])
            rastFunc.outputMaskedRaster(outputDir + "\\nso_" + regionAbb + "_s3", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 's3', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'multiply', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.166666667, 'power', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\nso_" + regionAbb + "_FCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fci', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'multiply', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.333333333, 'power', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\nso_" + regionAbb + "_HCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\nso_" + regionAbb + "_HCI_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        elif (species == 'marten'):
            if (region == 'Tillamook'):
                errorHandler(logFile, indexAreaFile, indexAreaFile2, "**Warning. No model variation for Tillamook")

            # Build large snag index for American Marten
            print "*** Processing large snag index"
            if (region == 'John Day'):
                nodeList = [ (0,0.0), (5,1.0) ]
            elif (region == 'West Cascades'):
                nodeList = [ (0,0.0), (3,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STPH_50_75'], 'STPH_50_75', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mrt_" + regionAbb + "_ls", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ls', headerList, area33List, area66List, printAreaHeader)

            # Build down log index for American Marten
            print "*** Processing down log index"
            if (region == 'John Day'):
                nodeList = [ (0,0.0), (2,1.0) ]
                indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['DLPH_GE_30'], 'DLPH_GE_30', nodeList, logFile))
            elif (region == 'West Cascades'):
                nodeList = [ (0,0.0), (160,1.0) ]
                indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['DVPH_GE_50'], 'DVPH_GE_50', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mrt_" + regionAbb + "_dls", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'dls', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'multiply', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.5, 'power', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mrt_" + regionAbb + "_NCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'nci', headerList, area33List, area66List, printAreaHeader)

            # Build ddi index for American Marten
            print "*** Processing ddi index"
            nodeList = [(0,0.0), (3.5,1.0)]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['DDI'], 'DDI', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mrt_" + regionAbb + "_ddi", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'ddi', headerList, area33List, area66List, printAreaHeader)

            # Get mean score for circular neighborhood surrounding the focal pixel
            print "*** Processing Neighborhood Average for Landscape Score"
            if (region == 'John Day'):
                rastFunc.circularWindowStats(indexStack, 1000, referenceDesc, 'mean', logFile)
            elif (region == 'West Cascades'):
                rastFunc.circularWindowStats(indexStack, 500, referenceDesc, 'mean', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mrt_" + regionAbb + "_FCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fci', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'multiply', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.5, 'power', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mrt_" + regionAbb + "_HCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\mrt_" + regionAbb + "_HCI_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        elif (species == 'murrelet'):
            if (region != 'Tillamook'):
                errorHandler(logFile, indexAreaFile, indexAreaFile2, "**Warning. No model variation outside Tillamook")

            # Build Stand age index for Marbled Murrelet
            print "*** Processing stand age index"
            nodeList = [ (0,0.0), (60,0.0), (175,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['AGE_DOM_BA'], 'AGE_DOM_BA', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mur_" + regionAbb + "_AGE", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'age', headerList, area33List, area66List, printAreaHeader)

            # Build Canopy Cover index for Marbled Murrelet
            print "*** Processing canopy cover index"
            nodeList = [ (0,0.0), (10,0.0), (26,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mur_" + regionAbb + "_CC", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'cc', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build Stand Height index for Marbled Murrelet
            print "*** Processing stand height index"
            nodeList = [ (0,0.0), (30,0.0), (50,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STNDHGT'], 'STNDHGT', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mur_" + regionAbb + "_HGT", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hgt', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build Mean Diameter index for Marbled Murrelet
            print "*** Processing mean diameter index"
            nodeList = [ (0,0.0), (40,0.0), (140,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['QMDA_GE_3'], 'QMDA_GE_3', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mur_" + regionAbb + "_QMD", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'qmd', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build Canopy Layers index for Marbled Murrelet
            print "*** Processing canopy layers index"
            nodeList = [ (0,0.0), (2,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['IMAP_LAYERS'], 'IMAP_LAYERS', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mur_" + regionAbb + "_LYR", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'lyr', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.2, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mur_" + regionAbb + "_HCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\mur_" + regionAbb + "_HCI_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        elif (species == 'goshawk'):
            if (region != 'John Day'):
                errorHandler(logFile, indexAreaFile, indexAreaFile2, "**Warning. No model variation outside John Day")

            # Build tph index for Goshawk
            print "*** Processing tph index"
            nodeList = [ (0,0.0), (200,0.0), (450,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['TPH_GE_3'], 'TPH_GE_3', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\gos_" + regionAbb + "_TPH", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'tph', headerList, area33List, area66List, printAreaHeader)

            # Build qmd index for Goshawk
            print "*** Processing qmd index"
            nodeList = [ (0,0.0), (12,0.0), (22,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['QMDA_GE_10'], 'QMDA_GE_10', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\gos_" + regionAbb + "_QMD", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'qmd', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build stand height index for Goshawk
            print "*** Processing stand height index"
            nodeList = [ (0,0.0), (5,0.0), (19,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['STNDHGT'], 'STNDHGT', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\gos_" + regionAbb + "_HGT", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hgt', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build canopy cover index for Goshawk
            print "*** Processing canopy cover index"
            nodeList = [ (0,0.0), (14,0.0), (48,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\gos_" + regionAbb + "_CC", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'cc', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.25, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\gos_" + regionAbb + "_HCI", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\gos_" + regionAbb + "_HCI_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        elif (species == 'muledeer'):
            # ************ Build foraging habitat suitability - size class index
            print "*** Processing Foraging habitat suitability - Size Class Index"
            nodeList = [ (0,1.0), (13,1.0), (76,0.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['AVGDBH_GE_3'], 'AVGDBH_GE_3', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_fsz", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fsz', headerList, area33List, area66List, printAreaHeader)

            # Build foraging habitat suitability - canopy cover index
            print "*** Processing Foraging habitat suitability - Canopy Cover Index"
            nodeList = [ (0,1.0), (10,1.0), (60,0.1) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_fcc", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fcc', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build foraging habitat suitability - canopy layers index
            print "*** Processing Foraging habitat suitability - Canopy Layers Index"
            nodeList = [ (0,1.0), (1,1.0), (2,0.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['IMAP_LAYERS'], 'IMAP_LAYERS', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_fly", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fly', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            if (region == 'John Day'):
                # Build foraging habitat suitability - Bitterbrush Cover
                print "*** Processing Foraging habitat suitability - Bitterbrush cover Index"
                nodeList = [ (0,0.0), (10,1.0) ]
                indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['PUTR_COV'], 'PUTR_COV', nodeList, logFile))

                rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_fbb", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
                rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fbb', headerList, area33List, area66List, printAreaHeader)

                rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
                rastFunc.mathFunc(indexStack, 'value', 0.25, 'multiply', logFile)
            else:
                rastFunc.mathFunc(indexStack, 'value', 0.333333, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_fh", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fh', headerList, area33List, area66List, printAreaHeader)

            # rastFunc.squareWindowStats(indexStack, 3, referenceDesc, 'max', logFile)
            indexStack[len(indexStack) - 1] = rastFunc.habitatFocalWindow(indexStack[len(indexStack) - 1], rastFunc.getPercentile(indexStack[len(indexStack) - 1], 75), 9999, 691, referenceDesc, logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_fhs", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'fhs', headerList, area33List, area66List, printAreaHeader)

            # weight final HSI score towards the foraging habitat by multiplying by 4
            rastFunc.mathFunc(indexStack, 'value', 4, 'multiply', logFile)
            # truncate foraging score to 0.5
            # rastFunc.mathFunc(indexStack, 'value', 0.5, 'minimum', logFile)

            # ********** Build hiding habitat suitability - size class index
            print "*** Processing hiding habitat suitability - Size Class Index"
            nodeList = [ (0,0.0), (13,1.0), (26,1.0), (50,0.25) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['AVGDBH_GE_3'], 'AVGDBH_GE_3', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_hsz", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hsz', headerList, area33List, area66List, printAreaHeader)

            # Build hiding habitat suitability - canopy cover index
            print "*** Processing hiding habitat suitability - Canopy Cover Index"
            nodeList = [ (0,0.0), (10,0.0), (60,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_hcc", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hcc', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build hiding habitat suitability - canopy layers index
            print "*** Processing hiding habitat suitability - Canopy Layers Index"
            nodeList = [ (0,0.0), (1,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['IMAP_LAYERS'], 'IMAP_LAYERS', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_hly", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hly', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            if (region == 'John Day'):
                # Build hiding habitat suitability - Bitterbrush Cover
                print "*** Processing hiding habitat suitability - Bitterbrush cover Index"
                nodeList = [ (0,0.0), (30,1.0) ]
                indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['PUTR_COV'], 'PUTR_COV', nodeList, logFile))

                rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_hbb", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
                rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hbb', headerList, area33List, area66List, printAreaHeader)

                rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
                rastFunc.mathFunc(indexStack, 'value', 0.25, 'multiply', logFile)
            else:
                rastFunc.mathFunc(indexStack, 'value', 0.333333, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_hh", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hh', headerList, area33List, area66List, printAreaHeader)

            # rastFunc.squareWindowStats(indexStack, 3, referenceDesc, 'max', logFile)
            indexStack[len(indexStack) - 1] = rastFunc.habitatFocalWindow(indexStack[len(indexStack) - 1], rastFunc.getPercentile(indexStack[len(indexStack) - 1], 75), 9999, 691, referenceDesc, logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_hhs", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hhs', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # ********** Build thermal cover suitability - size class index
            print "*** Processing thermal cover suitability - Size Class Index"
            nodeList = [ (0,0.0), (25,0.0), (38,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['AVGDBH_GE_3'], 'AVGDBH_GE_3', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_tsz", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'tsz', headerList, area33List, area66List, printAreaHeader)

            # Build thermal cover suitability - canopy cover index
            print "*** Processing thermal cover suitability - Canopy Cover Index"
            nodeList = [ (0,0.0), (60,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['CANCOV'], 'CANCOV', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_tcc", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'tcc', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)

            # Build thermal cover suitability - canopy layers index
            print "*** Processing thermal cover suitability - Canopy Layers Index"
            nodeList = [ (0,0.0), (2,1.0) ]
            indexStack.append(rastFunc.multipleSegmentIndex(strVarDict['IMAP_LAYERS'], 'IMAP_LAYERS', nodeList, logFile))

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_tly", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'tly', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(indexStack, 'value', 0.333333, 'multiply', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_tc", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'tc', headerList, area33List, area66List, printAreaHeader)

            # rastFunc.squareWindowStats(indexStack, 3, referenceDesc, 'max', logFile)
            indexStack[len(indexStack) - 1] = rastFunc.habitatFocalWindow(indexStack[len(indexStack) - 1], rastFunc.getPercentile(indexStack[len(indexStack) - 1], 75), 9999, 691, referenceDesc, logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_tcs", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'tcs', headerList, area33List, area66List, printAreaHeader)

            rastFunc.mathFunc(indexStack, 'image', 0, 'add', logFile)
            rastFunc.mathFunc(indexStack, 'value', 6, 'divide', logFile)

            rastFunc.outputMaskedRaster(outputDir + "\\mdr_" + regionAbb + "_hci", indexStack[len(indexStack) - 1], referenceDesc, logFile, 'index')
            rastFunc.outputClassRaster(outputDir + "\\mdr_" + regionAbb + "_hci_3c", indexStack[len(indexStack) - 1], referenceDesc, 0.3333, 0.6666, logFile)
            rastFunc.readIndexArea(indexStack[len(indexStack) - 1], 'hci', headerList, area33List, area66List, printAreaHeader)


        # write index area numbers to output files
        if (printAreaHeader == True):
            indexAreaFile.write('Year')
            indexAreaFile2.write('Year')

            for item in headerList:
                indexAreaFile.write(',' + item)
                indexAreaFile2.write(',' + item)

            indexAreaFile.write("\n")
            indexAreaFile2.write("\n")
            printAreaHeader = False

        indexAreaFile.write(lcYear)
        indexAreaFile2.write(lcYear)

        for value in area33List:
            indexAreaFile.write(',' + value)

        for value in area66List:
            indexAreaFile2.write(',' + value)

        indexAreaFile.write("\n")
        indexAreaFile2.write("\n")

    # write strVar averages to ouput files if strVar is being run
    if (len(averageValList) > 0):
        strVarOutputFile = open("C:\\temp\\Biodiversity\\" + regionDir + "\\TestSet.11.20\\" + species + "\\strVar_2920_3010_avg.csv", 'a')
        rastFunc.avgStrVarOuput(averageValList, strVarOutputFile, runNum)
        strVarOutputFile.close()

    logFile.write("\n\nEnd time: " + datetime.date.today().strftime("%B %d, %Y") + " - " + datetime.datetime.now().time().strftime("%H:%M:%S"))
    logFile.close()
    indexAreaFile.close()
    indexAreaFile2.close()

    # garbage collection
    del strVarDict, indexStack, preDict
    strVarDict = None
    indexStack = None
    preDict = None
    print "\nDone."


    ##
    ##import datetime
    ##from openpyxl.workbook import Workbook
    ##
    ##wb = Workbook()
    ##ws = wb.worksheets[0]
    ##
    ### set date using a Python datetime
    ##ws.cell('A1').value = datetime.datetime(2010, 7, 21)
    ##
    ##print ws.cell('A1').style.number_format.format_code # returns 'yyyy-mm-dd'
    ##
    ### set percentage using a string followed by the percent sign
    ##ws.cell('B1').value = '3.14%'
    ##
    ##print ws.cell('B1').value # returns 0.031400000000000004
    ##
    ##print ws.cell('B1').style.number_format.format_code # returns '0%'

if __name__ == '__main__':
        # Test for correct number of arguments
    if len(sys.argv) != 7:
        print "Usage: wildHCI_landcarb.py <region> <species> <runSet> <run> <rep> <year list>"
        sys.exit(1)

    try:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

    except Exception, e:
        print "\n\n" + sys.argv[2] + ": " + e.args[0]

    except:
        print "unhandled Error!!"

