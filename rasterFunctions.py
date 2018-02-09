# ---------------------------------------------------------------------------
# rasterFunctions.py
# Created on: 27 December 2012
#
# Description: These functions are required by model_strVar.py
# ---------------------------------------------------------------------------

import arcpy
import numpy as np
# import pandas as pd
from pandas.io import excel
# from openpyxl import load_workbook
from scipy.ndimage import filters
from math import pi, sqrt
# arcpy.env.overwriteOutput = "True"

nodataVal = -999

# Error handler
def errorFunctions(logFile, message):
    print message
    logFile.write("\n\n" + message)
    logFile.close()
    raise Exception(message)

# Function to load rasters and convert to masked array
def loadMaskedRaster (rasterInputName, referenceDesc, logFile):
    logFile.write("\n\nLoading: " + rasterInputName)
    rasterDesc = arcpy.Describe(rasterInputName)
    if (rasterDesc.meanCellHeight != referenceDesc.meanCellHeight):
        errorFunctions(logFile, "meanCellHeight (cell size) does not match reference.")

    if (rasterDesc.spatialReference.projectionCode != referenceDesc.spatialReference.projectionCode):
        errorFunctions(logfile, "spatialReference does not match reference. \nRaster: " + rasterDesc.spatialReference.exportToString() + "\nReference: " + referenceDesc.spatialReference.exportToString())

    if (rasterDesc.spatialReference.centralMeridian != referenceDesc.spatialReference.centralMeridian):
        errorFunctions(logfile, "centralMeridian does not match reference. \nRaster: " + rasterDesc.spatialReference.exportToString() + "\nReference: " + referenceDesc.spatialReference.exportToString())

    if (rasterDesc.spatialReference.linearUnitCode != referenceDesc.spatialReference.linearUnitCode):
        errorFunctions(logfile, "linearUnitCode does not match reference. \nRaster: " + rasterDesc.spatialReference.exportToString() + "\nReference: " + referenceDesc.spatialReference.exportToString())

    logFile.write("\nRaster matches reference projection and cell size")
    lowerLeftCorner = arcpy.Point(referenceDesc.Extent.XMin + (referenceDesc.meanCellWidth / 2), referenceDesc.Extent.YMin + (referenceDesc.meanCellHeight / 2))
    myArray = arcpy.RasterToNumPyArray(rasterInputName, lowerLeftCorner, referenceDesc.width, referenceDesc.height, nodataVal)
    logFile.write("\nMin value: " + str(myArray.min()) + " Max value: " + str(myArray.max()))
    return (np.ma.masked_values(myArray, nodataVal))

# -----------------------------------------------------------------------------------------
# write masked numpy array to raster with projection and vat
def outputMaskedRaster(outRasterName, outRaster, referenceDesc, logFile, type):
    if (type == 'index'):
        # convert float 0-1 to int 0-100
        outRaster = ((outRaster * 100) + 0.5).astype(np.int)

    logFile.write("\nWriting output raster: " + outRasterName)
    logFile.write("\nMin value: " + str(outRaster.min()) + " Max value: " + str(outRaster.max()))

    # set coordinates
    lowerLeftCorner = arcpy.Point(referenceDesc.Extent.XMin, referenceDesc.Extent.YMin)

    # convert the masked array to normal numpy array and fill the masked area with nodataVal
    rasterToSave = arcpy.NumPyArrayToRaster(outRaster.filled(nodataVal), lowerLeftCorner, referenceDesc.meanCellWidth, referenceDesc.meanCellHeight, nodataVal)
    rasterToSave.save(outRasterName)

    # define the projection based on the reference raster
    arcpy.DefineProjection_management(outRasterName, referenceDesc.spatialReference)
    # arcpy.BuildRasterAttributeTable_management(rasterOutputName, "Overwrite")

# -------------------------------------------------------------------------------------------
# write classified HCI score and write hectares to text file
def outputClassRaster(outRasterName, outRaster, referenceDesc, classBr1, classBr2, logFile):
    # classify the index raster
    outRaster = np.ma.where(outRaster >= classBr2, 3, outRaster)
    outRaster = np.ma.where( (outRaster < classBr2) & (outRaster >= classBr1), 2, outRaster)
    outRaster = np.ma.where(outRaster < classBr1, 1, outRaster)
    outRaster = outRaster.astype(np.int)

    # set coordinates
    lowerLeftCorner = arcpy.Point(referenceDesc.Extent.XMin, referenceDesc.Extent.YMin)

    # convert the masked array to normal numpy array and fill the masked area with nodataVal
    logFile.write("\nWriting output raster: " + outRasterName)
    logFile.write("\nClass raster with class break 1 (" + str(classBr1) + ") and class break 2 (" + str(classBr2) + ")")
    rasterToSave = arcpy.NumPyArrayToRaster(outRaster.filled(nodataVal), lowerLeftCorner, referenceDesc.meanCellWidth, referenceDesc.meanCellHeight, nodataVal)
    rasterToSave.save(outRasterName)

    # define the projection based on the reference raster
    arcpy.DefineProjection_management(outRasterName, referenceDesc.spatialReference)
    # arcpy.BuildRasterAttributeTable_management(rasterOutputName, "Overwrite")

# -------------------------------------------------------------------------------------------
def writeAreaFileHeader(headerString, file1, file2):
    file1.write(headerString)
    file2.write(headerString)

# -------------------------------------------------------------------------------------------
# write area of index over a given threshold to index area file
def outputIndexArea(outRaster, indexThreshold, indexAreaFile):

    # count hectares in classes and write to file
    classRas = outRaster.filled(nodataVal)
    indexAreaFile.write(',' + str(len(classRas[np.where(classRas >= indexThreshold)])))

# -------------------------------------------------------------------------------------------
# write area of index over a given threshold to index area file
def outputIndexArea2(outRaster, file1, file2):

    # count hectares in classes and write to file
    classRas = outRaster.filled(nodataVal)
    file1.write(',' + str(len(classRas[np.where(classRas >= 0.3333)])))
    file2.write(',' + str(len(classRas[np.where(classRas >= 0.6666)])))

# -------------------------------------------------------------------------------------------
# write area of index over a given threshold to index area list to be written to a file
def readIndexArea(outRaster, headerName, headerList, area33List, area66List, printAreaHeader):
    if (printAreaHeader == True):
        headerList.append(headerName)

    # count hectares in classes and write to file
    classRas = outRaster.filled(nodataVal)
    area33List.append(str(len(classRas[np.where(classRas >= 0.3333)])))
    area66List.append(str(len(classRas[np.where(classRas >= 0.6666)])))

# -------------------------------------------------------------------------------------------
# Perform index function with simple breakpoint
def singleBreakIndex(strVarDict, varName, breakValue, indexSlope, logFile):
    logFile.write("\n\nAdding image and processing " + varName + " single breakpoint index with slope: " + str(indexSlope) + " and breakpoint value: " + str(breakValue))

    if (strVarDict.has_key(varName)):
        if (indexSlope > 0):
            return(np.ma.where(strVarDict[varName] > breakValue, 1.0, strVarDict[varName] * indexSlope))
        else:
            return(np.ma.where(strVarDict[varName] > breakValue, 0, (strVarDict[varName] * indexSlope) + 1.0))
    else:
        errorFunctions(logFile, "Error!! variable " + varName + " not on stack")

# -------------------------------------------------------------------------------------------
# Perform index function with multiple inflections
def multipleSegmentIndex(strVarImage, varName, nodeList, logFile):
    logFile.write("\n\nAdding image and processing " + varName + " multiple segment index with nodes: " + str(nodeList))

    # create slope and intercept for each segment (write to log file)
    paramList = []
    for node in xrange(len(nodeList) - 1):
        slope = float(nodeList[node + 1][1] - nodeList[node][1]) / float(nodeList[node + 1][0] - nodeList[node][0])
        intercept = float(nodeList[node][1]) - float(slope * nodeList[node][0])
        paramList.append((slope, intercept))

    logFile.write("\nMultiple segment index with parameters: " + str(paramList))

    # create variable image with fcid image and variable lookup dictionary
    indexImage = np.ma.zeros(strVarImage.shape, dtype=np.float64)
    indexImage.mask = strVarImage.mask
    for y in xrange(strVarImage.shape[0]):
        for x in xrange(strVarImage.shape[1]):
            if (strVarImage.mask[y,x] == False):
                if (strVarImage.data[y,x] >= nodeList[len(nodeList) - 1][0]):
                    indexImage.data[y,x] = float(nodeList[len(nodeList) - 1][1])
                else:
                    for node in xrange(len(nodeList) - 1):
                        if strVarImage.data[y,x] >= nodeList[node][0] and strVarImage.data[y,x] < nodeList[node + 1][0]:
                            indexImage.data[y,x] = (strVarImage.data[y,x]) * paramList[node][0] + paramList[node][1]

    return(indexImage)

# -------------------------------------------------------------------------------------------
# Perform index function with simple math and single value
def mathFunc(indexStack, mathMode, value, operator, logFile):
    topIndex = len(indexStack) - 1
    if (mathMode == 'value'):
        logFile.write("\n\nProcessing image (" + str(topIndex) + "): " + operator + " with value: " + str(value))

        if (operator == 'add'):
            indexStack[topIndex] = indexStack[topIndex] + value
        elif (operator == 'multiply'):
            indexStack[topIndex] = indexStack[topIndex] * value
        elif (operator == 'divide'):
            indexStack[topIndex] = indexStack[topIndex] / value
        elif (operator == 'power'):
            indexStack[topIndex] = np.ma.power(indexStack[topIndex], value)
        elif (operator == 'minimum'):
            indexStack[topIndex] = np.ma.where(indexStack[topIndex] <= value, indexStack[topIndex], value)

    elif (mathMode == 'image'):
        logFile.write("\n\nProcessing two image math function (" + operator + ") with top images")

        if (operator == 'add'):
            indexStack[topIndex - 1] = indexStack[topIndex] + indexStack[topIndex - 1]
            indexStack.pop()
        elif (operator == 'multiply'):
            indexStack[topIndex - 1] = indexStack[topIndex] * indexStack[topIndex - 1]
            indexStack.pop()
        elif (operator == 'minimum'):
            indexStack[topIndex - 1] = np.ma.minimum(indexStack[topIndex], indexStack[topIndex - 1])
            indexStack.pop()
        elif (operator == 'copy'):
            indexStack.append(indexStack[topIndex])
        else:
            errorFunctions(logFile, "Error: mathMode operator not found in mathFunc: " + operator)

    else:
        errorFunctions(logFile, "Error: mathMode not found in mathFunc: " + mathMode)

# --------------------------------------------------------------------------------------------------------------
# find percentile from non-masked portion of raster image
def getPercentile (image, percentileValue):
    img = image.filled(nodataVal)
    return (np.percentile(img[np.where( img >= 0 )], percentileValue))

# --------------------------------------------------------------------------------------------------------------
# implement square buffer filter to generate the mean or max in a given area size
def squareWindowStats (indexStack, pixDim, referenceDesc, statName, logFile):
    topImage = len(indexStack) - 1
    saveMask = indexStack[topImage].mask

    # build kernel from pixel dimention
    kernelFP = np.array([], dtype=bool)
    for row in xrange(pixDim):
        for col in xrange(pixDim):
            kernelFP = np.append(kernelFP, [True])

    kernelFP.shape = (pixDim, pixDim)
    logFile.write("\nProcessing square window " + statName + " for " + str(pixDim) + " rows and cols")

    if (statName == 'mean'):
        outImage = filters.generic_filter(indexStack[topImage].filled(nodataVal), meanFunc, footprint=kernelFP, mode='constant', cval=nodataVal)
    elif (statName == 'max'):
        outImage = filters.generic_filter(indexStack[topImage].filled(nodataVal), maxFunc, footprint=kernelFP, mode='constant', cval=nodataVal)
    elif (statName == 'sum'):
        outImage = filters.generic_filter(indexStack[topImage].filled(nodataVal), sumFunc, footprint=kernelFP, mode='constant', cval=nodataVal)

    indexStack[topImage] = np.ma.asarray(outImage)
    indexStack[topImage].mask = saveMask

# --------------------------------------------------------------------------------------------------------------
# implement circular buffer filter to generate the mean in a given area size
def circularWindowStats (indexStack, areaHa, referenceDesc, statName, logFile):
    topImage = len(indexStack) - 1
    saveMask = indexStack[topImage].mask

    # build kernel from radius and input raster cell size
    radiusMeters = sqrt(areaHa * 10000 / pi)
    kernelFP = buildCircularKernel(radiusMeters, referenceDesc, logFile)

    logFile.write("\nProcessing circular window " + statName + " for " + str(areaHa) + " ha")

    if (statName == 'mean'):
        outImage = filters.generic_filter(indexStack[topImage].filled(nodataVal), meanFunc, footprint=kernelFP, mode='constant', cval=nodataVal)
    elif (statName == 'max'):
        outImage = filters.generic_filter(indexStack[topImage].filled(nodataVal), maxFunc, footprint=kernelFP, mode='constant', cval=nodataVal)
    elif (statName == 'sum'):
        outImage = filters.generic_filter(indexStack[topImage].filled(nodataVal), sumFunc, footprint=kernelFP, mode='constant', cval=nodataVal)

    indexStack[topImage] = np.ma.asarray(outImage)
    indexStack[topImage].mask = saveMask

# -----------------------------------------------------------------------------------------------------------------
# function to generate average value in kernel that is not NODATA
def meanFunc(kernel):
    if kernel[int(kernel.size / 2)] == nodataVal:
        return nodataVal
    else:
        # remove backgroud values from array
        return np.average(kernel[np.where( kernel != nodataVal )])

# -----------------------------------------------------------------------------------------------------------------
# function to generate average value in kernel that is not NODATA
def maxFunc(kernel):
    if kernel[int(kernel.size / 2)] == nodataVal:
        return nodataVal
    else:
        # remove backgroud values from array
        return np.max(kernel[np.where( kernel != nodataVal )])

# -----------------------------------------------------------------------------------------------------------------
# function to generate average value in kernel that is not NODATA
def sumFunc(kernel):
    if kernel[int(kernel.size / 2)] == nodataVal:
        return nodataVal
    else:
        # remove backgroud values from array
        return np.sum(kernel[np.where( kernel != nodataVal )])

# ----------------------------------------------------------------------------------------------------------------
# funciton to construct a circular kernal from a given radius
def buildCircularKernel(radiusMeters, referenceDesc, logFile):
    kernel = np.array([], dtype=bool)
    # build kernel from radius and input raster cell size
    window = int(radiusMeters / referenceDesc.meanCellHeight)
    for row in range(-window, window + 1):
        for col in range(-window, window + 1):
            Xdist = abs(col * referenceDesc.meanCellHeight)
            Ydist = abs(row * referenceDesc.meanCellHeight)
            totalDist = pow(Xdist**2 + Ydist**2, 0.5)
            if totalDist > radiusMeters:
                kernel = np.append(kernel, [False])
            else:
                kernel = np.append(kernel, [True])

    kernel.shape = (pow(kernel.size,0.5),pow(kernel.size,0.5))
    print str(kernel.size) + " pixels in circular kernel"

    logFile.write("\n\nFocal window kernel created with radius: " + str(radiusMeters) + "m / kernel total size: " + str(kernel.size) + " pixels / kernel true size: " + str(kernel[np.where(kernel == [True])].size) + " pixels")

    return(kernel)

# ----------------------------------------------------------------------------------------------------------------
# function to create image of proportion of selected habitat within a focal window
def habitatFocalWindow(strVarImage, lowerBreak, upperBreak, radiusMeters, referenceDesc, logFile):
    kernelMeanFP = buildCircularKernel(radiusMeters, referenceDesc, logFile)

    logFile.write("\nProcessing habitat focal window for radius: " + str(radiusMeters) + " and lowerBreak " + str(lowerBreak) + " / upperBreak " + str(upperBreak))
    print ("Processing habitat focal window for radius: " + str(radiusMeters) + " and lowerBreak " + str(lowerBreak) + " / upperBreak " + str(upperBreak) + "\n")

    habitatImage = np.ma.where( (strVarImage >= lowerBreak) & (strVarImage < upperBreak), 1.0, 0.0)
    outImage = filters.generic_filter(habitatImage.filled(nodataVal), meanFunc, footprint=kernelMeanFP, mode='constant', cval=nodataVal)
    outImageMasked = np.ma.asarray(outImage)
    outImageMasked.mask = strVarImage.mask
    return(outImageMasked)

# ----------------------------------------------------------------------------------------------------------------
# implement edge contrast function for Olive-sided Flycatcher
def edgeContrastIndex(strVarDict, region, radiusMeters, outFileName, referenceDesc, logFile):
    # build kernel from radius and input raster cell size
    kernelFP = np.array([[False, True, False], [True, True, True], [False, True, False]], dtype=bool)
    print str(kernelFP)

    # build kernel from radius and input raster cell size
    kernelMeanFP = buildCircularKernel(radiusMeters, referenceDesc, logFile)

    logFile.write("\n\nProcessing edge contrast function for: " + region)
    edgeDict = {}
    # all masked areas are fed into the filter with their nodata values. This forces masked areas to act as open or 'early' patches
    if (region == "John Day"):
        argsFilter = [0,0,40,100,0.016666667,-0.66666667]
        edgeImage = filters.generic_filter(strVarDict['CANCOV'].filled(nodataVal), edgeContrastFunc, footprint=kernelFP, mode='constant', cval=nodataVal, extra_arguments=(edgeDict, argsFilter))
    else:
        # null hyppthesis argsFilter = [0,0,20,90,0.0142857,-0.2857143]
        argsFilter = [0,0,20,65,0.0222222,-0.444444]
        edgeImage = filters.generic_filter(strVarDict['STNDHGT'].filled(nodataVal), edgeContrastFunc, footprint=kernelFP, mode='constant', cval=nodataVal, extra_arguments=(edgeDict, argsFilter))

    print "Edges identified: " + str(argsFilter[0]) + ", with max edge contrast value: " + str(argsFilter[1])
    logFile.write("\nEdges identified: " + str(argsFilter[0]) + ", with max edge contrast value: " + str(argsFilter[1]))
    outputMaskedRaster(outFileName, np.ma.asarray(edgeImage), referenceDesc, logFile, 'float')
    outImage = filters.generic_filter(edgeImage, edgeContrastMean, footprint=kernelMeanFP, mode='constant', cval=nodataVal, extra_arguments=(edgeDict, "touple"))
    return(np.ma.masked_where(strVarDict['STNDHGT'].mask, outImage))

# -----------------------------------------------------------------------------------------------------------------
# function to generate edge contrast index in kernel that is not NODATA
def edgeContrastMean(kernel, edgeDict, junk):
    # remove backgroud values from array
    kernelDataList = kernel[np.where( kernel != nodataVal )]
    if len(kernelDataList) > 0:
        edgeContrastValList = []
        for kernelData in kernelDataList:
            edgeContrastValList.extend(edgeDict[kernelData])

        return np.average(edgeContrastValList)
    else:
        return 0.0

# -----------------------------------------------------------------------------------------------------------------
# function to generate edge contrast index in kernel that is not NODATA
def edgeContrastFunc(kernel, edgeDict, argsFilter):
    # if center pixel is bigger than 'early' threshold, skip
    # nodata pixels are still scanned in case they have large neighbor which would indicate an edge
    # for this exercise, nodata is assumed to be open or 'early'
    if (kernel[int(kernel.size / 2)] > argsFilter[2]):
        return nodataVal
    else:
        # search early pixel for edges and record in early pixel dictionary
        edgeContrastList = []
        for adjPixel in (0, 1, 3 ,4):
            if (kernel[adjPixel] > argsFilter[3]):
                edgeContrast = 1.0
            elif (kernel[adjPixel] > argsFilter[2]):
                edgeContrast = (argsFilter[4] * kernel[adjPixel]) + argsFilter[5]
            else:
                edgeContrast = -1.0

            # edge found! add to list that will be inserted into the dictionary and calc stats
            if (edgeContrast > -1.0):
                edgeContrastList.append(edgeContrast)
                argsFilter[0] += 1
                if edgeContrast > argsFilter[1]:
                    argsFilter[1] = edgeContrast

        if (len(edgeContrastList) == 0):
            return nodataVal
        else:
            earlyPatchNum = len(edgeDict) + 1
            edgeDict[earlyPatchNum] = edgeContrastList
            return earlyPatchNum

# -------------------------------------------------------------------------------------------------
# load input predictors
def loadPredictors(landCarbDir, lcYear, referenceDesc, logFile, preDict):
    print "* Loading predictor variables"
    # load input raster into numpy array
    preDict['pRas_tBio_mgha'] = loadMaskedRaster(landCarbDir + r'\LiveMassAboveGround_' + lcYear + '.img', referenceDesc, logFile)
    preDict['pRas_sBio_mgha'] = loadMaskedRaster(landCarbDir + r'\SnagMass_' + lcYear + '.img', referenceDesc, logFile)
    preDict['pRas_lBio_mgha'] = loadMaskedRaster(landCarbDir + r'\LogMass_' + lcYear + '.img', referenceDesc, logFile)
    preDict['pRas_ageDom'] = loadMaskedRaster(landCarbDir + r'\Age_' + lcYear + '.img', referenceDesc, logFile)
    preDict['pRas_stndHgt'] = loadMaskedRaster(landCarbDir + r'\HeightUpperTree_' + lcYear + '.img', referenceDesc, logFile)

    # load static predictor variables
    #preDict['pRas_elev'] = loadMaskedRaster(staticPredictorDir + r'\dem_or100cc', referenceDesc, logFile)
    #preDict['pRas_slopeDeg'] = loadMaskedRaster(staticPredictorDir + r'\slope_deg100a', referenceDesc, logFile)
    #preDict['pRas_annPre30'] = loadMaskedRaster(staticPredictorDir + r'\annpre30_100', referenceDesc, logFile)
    # preDict['pRas_smrTp30'] = loadMaskedRaster(staticPredictorDir + r'\smrtp30_100', referenceDesc, logFile)

    # convert carbon mass from LandCarb to biomass for regression model input
    preDict['pRas_tBio_mgha'] = preDict['pRas_tBio_mgha'] * 2.0
    preDict['pRas_sBio_mgha'] = preDict['pRas_sBio_mgha'] * 2.0
    preDict['pRas_lBio_mgha'] = preDict['pRas_lBio_mgha'] * 2.0

    # truncate snag and log biomass to 1 so ln function returns valid values. Also, truncate the other predictors to 0
    if (preDict['pRas_tBio_mgha'].min() < 0):
        preDict['pRas_tBio_mgha'] = np.ma.where(preDict['pRas_tBio_mgha'] < 0, 0, preDict['pRas_tBio_mgha'])
    if (preDict['pRas_sBio_mgha'].min() < 0):
        preDict['pRas_sBio_mgha'] = np.ma.where(preDict['pRas_sBio_mgha'] < 0, 0, preDict['pRas_sBio_mgha'])
    if (preDict['pRas_lBio_mgha'].min() < 0):
        preDict['pRas_lBio_mgha'] = np.ma.where(preDict['pRas_lBio_mgha'] < 0, 0, preDict['pRas_lBio_mgha'])
    if (preDict['pRas_ageDom'].min() < 0):
        preDict['pRas_ageDom'] = np.ma.where(preDict['pRas_ageDom'] < 0, 0, preDict['pRas_ageDom'])
    if (preDict['pRas_stndHgt'].min() < 0):
        preDict['pRas_stndHgt'] = np.ma.where(preDict['pRas_stndHgt'] < 0, 0, preDict['pRas_stndHgt'])

# -------------------------------------------------------------------------------------------------
# load input GNN images
def loadGNNVars(gnnDir, varList, referenceDesc, logFile, strVarDict, landcarbMask):
    # load input raster into numpy array
    gnnFCIDimage = loadMaskedRaster(gnnDir + r'\fcid_or2_100m', referenceDesc, logFile)
    gnnFCIDimage.mask = landcarbMask

    # load GNN structure variable table into a dictionary
    for varName in varList:
        varNameGNN = varName
        #varNameGNN = varNameGNN.replace("SPH", "STPH")
        varNameGNN = varNameGNN.replace("_GT_", "_GE_")
        varNameGNN = varNameGNN.replace("LVPH", "DVPH")
        varNameGNN = varNameGNN.replace("QMDALL", "QMDA_GE_3")
        print "** Loading variable: " + varName + " (" + varNameGNN + ")"
        fcidDict = {}
        for row in arcpy.da.SearchCursor(gnnDir + r'\fcid_or2_100m', ["FCID", varNameGNN]):
            if (row[0] > 0):
                fcidDict[row[0]] = row[1]
            else:
                fcidDict[row[0]] = nodataVal

        # create variable image with fcid image and variable lookup dictionary
        strVarDict[varName] = np.ma.zeros((referenceDesc.height, referenceDesc.width), dtype=np.float64)
        strVarDict[varName].mask = landcarbMask
        for y in xrange(referenceDesc.height):
            for x in xrange(referenceDesc.width):
                if (gnnFCIDimage.data[y,x] != nodataVal):
                    strVarDict[varName].data[y,x] = fcidDict[gnnFCIDimage.data[y,x]]

        # set missing or nonforest data in GNN tables to nodata
        strVarDict[varName].mask = np.where(strVarDict[varName].data < -998.0, True, landcarbMask)

# ----------------------------------------------------------------------------------------------------
# calculate structure variables from regression equations
def calcStructureVars(varList, strVarDict, preDict, region, outputDir, referenceDesc, logFile):
    for varName in varList:
        # run logistic regression models and output presence (1) and absence (0) rasters for each variable
        logFile.write("\n\n** Processing variable: " + varName)
        print "** Processing variable: " + varName

        wb = excel.read_excel(r'R:\Biodiversity\JD_EC_WC_TK_logistic_models_January2014.xlsx', sheetname='Logistic_mods')
        wb = wb.fillna(0.0)
        coeffList = wb[(wb['AREA'] == region) & (wb['DEPVAR'] == varName + '_pa')]

        if (len(coeffList) == 1):
            logFile.write("\nRead logistic coefficients (" + str(coeffList.iloc[0]))

            # Add 1.0 to log equation to avoid error (R. Pabst - 6 Jan 2014)
            presenceTerm = np.ma.zeros((referenceDesc.height, referenceDesc.width), dtype=np.float64)
            presenceTerm += coeffList.iloc[0]['INTERCEPT']
            presenceTerm += coeffList.iloc[0]['sqrtTBIO_MGHA']*np.ma.sqrt(preDict['pRas_tBio_mgha'])
            presenceTerm += coeffList.iloc[0]['sqrtSBIO_MGHA']*np.ma.sqrt(preDict['pRas_sBio_mgha'])
            presenceTerm += coeffList.iloc[0]['lnSBIO_MGHA']*np.ma.log(preDict['pRas_sBio_mgha'] + 1.0)
            presenceTerm += coeffList.iloc[0]['lnLBIO_MGHA']*np.ma.log(preDict['pRas_lBio_mgha'] + 1.0)
            presenceTerm += coeffList.iloc[0]['sqrtAGE_DOM_BA']*np.ma.sqrt(preDict['pRas_ageDom'])
            presenceTerm += coeffList.iloc[0]['STNDHGT']*preDict['pRas_stndHgt']
            presenceTerm += coeffList.iloc[0]['pow2STNDHGT']*np.ma.power(preDict['pRas_stndHgt'], 2.0)

            #presenceTerm += coeffList[8]*np.ma.sqrt(preDict['pRas_stndHgt'])
            #presenceTerm += coeffList[9]*preDict['pRas_elev']
            #presenceTerm += coeffList[10]*np.ma.sqrt(preDict['pRas_slopeDeg'])
            #presenceTerm += coeffList[11]*preDict['pRas_annPre30']

            presenceRas = 1.0 / (1.0 + np.ma.exp(presenceTerm * -1.0))
            # outputMaskedRaster(outputDir + "\\" + varName + "pr", presenceRas, referenceDesc, logFile)

            presenceRas = np.ma.where(presenceRas >= 0.5, 1.0, 0.0)
            # outputMaskedRaster(outputDir + "\\" + varName + "pa", presenceRas.astype(np.int), referenceDesc, logFile)

        elif (len(coeffList) == 0):
            presenceRas = np.ma.ones((referenceDesc.height, referenceDesc.width), dtype=np.float64)
        else:
            errorFunctions(logFile, "Error. logistic regression coefficients duplicated")

        # run linear regression models to predict forest structure variable

        wb = excel.read_excel(r'R:\Biodiversity\JD_EC_WC_TK_linear_models_January26_2014.xlsx', sheetname='Linear_mods', skiprows=2)
        wb = wb.fillna(0.0)
        coeffList = wb[wb['Area'] == region]

        coeffRow = []
        for row in range(0, len(coeffList)):
            if (varName in coeffList.iloc[row]['DEPVAR']):
                if (coeffList.iloc[row]['DEPVAR'][0:2] == 'ln'):
                    varMath = 'log'
                    tableVarName = coeffList.iloc[row]['DEPVAR'][2:]
                elif (coeffList.iloc[row]['DEPVAR'][0:4] == 'sqrt'):
                    varMath = 'sqrt'
                    tableVarName = coeffList.iloc[row]['DEPVAR'][4:]
                elif (coeffList.iloc[row]['DEPVAR'][0:4] == 'pow2'):
                    varMath = 'pow2'
                    tableVarName = coeffList.iloc[row]['DEPVAR'][4:]
                else:
                    varMath = 'none'
                    tableVarName = coeffList.iloc[row]['DEPVAR']

                if (tableVarName == varName):
                    coeffRow = coeffList.iloc[row]
                    # exit for loop early
                    break


        if (len(coeffRow) > 0):
            logFile.write("\nRead linear regression coefficients (" + str(coeffRow))

            # Add 1.0 to log equation to avoid error (R. Pabst - 6 Jan 2014)
            vegTerm = np.ma.zeros((referenceDesc.height, referenceDesc.width), dtype=np.float64)
            vegTerm += coeffRow['Intercept']
            vegTerm += coeffRow['sqrtTBIO_MGHA']*np.ma.sqrt(preDict['pRas_tBio_mgha'])
            vegTerm += coeffRow['pow2TBIO_MGHA']*np.ma.power(preDict['pRas_tBio_mgha'], 2.0)
            vegTerm += coeffRow['sqrtSBIO_MGHA']*np.ma.sqrt(preDict['pRas_sBio_mgha'])
            vegTerm += coeffRow['lnSBIO_MGHA']*np.ma.log(preDict['pRas_sBio_mgha'] + 1.0)
            vegTerm += coeffRow['lnLBIO_MGHA']*np.ma.log(preDict['pRas_lBio_mgha'] + 1.0)
            vegTerm += coeffRow['sqrtAGE_DOM_BA']*np.ma.sqrt(preDict['pRas_ageDom'])
            vegTerm += coeffRow['STNDHGT']*preDict['pRas_stndHgt']
            vegTerm += coeffRow['pow2STNDHGT']*np.ma.power(preDict['pRas_stndHgt'], 2.0)
            vegTerm += coeffRow['sqrtSTNDHGT']*np.ma.sqrt(preDict['pRas_stndHgt'])
##            vegTerm += coeffList[9]*preDict['pRas_elev']
##            vegTerm += coeffList[10]*np.ma.sqrt(preDict['pRas_elev'])
##            vegTerm += coeffList[11]*np.ma.sqrt(preDict['pRas_slopeDeg'])
##            vegTerm += coeffList[12]*preDict['pRas_annPre30']
##            vegTerm += coeffList[13]*preDict['pRas_smrTp30']

            # resolve the dependant variable
            if (varMath == 'log'):
                logFile.write("\nprocessing expotential function for veg structure")
                vegStrRas = np.ma.exp(vegTerm)
            elif (varMath == 'sqrt'):
                logFile.write("\nprocessing power (2) function for veg structure")
                vegStrRas = np.ma.power(vegTerm, 2)
            elif (varMath == 'pow2'):
                logFile.write("\nprocessing square root function for veg structure")
                vegStrRas = np.ma.sqrt(vegTerm)
            else: vegStrRas = vegTerm

            # multiply by presenceRas to incorporate the presence/asbsence model and truncate to zero
            vegStructureRas = vegStrRas * presenceRas
            if (vegStructureRas.min() < 0):
                vegStructureRas = np.ma.where(vegStructureRas < 0, 0, vegStructureRas)

            if (varName == "CANCOV"):
                # fix for regression equation so high biomass doesn't create low cover. - R. Pabst 20 July 2016
                strVarDict[varName] = np.ma.where(preDict['pRas_tBio_mgha'] > 550, 100, vegStructureRas)
            else:
                strVarDict[varName] = vegStructureRas
            # outputMaskedRaster(outputDir + "\\" + varName, vegStructureRas, referenceDesc, logFile, 'float')

        elif (varName == "STNDHGT"):
            strVarDict[varName] = preDict['pRas_stndHgt']
        elif (varName == "AGE_DOM_BA"):
            strVarDict[varName] = preDict['pRas_ageDom']
        elif (varName == "LC_TBIO_MGHA"):
            strVarDict[varName] = preDict['pRas_tBio_mgha']
        elif (varName == "LC_SBIO_MGHA"):
            strVarDict[varName] = preDict['pRas_sBio_mgha']
        elif (varName == "LC_LBIO_MGHA"):
            strVarDict[varName] = preDict['pRas_lBio_mgha']

        else:
            errorFunctions(logFile, "Error. linear regression coefficients not found or duplicated")

        del wb
        wb = None

#------------------------------------------------------------------------------------------------------
def avgStrVarOuput(averageValList, strVarOutputFile, runNum):
    avgList = []
    for fieldNum in range(len(averageValList[0])):
        avgList.append(averageValList[0][fieldNum][1])

    for rowNum in range(1, len(averageValList)):
        for fieldNum in range(len(averageValList[rowNum])):
            if (averageValList[0][fieldNum][0] == averageValList[rowNum][fieldNum][0]):
                avgList[fieldNum] += averageValList[rowNum][fieldNum][1]
            else:
                print "Error!!!" + str(averageValList[rowNum][fieldNum][0]) + " does not match"
                exit(1)

    strVarOutputFile.write(runNum)

    for fieldNum in range(len(averageValList[0])):
        strVarOutputFile.write("," + averageValList[0][fieldNum][0] + "," + str(avgList[fieldNum] / len(averageValList)))

    strVarOutputFile.write("\n")
