import methodLibary as ml
import os
import sys
import argparse


def diseaseClusterNetwork(directory, tsi, pValue):
    """
    A protocol for an individual network by disease clusters
    """
    print("\n --------- diseaseclusterNetwork Protocol ------------ \n")

    if not os.path.exists(directory):
        os.makedirs(directory)
    newDir = directory + ("/diseaseClusterNetwork")

    if not os.path.exists(newDir):
        os.makedirs(newDir)

    allDiseaseMirnaInteractions = ml.readMiriadData()
    matrix = ml.Matrix(allDiseaseMirnaInteractions)
    numberOfClusters = matrix.createClusters()
    mirnaTissueInteractions = ml.readTissueMatrix(ml.readMirnaTSI(tsi))

    #network for all clusters
    print("\n All Clusters \n")
    diseaseMirnaInteractions = matrix.getClusterInteractions()
    mirnaInteractions = ml.mirnaFullNames(matrix.getClusterRows())
    mirnaNames = [interaction[1] for interaction in mirnaInteractions]
    mirnaTargetInterations = ml.readMirTarBase(mirnaNames)
    #create edgeList
    edgeList = ml.EdgeList(mirnaInteractions)
    edgeList.addInteractions(diseaseMirnaInteractions)
    edgeList.addInteractions(mirnaTargetInterations)
    edgeList.addInteractions(mirnaTissueInteractions)
    # Lists for node selection in cytoscape
    edgeList.writeLists(newDir + "/list")
    # reduce network to only through interactions
    edgeList.removeIncompleteLinks() ##################################
    #edgeList.simplifyMirnaNames()
    #All groups
    edgeList.writeEdgeList(newDir + "/EdgeListAll.txt")

    # No Tissues
    tissues = edgeList.popInteractions()

    # No Pathways
    pathways = edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoTargets.txt")

    # No Diseases
    edgeList.popInteractions()
    edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.addInteractions(pathways)
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoDiseases.txt")

    # network for each cluster
    allClusterEdgeList = open(newDir + "/allClusterEdgeList.txt", "w")
    allClusterEdgeList.write("Source\tTarget\n")
    for i in range(numberOfClusters):
        print("\n Cluster " + str(i + 1) + "\n")
        diseaseMirnaInteractions = matrix.getClusterInteractions(i)
        mirnaInteractions = ml.mirnaFullNames(matrix.getClusterRows(i))
        mirnaNames = [interaction[1] for interaction in mirnaInteractions]
        mirnaTargetInterations = ml.readMirTarBase(mirnaNames)
        #create edgeList
        edgeList = ml.EdgeList(mirnaInteractions)
        edgeList.addInteractions(mirnaTargetInterations)
        edgeList.addInteractions(mirnaTissueInteractions)
        edgeList.addInteractions(diseaseMirnaInteractions)
        # reduce network to only through interactions
        edgeList.removeIncompleteLinks() ################################
        #edgeList.simplifyMirnaNames()
        #All groups
        allClusterEdgeList.write(edgeList.writeEdgeList(newDir + "/EdgeList"
        + str(i + 1) + ".txt"))

        # No Tissues
        tissues = edgeList.popInteractions()
        edgeList.writeEdgeList(newDir + "/EdgeListNoTissues" + str(i + 1)
        + ".txt")

        # No Pathways
        pathways = edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.writeEdgeList(newDir + "/EdgeListNoTargets" + str(i + 1)
        + ".txt")

        # No Diseases
        edgeList.popInteractions()
        edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.addInteractions(pathways)
        edgeList.writeEdgeList(newDir +
        "/EdgeListDiseaseNoDiseases"
        + str(i + 1) + ".txt")

    allClusterEdgeList.close()

def tissueClusterNetwork(directory, tsi, pValue):
    """
    A protocol for an individual network by tissue clusters
    """
    print("\n --------- tissueClusterNetwork Protocol ------------ \n")

    if not os.path.exists(directory):
        os.makedirs(directory)
    newDir = directory + ("/tissueClusterNetwork")

    if not os.path.exists(newDir):
        os.makedirs(newDir)

    allTissueMirnaInteractions = ml.readTissueMatrix(ml.readMirnaTSI(tsi))
    matrix = ml.Matrix(allTissueMirnaInteractions)
    numberOfClusters = matrix.createClusters()
    mirnaDiseaseInteractions = ml.readMiriadData()

    #network for all clusters
    print("\n All Clusters \n")
    tissueMirnaInteractions = matrix.getClusterInteractions()
    mirnaInteractions = ml.simplifyMirnaNames(matrix.getClusterRows())
    mirnaTargetInterations = ml.readMirTarBase(matrix.getClusterRows())
    #create edgeList
    edgeList = ml.EdgeList(mirnaInteractions)
    edgeList.addInteractions(tissueMirnaInteractions)
    edgeList.addInteractions(mirnaTargetInterations)
    edgeList.addInteractions(mirnaDiseaseInteractions)
    # Lists for node selection in cytoscape
    edgeList.writeLists(newDir + "/list")
    # reduce network to only through interactions
    edgeList.removeIncompleteLinks()##############################
    #edgeList.simplifyMirnaNames()
    #All groups
    edgeList.writeEdgeList(newDir + "/EdgeListAll.txt")

    # No Tissues
    tissues = edgeList.popInteractions()
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoTissues.txt")

    # No Pathways
    pathways = edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoTargets.txt")

    # No Diseases
    edgeList.popInteractions()
    edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.addInteractions(pathways)
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoDiseases.txt")

    # network for each cluster
    allClusterEdgeList = open(newDir + "/allClusterEdgeList.txt", "w")
    allClusterEdgeList.write("Source\tTarget\n")
    for i in range(numberOfClusters):
        print("\n Cluster " + str(i + 1) + "\n")
        tissueMirnaInteractions = matrix.getClusterInteractions(i)
        mirnaInteractions = ml.simplifyMirnaNames(matrix.getClusterRows(i))
        mirnaTargetInterations = ml.readMirTarBase(matrix.getClusterRows(i))
        #create edgeList
        edgeList = ml.EdgeList(mirnaInteractions)
        edgeList.addInteractions(mirnaTargetInterations)
        edgeList.addInteractions(mirnaDiseaseInteractions)
        edgeList.addInteractions(tissueMirnaInteractions)
        # reduce network to only through interactions
        edgeList.removeIncompleteLinks() ###########################
        #edgeList.simplifyMirnaNames()
        #All groups
        allClusterEdgeList.write(edgeList.writeEdgeList(newDir + "/EdgeList" + str(i + 1) + ".txt"))

        # No Tissues
        tissues = edgeList.popInteractions()
        edgeList.writeEdgeList(newDir + "/EdgeListNoTissues" + str(i + 1)
        + ".txt")

        # No Pathways
        pathways = edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.writeEdgeList(newDir + "/EdgeListNoTargets" + str(i + 1)
        + ".txt")

        # No Diseases
        edgeList.popInteractions()
        edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.addInteractions(pathways)
        edgeList.writeEdgeList(newDir + "/EdgeListNoDiseases"
        + str(i + 1) + ".txt")

    allClusterEdgeList.close()


def targetClusterNetwork(directory, tsi, pValue):
    """
    A protocol for an individual network by target clusters
    """
    print("\n --------- targetClusterNetwork Protocol ------------ \n")

    if not os.path.exists(directory):
        os.makedirs(directory)
    newDir = directory + ("/targetClusterNetwork")

    if not os.path.exists(newDir):
        os.makedirs(newDir)

    allMirnaTargetInteractions = ml.readMirTarBase()
    matrix = ml.Matrix(allMirnaTargetInteractions)
    numberOfClusters = matrix.createClusters()
    mirnaTissueInteractions = ml.readTissueMatrix(ml.readMirnaTSI(tsi))
    mirnaDiseaseInteractions = ml.readMiriadData()

    #network for all clusters
    print("\n All Clusters \n")
    targetMirnaInteractions = matrix.getClusterInteractions()
    mirnaInteractions = ml.simplifyMirnaNames(matrix.getClusterRows())
    #create edgeList
    edgeList = ml.EdgeList(mirnaInteractions)
    edgeList.addInteractions(targetMirnaInteractions)
    edgeList.addInteractions(mirnaDiseaseInteractions)
    edgeList.addInteractions(mirnaTissueInteractions)
    # Lists for node selection in cytoscape
    edgeList.writeLists(newDir + "/list")
    # reduce network to only through interactions
    edgeList.removeIncompleteLinks()###############################
    #edgeList.simplifyMirnaNames()
    #All groups
    edgeList.writeEdgeList(newDir + "/EdgeListAll.txt")

    # No Tissues
    tissues = edgeList.popInteractions()
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoTissues.txt")

    # No Pathways
    pathways = edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.writeEdgeList(newDir + "/EdgeAllNoTargets.txt")

    # No Diseases
    edgeList.popInteractions()
    edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.addInteractions(pathways)
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoDiseases.txt")

    # network for each cluster
    allClusterEdgeList = open(newDir + "/allClusterEdgeList.txt", "w")
    allClusterEdgeList.write("Source\tTarget\n")
    for i in range(numberOfClusters):
        print("\n Cluster " + str(i + 1) + "\n")
        targetMirnaInteractions = matrix.getClusterInteractions(i)
        mirnaInteractions = ml.simplifyMirnaNames(matrix.getClusterRows(i))
        #create edgeList
        edgeList = ml.EdgeList(mirnaInteractions)
        edgeList.addInteractions(mirnaDiseaseInteractions)
        edgeList.addInteractions(mirnaTissueInteractions)
        edgeList.addInteractions(targetMirnaInteractions)
        # reduce network to only through interactions
        edgeList.removeIncompleteLinks() ####
        #edgeList.simplifyMirnaNames()

        #All groups
        allClusterEdgeList.write(edgeList.writeEdgeList(newDir + "/EdgeList" + str(i + 1)
        + ".txt"))

        # No Tissues
        tissues = edgeList.popInteractions()
        edgeList.writeEdgeList(newDir + "/EdgeListNoTissues" + str(i + 1)
        + ".txt")

        # No Pathways
        pathways = edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.writeEdgeList(newDir + "/EdgeListNoTargets" + str(i + 1)
        + ".txt")

        # No Diseases
        edgeList.popInteractions()
        edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.addInteractions(pathways)
        edgeList.writeEdgeList(newDir +
        "/EdgeListNoDiseases"
        + str(i + 1) + ".txt")

    allClusterEdgeList.close()


def diseaseClusterNetworkPA(directory, tsi, pValue):
    """
    A protocol for an individual network by disease clusters using reacotme
    pathways
    """
    print("\n --------- diseaseclusterNetworkPA Protocol ------------ \n")

    if not os.path.exists(directory):
        os.makedirs(directory)
    newDir = directory + ("/diseaseClusterNetworkPA")
    if not os.path.exists(newDir):
        os.makedirs(newDir)

    allDiseaseMirnaInteractions = ml.readMiriadData()
    matrix = ml.Matrix(allDiseaseMirnaInteractions)
    numberOfClusters = matrix.createClusters()
    mirnaTissueInteractions = ml.readTissueMatrix(ml.readMirnaTSI(tsi))

    #network for all clusters
    print("\n All Clusters \n")
    diseaseMirnaInteractions = matrix.getClusterInteractions()
    mirnaInteractions = ml.mirnaFullNames(matrix.getClusterRows())
    mirnaNames = [interaction[1] for interaction in mirnaInteractions]
    mirnaTargetInterations = ml.readMirTarBase(mirnaNames, useGeneClass=True)
    mirnaPathwayInterations = ml.reactomePA(mirnaTargetInterations, pValue)
    #create edgeList
    edgeList = ml.EdgeList(mirnaInteractions)
    edgeList.addInteractions(diseaseMirnaInteractions)
    edgeList.addInteractions(mirnaPathwayInterations)
    edgeList.addInteractions(mirnaTissueInteractions)
    # Lists for node selection in cytoscape
    edgeList.writeLists(newDir + "/list")
    # reduce network to only through interactions
    edgeList.removeIncompleteLinks()
    #edgeList.simplifyMirnaNames()
    #All groups
    edgeList.writeEdgeList(newDir + "/EdgeListAll.txt")

    # No Tissues
    tissues = edgeList.popInteractions()

    # No Pathways
    pathways = edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoTargets.txt")

    # No Diseases
    edgeList.popInteractions()
    edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.addInteractions(pathways)
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoDiseases.txt")

    # network for each cluster
    for i in range(numberOfClusters):
        print("\n Cluster " + str(i + 1) + "\n")
        diseaseMirnaInteractions = matrix.getClusterInteractions(i)
        mirnaInteractions = ml.mirnaFullNames(matrix.getClusterRows(i))
        mirnaNames = [interaction[1] for interaction in mirnaInteractions]
        mirnaTargetInterations = ml.readMirTarBase(mirnaNames)
        #create edgeList
        edgeList = ml.EdgeList(mirnaInteractions)
        edgeList.addInteractions(mirnaPathwayInterations)
        edgeList.addInteractions(mirnaTissueInteractions)
        edgeList.addInteractions(diseaseMirnaInteractions)
        # reduce network to only through interactions
        edgeList.removeIncompleteLinks()
        #edgeList.simplifyMirnaNames()
        #All groups
        edgeList.writeEdgeList(newDir + "/EdgeList" + str(i + 1)
        + ".txt")

        # No Tissues
        tissues = edgeList.popInteractions()
        edgeList.writeEdgeList(newDir + "/EdgeListNoTissues" + str(i + 1)
        + ".txt")

        # No Pathways
        pathways = edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.writeEdgeList(newDir + "/EdgeListNoTargets" + str(i + 1)
        + ".txt")

        # No Diseases
        edgeList.popInteractions()
        edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.addInteractions(pathways)
        edgeList.writeEdgeList(newDir +
        "/EdgeListDiseaseNoDiseases"
        + str(i + 1) + ".txt")


def tissueClusterNetworkPA(directory, tsi, pValue):
    """
    A protocol for an individual network by disease clusters using reacotme
    pathways
    """
    print("\n --------- tissueClusterNetworkPA Protocol ------------ \n")

    if not os.path.exists(directory):
        os.makedirs(directory)
    newDir = directory + ("/tissueClusterNetworkPA")

    if not os.path.exists(newDir):
        os.makedirs(newDir)

    allTissueMirnaInteractions = ml.readTissueMatrix(ml.readMirnaTSI(tsi))
    matrix = ml.Matrix(allTissueMirnaInteractions)
    numberOfClusters = matrix.createClusters()
    mirnaDiseaseInteractions = ml.readMiriadData()

    #network for all clusters
    print("\n All Clusters \n")
    tissueMirnaInteractions = matrix.getClusterInteractions()
    mirnaInteractions = ml.simplifyMirnaNames(matrix.getClusterRows())
    mirnaTargetInterations = ml.readMirTarBase(matrix.getClusterRows(),
         useGeneClass=True)
    mirnaPathwayInterations = ml.reactomePA(mirnaTargetInterations, pValue)
    #create edgeList
    edgeList = ml.EdgeList(mirnaInteractions)
    edgeList.addInteractions(tissueMirnaInteractions)
    edgeList.addInteractions(mirnaPathwayInterations)
    edgeList.addInteractions(mirnaDiseaseInteractions)
    # Lists for node selection in cytoscape
    edgeList.writeLists(newDir + "/list")
    # reduce network to only through interactions
    edgeList.removeIncompleteLinks()
    #edgeList.simplifyMirnaNames()
    #All groups
    edgeList.writeEdgeList(newDir + "/EdgeListAll.txt")

    # No Tissues
    tissues = edgeList.popInteractions()
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoTissues.txt")

    # No Pathways
    pathways = edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoTargets.txt")

    # No Diseases
    edgeList.popInteractions()
    edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.addInteractions(pathways)
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoDiseases.txt")

    # network for each cluster
    for i in range(numberOfClusters):
        print("\n Cluster " + str(i + 1) + "\n")
        tissueMirnaInteractions = matrix.getClusterInteractions(i)
        mirnaInteractions = ml.simplifyMirnaNames(matrix.getClusterRows(i))
        mirnaTargetInterations = ml.readMirTarBase(matrix.getClusterRows(i))
        #create edgeList
        edgeList = ml.EdgeList(mirnaInteractions)
        edgeList.addInteractions(mirnaPathwayInterations)
        edgeList.addInteractions(mirnaDiseaseInteractions)
        edgeList.addInteractions(tissueMirnaInteractions)
        # reduce network to only through interactions
        edgeList.removeIncompleteLinks()
        #edgeList.simplifyMirnaNames()
        #All groups
        edgeList.writeEdgeList(newDir + "/EdgeList" + str(i + 1) + ".txt")

        # No Tissues
        tissues = edgeList.popInteractions()
        edgeList.writeEdgeList(newDir + "/EdgeListNoTissues" + str(i + 1)
        + ".txt")

        # No Pathways
        pathways = edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.writeEdgeList(newDir + "/EdgeListNoTargets" + str(i + 1)
        + ".txt")

        # No Diseases
        edgeList.popInteractions()
        edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.addInteractions(pathways)
        edgeList.writeEdgeList(newDir + "/EdgeListNoDiseases"
        + str(i + 1) + ".txt")


def targetClusterNetworkPA(directory, tsi, pValue):
    """
    A protocol for an individual network by disease clusters using reacotme
    pathways
    """
    print("\n --------- targetClusterNetworkPA Protocol ------------ \n")

    if not os.path.exists(directory):
        os.makedirs(directory)
    newDir = directory + ("/targetClusterNetworkPA")

    if not os.path.exists(newDir):
        os.makedirs(newDir)

    allMirnaTargetInteractions = ml.readMirTarBase()
    matrix = ml.Matrix(allMirnaTargetInteractions)
    numberOfClusters = matrix.createClusters()
    mirnaTissueInteractions = ml.readTissueMatrix(ml.readMirnaTSI(tsi))
    mirnaDiseaseInteractions = ml.readMiriadData()

    #network for all clusters
    print("\n All Clusters \n")
    targetMirnaInteractions = matrix.getClusterInteractions()

    entrezDict = ml.geneNameEntrezDict()
    newTargetMirnaInteractions = []
    for interaction in targetMirnaInteractions:
        newTargetMirnaInteractions.append([interaction[0],
        ml.Gene(interaction[1], entrezDict[interaction[1]])])

    mirnaPathwayInterations = ml.reactomePA(newTargetMirnaInteractions, pValue)
    mirnaInteractions = ml.simplifyMirnaNames(matrix.getClusterRows())
    #create edgeList
    edgeList = ml.EdgeList(mirnaInteractions)
    edgeList.addInteractions(mirnaPathwayInterations)
    edgeList.addInteractions(mirnaDiseaseInteractions)
    edgeList.addInteractions(mirnaTissueInteractions)
    # Lists for node selection in cytoscape
    edgeList.writeLists(newDir + "/list")
    # reduce network to only through interactions
    edgeList.removeIncompleteLinks()
    #edgeList.simplifyMirnaNames()
    #All groups
    edgeList.writeEdgeList(newDir + "/EdgeListAll.txt")

    # No Tissues
    tissues = edgeList.popInteractions()
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoTissues.txt")

    # No Pathways
    pathways = edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.writeEdgeList(newDir + "/EdgeAllNoTargets.txt")

    # No Diseases
    edgeList.popInteractions()
    edgeList.popInteractions()
    edgeList.addInteractions(tissues)
    edgeList.addInteractions(pathways)
    edgeList.writeEdgeList(newDir + "/EdgeListAllNoDiseases.txt")

    # network for each cluster
    for i in range(numberOfClusters):
        print("\n Cluster " + str(i + 1) + "\n")
        targetMirnaInteractions = matrix.getClusterInteractions(i)

        entrezDict = ml.geneNameEntrezDict()
        newTargetMirnaInteractions = []
        for interaction in targetMirnaInteractions:
            newTargetMirnaInteractions.append([interaction[0],
            ml.Gene(interaction[1], entrezDict[interaction[1]])])
        mirnaPathwayInterations = ml.reactomePA(newTargetMirnaInteractions, pValue)
        mirnaInteractions = ml.simplifyMirnaNames(matrix.getClusterRows(i))
        #create edgeList
        edgeList = ml.EdgeList(mirnaInteractions)
        edgeList.addInteractions(mirnaDiseaseInteractions)
        edgeList.addInteractions(mirnaTissueInteractions)
        edgeList.addInteractions(mirnaPathwayInterations)
        # reduce network to only through interactions
        edgeList.removeIncompleteLinks()
        #edgeList.simplifyMirnaNames()

        #All groups
        edgeList.writeEdgeList(newDir + "/EdgeList" + str(i + 1)
        + ".txt")

        # No Tissues
        tissues = edgeList.popInteractions()
        edgeList.writeEdgeList(newDir + "/EdgeListNoTissues" + str(i + 1)
        + ".txt")

        # No Pathways
        pathways = edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.writeEdgeList(newDir + "/EdgeListNoTargets" + str(i + 1)
        + ".txt")

        # No Diseases
        edgeList.popInteractions()
        edgeList.popInteractions()
        edgeList.addInteractions(tissues)
        edgeList.addInteractions(pathways)
        edgeList.writeEdgeList(newDir +
        "/EdgeListNoDiseases"
        + str(i + 1) + ".txt")


def tissueTargetCircos(directory, tsi, empty):
    print("\n --------- tissueTargetCircos Protocol  ------------ \n")
    """
    Create edgelist for network clustered by tissue-miRNA.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

    allMirnaTissueInteractions = ml.readTissueMatrix(ml.readMirnaTSI(tsi))
    tissueMatrix = ml.Matrix(allMirnaTissueInteractions)
    tisClustnum = tissueMatrix.createClusters()

    allMirnaTargetInteractions = ml.readMirTarBase()
    targetMatrix = ml.Matrix(allMirnaTargetInteractions)
    tarClustnum = targetMatrix.createClusters()

    num = 0
    interactionDict = {}

    mirnaSet = set()
    tissueSet = set()
    targetSet = set()
    for i in range(tisClustnum):
        mirnaTissueInteractions = tissueMatrix.getClusterInteractions(i)
        tissueClusterName = ""
        for name in tissueMatrix.getClusterColumns(i):
            tissueClusterName += name
        for j in range(tarClustnum):
            mirnaTargetInteractions = targetMatrix.getClusterInteractions(j)
            #print(mirnaTargetInteractions)
            targetClusterName = ""
            for name in targetMatrix.getClusterColumns(j):
                targetClusterName += name
            print("\n")
            for tarItems, tisItems in [(tarItems, tisItems) for tarItems in
            mirnaTargetInteractions for tisItems in mirnaTissueInteractions]:
                #print(tisItems)
                if tarItems[0] == tisItems[0]:
                    print(tarItems[0])
                    #print("+")
                    #print("\n\nCluster: ", i + 1)
                    #print("miRNAs: ", + len(tissueMatrix.getClusterRows(i)))
                    #print("tissues: ", + len(tissueMatrix.getClusterColumns(i)))
                    tissueClusterSize = (len(tissueMatrix.getClusterRows(i)) + len(tissueMatrix.getClusterColumns(i)))
                    targetClusterSize = (len(targetMatrix.getClusterRows(j)) + len(targetMatrix.getClusterColumns(j)))
                    key = ("tissueCluster" + str(i + 1), "targetCluster"
                    + str(j + 1), tissueClusterSize, targetClusterSize)
                    #print("\n\nCluster: ", i + 1, "   Size: ", tissueClusterSize)
                    if key not in interactionDict:
                        interactionDict[key] = 1
                    else:
                        interactionDict[key] += 1
                    num += 1
                    for mirna in targetMatrix.getClusterRows(i):
                        mirnaSet.add(mirna)
                    for mirna in tissueMatrix.getClusterRows(i):
                        mirnaSet.add(mirna)
                    for target in targetMatrix.getClusterColumns(i):
                        targetSet.add(target)
                    for tissue in tissueMatrix.getClusterColumns(i):
                        tissueSet.add(tissue)
    print("tissueTarget mirna total: ", len(mirnaSet))
    print("tissueTarget target total: ", len(targetSet))
    print("tissueTarget tissue total: ", len(tissueSet), "\n")
    interactions = []
    for key, value in list(interactionDict.items()):
        interactions.append(key + (value,))
        print(key, "\t", value)
    edgeList = ml.EdgeList([])
    edgeList.addInteractions(interactions)
    edgeList.circos(directory + "/tissueTarget")


def tissueDiseaseCircos(directory, tsi, empty):
    print("\n --------- tissueDiseaseCircos Protocol ------------ \n")
    """
    Create edgelist for network clustered by tissue-miRNA.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

    allMirnaTissueInteractions = ml.readTissueMatrix(ml.readMirnaTSI(tsi))
    tissueMatrix = ml.Matrix(allMirnaTissueInteractions)
    tisClustnum = tissueMatrix.createClusters()

    allMirnaDiseaseInteractions = ml.readMiriadData()
    diseaseMatrix = ml.Matrix(allMirnaDiseaseInteractions)
    diseaseClustnum = diseaseMatrix.createClusters()
    interactionDict = {}
    interactions = []
    num = 0
    mirnaSet = set()
    diseaseSet = set()
    tissueSet = set()
    for i in range(tisClustnum):
        mirnaTissueInteractions = tissueMatrix.getClusterInteractions(i)
        for j in range(diseaseClustnum):
            mirnaDiseaseInteractions = diseaseMatrix.getClusterInteractions(j)
            for disItems, tisItems in [(disItems, tisItems) for disItems in
            mirnaDiseaseInteractions for tisItems in mirnaTissueInteractions]:
                if disItems[0] == ml.simplifyMirnaNames([tisItems[0]])[0][0]:
                    key = ("tissueCluster" + str(i + 1), "diseaseCluster"
                    + str(j + 1), len(tissueMatrix.getClusterRows(i))
                    + len(tissueMatrix.getClusterColumns(i)),
                    len(diseaseMatrix.getClusterRows(j))
                    + len(diseaseMatrix.getClusterColumns(j)))
                    if key not in interactionDict:
                        interactionDict[key] = 1
                    else:
                        interactionDict[key] += 1
                    num += 1

                    for mirna in tissueMatrix.getClusterRows(i):
                        mirnaSet.add(mirna)
                    for mirna in diseaseMatrix.getClusterRows(i):
                        mirnaSet.add(mirna)
                    for tissue in tissueMatrix.getClusterColumns(i):
                        #print(tissue)
                        tissueSet.add(tissue)
                    for disease in diseaseMatrix.getClusterColumns(i):
                        #print(disease)
                        diseaseSet.add(disease)
    print("tissueDisease mirna total: ", len(mirnaSet))
    print("tissueDisease tissue total: ", len(tissueSet))
    print("tissueDisease disease total: ", len(diseaseSet), "\n")
    interactions = []
    for key, value in list(interactionDict.items()):
        interactions.append(key + (value,))
        #print(key, "\t", value)
    edgeList = ml.EdgeList([])
    edgeList.addInteractions(interactions)
    edgeList.circos(directory + "/tissueDisease")


def targetDiseaseCircos(directory, empty, empty2):
    print("\n ---------targetDiseaseCircos Protocol ------------ \n")
    """
    Create edgelist for network clustered by tissue-miRNA.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

    allMirnaTargetInteractions = ml.readMirTarBase()
    targetMatrix = ml.Matrix(allMirnaTargetInteractions)
    tarClustnum = targetMatrix.createClusters()

    allMirnaDiseaseInteractions = ml.readMiriadData()
    diseaseMatrix = ml.Matrix(allMirnaDiseaseInteractions)
    diseaseClustnum = diseaseMatrix.createClusters()
    interactionDict = {}
    num = 0
    mirnaSet = set()
    diseaseSet = set()
    targetSet = set()
    for i in range(tarClustnum):
        mirnaTargetInteractions = targetMatrix.getClusterInteractions(i)
        for j in range(diseaseClustnum):
            mirnaDiseaseInteractions = diseaseMatrix.getClusterInteractions(j)
            for disItems, tarItems in [(disItems, tarItems) for disItems in
            mirnaDiseaseInteractions for tarItems in mirnaTargetInteractions]:
                if disItems[0] == ml.simplifyMirnaNames([tarItems[0]])[0][0]:
                    key = ("targetCluster" + str(i + 1), "diseaseCluster"
                    + str(j + 1), len(targetMatrix.getClusterRows(i))
                    + len(targetMatrix.getClusterColumns(i)),
                     len(diseaseMatrix.getClusterRows(j))
                     + len(diseaseMatrix.getClusterColumns(j)))
                    if key not in interactionDict:
                        interactionDict[key] = 1
                    else:
                        interactionDict[key] += 1
                    num += 1

                    for mirna in targetMatrix.getClusterRows(i):
                        mirnaSet.add(mirna)
                    for mirna in diseaseMatrix.getClusterRows(i):
                        mirnaSet.add(mirna)
                    for target in targetMatrix.getClusterColumns(i):
                        targetSet.add(target)
                    for disease in diseaseMatrix.getClusterColumns(i):
                        diseaseSet.add(disease)
    print("targetDisease mirna total: ", len(mirnaSet))
    print("targetDisease target total: ", len(targetSet))
    print("targetDisease disease total: ", len(diseaseSet), "\n")
    interactions = []
    for key, value in list(interactionDict.items()):
        interactions.append(key + (value,))
        #print(key, "\t", value)
    edgeList = ml.EdgeList([])
    edgeList.addInteractions(interactions)
    edgeList.circos(directory + "/targetDisease")


def clusterProtocol(directory, tsi, empty):

    print("Cluster lists Protocol")

    clusterDirectory = directory + "/clusters"

    if not os.path.exists(clusterDirectory):
        os.makedirs(clusterDirectory)

    allMirnaTargetInteractions = ml.readMirTarBase()
    #print(allMirnaTargetInteractions)
    targetMatrix = ml.Matrix(allMirnaTargetInteractions)
    #print(targetMatrix.getColumns())
    tarClustnum = targetMatrix.createClusters()
    #print(targetMatrix.getClusterColumns())

    allMirnaDiseaseInteractions = ml.readMiriadData()
    diseaseMatrix = ml.Matrix(allMirnaDiseaseInteractions)
    disClustnum = diseaseMatrix.createClusters()

    allMirnaTissueInteractions = ml.readTissueMatrix(ml.readMirnaTSI(tsi))
    tissueMatrix = ml.Matrix(allMirnaTissueInteractions)
    tisClustnum = tissueMatrix.createClusters()

    targetDirectory = clusterDirectory + "/target"
    if not os.path.exists(targetDirectory):
        os.makedirs(targetDirectory)
    allInteractions = []
    clusterTxtList = open(clusterDirectory + "/targetClusterList.txt", "w")
    for i in range(tarClustnum):
        edgeList = ml.EdgeList([])
        edgeList.addInteractions(targetMatrix.getClusterInteractions(i))
        edgeList.writeEdgeList(targetDirectory + "/targetCluster" + str(i + 1)
        + ".txt")
        #print(targetMatrix.getClusterInteractions(i))
        allInteractions.extend(targetMatrix.getClusterInteractions(i))
        clusterTxtList.write("Target Cluster " + str(i + 1) + "\n\n")
        targetList = []
        for target in targetMatrix.getClusterColumns(i):
            targetList.append(target)
            #print(target)
        targetList.sort()
        for target in targetList:
            clusterTxtList.write(target + "\n")
        clusterTxtList.write("\n")
        mirnaList = []
        for mirna in targetMatrix.getClusterRows(i):
            mirnaList.append(mirna)
        mirnaList.sort()
        for mirna in mirnaList:
            clusterTxtList.write(mirna + "\n")
        clusterTxtList.write("\n\n")
    clusterTxtList.close()

    targetList = ml.EdgeList([])
    targetList.addInteractions(allInteractions)
    targetList.writeLists(targetDirectory + "/targetList", False)

    tissueDirectory = clusterDirectory + "/tissue"
    allInteractions = []
    clusterTxtList = open(clusterDirectory + "/tissueClusterList.txt", "w")
    if not os.path.exists(tissueDirectory):
        os.makedirs(tissueDirectory)
    for i in range(tisClustnum):
        edgeList = ml.EdgeList([])
        edgeList.addInteractions(tissueMatrix.getClusterInteractions(i))
        edgeList.writeEdgeList(tissueDirectory + "/tissueCluster" + str(i + 1)
        + ".txt")
        allInteractions.extend(tissueMatrix.getClusterInteractions(i))
        clusterTxtList.write("Tissue Cluster " + str(i + 1) + "\n\n")
        tissueList = []
        for tissue in tissueMatrix.getClusterColumns(i):
            tissueList.append(tissue)
        tissueList.sort()
        for tissue in tissueList:
            clusterTxtList.write(tissue + "\n")
        clusterTxtList.write("\n")
        mirnaList = []
        for mirna in tissueMatrix.getClusterRows(i):
            mirnaList.append(mirna)
        mirnaList.sort()
        for mirna in mirnaList:
            clusterTxtList.write(mirna + "\n")
        clusterTxtList.write("\n\n")
    clusterTxtList.close()

    #print(allInteractions)
    tissueList = ml.EdgeList([])
    tissueList.addInteractions(allInteractions)
    tissueList.writeLists(tissueDirectory + "/tissueList", False)

    diseaseDirectory = clusterDirectory + "/disease"
    allInteractions = []
    clusterTxtList = open(clusterDirectory + "/diseaseClusterList.txt", "w")
    if not os.path.exists(diseaseDirectory):
        os.makedirs(diseaseDirectory)
    for i in range(disClustnum):
        edgeList = ml.EdgeList([])
        edgeList.addInteractions(diseaseMatrix.getClusterInteractions(i))
        edgeList.writeEdgeList(diseaseDirectory + "/diseaseCluster" + str(i + 1)
        + ".txt")
        allInteractions.extend(diseaseMatrix.getClusterInteractions(i))
        clusterTxtList.write("Disease Cluster " + str(i + 1) + "\n\n")
        diseaseList = []
        for disease in diseaseMatrix.getClusterColumns(i):
            diseaseList.append(disease)
        diseaseList.sort()
        for disease in diseaseList:
            clusterTxtList.write(disease + "\n")
        clusterTxtList.write("\n")
        mirnaList = []
        for mirna in diseaseMatrix.getClusterRows(i):
            mirnaList.append(mirna)
        mirnaList.sort()
        for mirna in mirnaList:
            clusterTxtList.write(mirna + "\n")
        clusterTxtList.write("\n\n")
    clusterTxtList.close()

    diseaseList = ml.EdgeList([])
    diseaseList.addInteractions(allInteractions)
    diseaseList.writeLists(diseaseDirectory + "/diseaseList", False)

    pathwayDirectory = clusterDirectory + "/pathway"
    if not os.path.exists(pathwayDirectory):
        os.makedirs(pathwayDirectory)
    allInteractions = []
    clusterTxtList = open(clusterDirectory + "/pathwayClusterList.txt", "w")
    for i in range(tarClustnum):
        targetMirnaInteractions = targetMatrix.getClusterInteractions(i)
        entrezDict = ml.geneNameEntrezDict()
        newTargetMirnaInteractions = []
        for interaction in targetMirnaInteractions:
            newTargetMirnaInteractions.append([interaction[0],
            ml.Gene(interaction[1], entrezDict[interaction[1]])])
        mirnaPathwayInterations = ml.reactomePA(newTargetMirnaInteractions, 0.05)
        edgeList = ml.EdgeList([])
        edgeList.addInteractions(mirnaPathwayInterations)
        edgeList.writeEdgeList(targetDirectory + "/pathwayCluster" + str(i + 1)
        + ".txt")
        allInteractions.extend(targetMatrix.getClusterInteractions(i))
        clusterTxtList.write("Pathway Cluster " + str(i + 1) + "\n\n")
        pathwaySet = set()
        mirnaSet = set()
        for interaction in mirnaPathwayInterations:
            pathwaySet.add(interaction[1])
            mirnaSet.add(interaction[0])
        pathwayList = sorted(list(pathwaySet))
        mirnaList = sorted(list(mirnaSet))

        for pathway in pathwayList:
            clusterTxtList.write(pathway + "\n")
        clusterTxtList.write("\n")

        for mirna in mirnaList:
            clusterTxtList.write(mirna + "\n")
        clusterTxtList.write("\n\n")
    clusterTxtList.close()





    """
    pathwayDirectory = directory + "/pathway"
    if not os.path.exists(pathwayDirectory):
        os.makedirs(pathwayDirectory)
    allInteractions = []
    for i in range(tarClustnum):
        entrezDict = ml.geneNameEntrezDict()
        geneInteractions = []
        for item in (targetMatrix.getClusterInteractions(i)):
            geneInteractions.append([item[0], ml.Gene(item[1],
            entrezDict[item[1]])])

        edgeList = ml.EdgeList([])
        edgeList.addInteractions(ml.reactomePA(geneInteractions, 0.03))
        edgeList.writeEdgeList(pathwayDirectory, "/pathwayCluster" + str(i + 1)
        + ".txt")
        allInteractions.extend(targetMatrix.getClusterInteractions(i))
    #print(allInteractions)
    targetList = ml.EdgeList([])
    targetList.addInteractions(allInteractions)
    targetList.writeLists(targetDirectory + "/pathwayList", False)
    """


def mirnaNumber(empty, tsi, empty2):
    mirnas = set()

    print("cluster mirnaTargets")
    allMirnaTargetInteractions = ml.readMirTarBase()
    targetMatrix = ml.Matrix(allMirnaTargetInteractions)
    tarClustnum = targetMatrix.createClusters()

    print("cluster mirnaDiseases")
    allMirnaDiseaseInteractions = ml.readMiriadData()
    diseaseMatrix = ml.Matrix(allMirnaDiseaseInteractions)
    disClustnum = diseaseMatrix.createClusters()

    print("cluster mirnaTissues")
    allMirnaTissueInteractions = ml.readTissueMatrix(ml.readMirnaTSI(tsi))
    tissueMatrix = ml.Matrix(allMirnaTissueInteractions)
    tisClustnum = tissueMatrix.createClusters()

    print("Calculate")
    count = 0
    for i in range(tisClustnum):
        clusterMirnas = tissueMatrix.getClusterRows(i)
        for mirna in clusterMirnas:
            mirnas.add(mirna)
            count += 1
    print("number of miRNAs in tissue clusters: ", count)

    tissues = set()
    for i in range(tisClustnum):
        clusterMirnas = tissueMatrix.getClusterColumns(i)
        for tissue in clusterMirnas:
            tissues.add(tissue)
    print("number of tissues: ", len(tissues))

    count = 0
    for i in range(tarClustnum):
        clusterMirnas = targetMatrix.getClusterRows(i)
        for mirna in clusterMirnas:
            mirnas.add(mirna)
            count += 1
    print("number of miRNAs in target clusters: ", count)

    targets = set()
    for i in range(tarClustnum):
        clusterMirnas = targetMatrix.getClusterColumns(i)
        for target in clusterMirnas:
            targets.add(target)
    print("number of targets: ", len(targets))

    count = 0
    for i in range(disClustnum):
        clusterMirnas = diseaseMatrix.getClusterRows(i)
        for mirna in clusterMirnas:
            count += 1
    print("number of miRNAs in disease clusters: ", count)

    diseases = set()
    for i in range(disClustnum):
        clusterMirnas = diseaseMatrix.getClusterColumns(i)
        for disease in clusterMirnas:
            diseases.add(disease)
    print("number of diseases: ", len(diseases))

    count = 0
    unusedClusters = [31, 34, 37, 42, 44, 45, 46, 48, 49, 54, 55, 57, 59, 60,
    61, 62, 63, 64, 65]
    disMirnas = set()
    for i in range(disClustnum):
        if i + 1 not in unusedClusters:
            clusterMirnas = diseaseMatrix.getClusterRows(i)
            for mirna in clusterMirnas:
                disMirnas.add(mirna)
                count += 1
    print("number of miRNAs used in disease clusters: ", count)

    simpleMirnas = [m[0] for m in ml.simplifyMirnaNames(mirnas)]

    totalMirnas = []
    for disMirna in disMirnas:
        if disMirna not in simpleMirnas:
            totalMirnas.append(disMirna)

    totalMirnas.extend(list(mirnas))

    print("number total of mirnas used: " + str(len(totalMirnas)))

    #check through miRNAs
    tarMirnas = targetMatrix.getClusterRows()
    tarTisMirnas = []
    tarTisDisMirnas = []
    for mirna in tarMirnas:
        if mirna in tissueMatrix.getClusterRows():
            tarTisMirnas.append(mirna)
    print(len(tarTisMirnas))
    simpleMirnas = [m[0] for m in ml.simplifyMirnaNames(tarTisMirnas)]
    print(simpleMirnas)
    for m in diseaseMatrix.getClusterRows():
        if m in simpleMirnas:
            tarTisDisMirnas.append(m)
    print("\nNumber of through miRNAs: ", len(tarTisDisMirnas))
    print(tarTisDisMirnas)

    tisMirnas = []
    tissues = set()
    for i in range(tisClustnum):
        clusterMirnas = tissueMatrix.getClusterRows(i)
        for mirna in clusterMirnas:
            if mirna in tarTisMirnas:
                tisMirnas.extend(clusterMirnas)
                print(i + 1)
                for tissue in tissueMatrix.getClusterColumns(i):
                    tissues.add(tissue)
                break
    print("\nNumber of cluster miRNAs from through miRNAs: ", len(tisMirnas))
    print("Number of cluster tissues from through miRNAs: ", len(tissues))


    tarMirnas = []
    targets = set()
    for i in range(tisClustnum):
        clusterMirnas = targetMatrix.getClusterRows(i)
        for mirna in clusterMirnas:
            if mirna in tarTisMirnas:
                tarMirnas.extend(clusterMirnas)
                print(i + 1)
                for target in targetMatrix.getClusterColumns(i):
                    targets.add(target)
                break
    print("\nNumber of cluster miRNAs from through miRNAs: ", len(tarMirnas))
    print("Number of cluster targets from through miRNAs: ", len(targets))


    disMirnas = []
    diseases = set()
    #simpleMirnas = [m[0] for m in ml.simplifyMirnaNames(tarTisMirnas)]
    for i in range(tisClustnum):
        clusterMirnas = diseaseMatrix.getClusterRows(i)
        for mirna in clusterMirnas:
            if mirna in tarTisDisMirnas:
                    disMirnas.extend(clusterMirnas)
                    print(i + 1)
                    for disease in diseaseMatrix.getClusterColumns(i):
                        diseases.add(disease)
                    break
    print("\nNumber of cluster miRNAs from through miRNAs: ", len(disMirnas))
    print("Number of cluster diseases from through miRNAs: ", len(diseases))

    totalMirnas = set()
    totalMirnas.update(tarMirnas)
    print("size1: ", len(totalMirnas))
    totalMirnas.update(tisMirnas)
    print("size2: ", len(totalMirnas))
    final = set()

    simpleMirnas = [m[0] for m in ml.simplifyMirnaNames(totalMirnas)]
    print(simpleMirnas)
    for mirna in disMirnas:
            if mirna not in simpleMirnas:
                final.add(mirna)
            #print(m)
    print("size3: ", len(final))
    final.update(totalMirnas)

    print("\ntotal miRNAs: ", len(final))
    #print(final)

    #check how prevalent dis/tis/tar are

    print("\n\n")
    tissueDict = {}
    for i in range(tisClustnum):
        for tissue in tissueMatrix.getClusterColumns(i):
            tissueDict[tissue] = tissueDict.setdefault(tissue, 0) + 1
    for key in sorted(tissueDict, key=tissueDict.get):
        print(key, ": ", tissueDict[key])

    print("\n\n")
    diseaseDict = {}
    for i in range(disClustnum):
        for disease in diseaseMatrix.getClusterColumns(i):
            diseaseDict[disease] = diseaseDict.setdefault(disease, 0) + 1
    for key in sorted(diseaseDict, key=diseaseDict.get):
        print(key, ": ", diseaseDict[key])

    print("\n\n")
    targetDict = {}
    for i in range(tarClustnum):
        for target in targetMatrix.getClusterColumns(i):
            targetDict[target] = targetDict.setdefault(target, 0) + 1
    for key in sorted(targetDict, key=targetDict.get):
        print(key, ": ", targetDict[key])



################
    for key in sorted(diseaseDict, key=diseaseDict.get):
        print("\n\n" + key)
        for i in range(disClustnum):
            for disease in diseaseMatrix.getClusterColumns(i):
                if key == disease:
                    print("cluster: ", i + 1)

    for key in sorted(tissueDict, key=tissueDict.get):
         print("\n\n" + key)
         for i in range(tisClustnum):
            for tissue in tissueMatrix.getClusterColumns(i):
                if tissue == key:
                     print("cluster: ", i + 1)

    for key in sorted(targetDict, key=targetDict.get):
         print("\n\n" + key)
         for i in range(tarClustnum):
            for target in targetMatrix.getClusterColumns(i):
                if target == key:
                     print("cluster: ", i + 1)

def main():
    protocolDict = {"1": tissueTargetCircos, "2": tissueDiseaseCircos,
                    "3": targetDiseaseCircos, "4": clusterProtocol,
                    "5": diseaseClusterNetwork, "6": tissueClusterNetwork,
                    "7": targetClusterNetwork, "8": diseaseClusterNetworkPA,
                    "9": tissueClusterNetworkPA, "10": targetClusterNetworkPA,
                    "11": mirnaNumber}

    parser = argparse.ArgumentParser(description=
    "This program contains methods for creating miRNA networks with disease,\n"
    "tissue and target interactions. To use this program specify the protocol\n"
    "numbers to use and follow with optional arguments where required.\n\n"
    "Choose one or more protocols from the following list:\n"
    "1 :  tissueTargetCircos\n"
    "2 :  tissueDiseaseCircos\n"
    "3 :  targetDiseaseCircos\n"
    "4 :  clusterProtocol\n"
    "5 :  diseaseclusterNetwork\n"
    "6 :  tissueClusterNetwork\n"
    "7 :  targetClusterNetwork\n"
    "8 :  diseaseClusterNetworkPA\n"
    "9 :  tissueClusterNetworkPA\n"
    "10:  targetClusterNetworkPA\n\n"
    "To run all protocols type 'all'/n",
    epilog="References:\n\n"
    "mirTarBase\n"
    "miriad\n"
    "miTissueAtlas\n",
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("protocols", nargs='*',
         help="""Choose your desired protocols from the list above, seperate
         multiple protocols witha space""")

    parser.add_argument(
         "-d", "--directory", default="../output",
         help="""Choose your directory in relation to the directory of this
         program. Deafulat directory is ../output""")
    parser.add_argument(
         "-t", "--miRNA_TSI", type=float, default=0.85000001, #maybe change to 0.85 - includes one more miRNA
         help="""Choose cut-off value for miRNA tissue specificity index.
         Default is 0.85""")
    parser.add_argument("-p", "--p_value", type=float, default=0.005,
         help="""Choose p-value for reactome pathway. Default is 0.005""")

    args = parser.parse_args()

    if len(args.protocols) < 1:
        print("""\nPlease specify a protocol.
        Use -h flag for more information\n""")
        sys.exit()

    for arg in args.protocols:
        try:
            print("\nProtocol: ", arg, "\n")
            protocolDict[arg](args.directory, args.miRNA_TSI, args.p_value)
        except KeyError:
            if arg == "all" or arg == "All":
                for key, value in sorted(list(protocolDict.items())):
                    value(args.directory, args.miRNA_TSI, args.p_value)
            else:
                print("No corresponding protocol")


if __name__ == "__main__":
    main()


