import numpy as np
import subprocess
import os
import statistics as st

class Cluster():
    """
    Represent a cluster found from a Matrix object using BCrepBimax() in R
    """
    def __init__(self, rows, columns):
        super(Cluster, self).__init__()
        self.rows = rows
        self.columns = columns

    def getRows(self):
        return self.rows

    def getColumns(self):
        return self.columns

    def compare(self, other):
        """
        Compare if another cluster is the same
        """
        if (len(self.rows) == len(other.getRows()) and len(self.columns) ==
        len(other.columns)):
            if (len(set(self.rows).intersection(other.getRows())) ==
            len(self.rows)):
                if (len(set(self.columns).intersection(other.getColumns())) ==
                 len(self.columns)):
                    return True
        return False


class Matrix():
    """
    Representation of matrix.
    Pass in a double list of interations [rows, columns] as argument
    """
    def __init__(self, interactions):
        super(Matrix, self).__init__()
        self.interactions = interactions
        self.rowList = self.createRowList()
        self.columnList = self.createColumnList()
        self.matrix = self.createMatrix()
        self.clusters = []

    def getRows(self):
        return self.rowList

    def getColumns(self):
        return self.columnList

    def getInteractions(self):
        return self.interations

    def getClusters(self):
        return self.clusters

    def createRowList(self):
        """
        Row list serves as labels and index for matrix rows
        """
        rows = set()
        for items in self.interactions:
            rows.add(items[0])
        return sorted(list(rows))

    def createColumnList(self):
        """
        Column list serves as labels and index for matrix columns
        """
        columns = set()
        for items in self.interactions:
            columns.add(items[1])
        return sorted(list(columns))

    def createMatrix(self):
        """Creates the a matrix object as a 2d array in relation to rows and
        columns"""
        # the binary matrix to set
        matrix = np.zeros((len(self.rowList), len(self.columnList)))
        #set correct 1s in binary matirx
        for items in self.interactions:
            i = 0
            j = 0
            while items[0] != self.rowList[i]:
                i += 1
            while items[1] != self.columnList[j]:
                j += 1
            matrix[i][j] = 1
        return matrix

    def writeMatirx(self, fileName):
        """
        Writes the matrix to a csv file with row and column labels
        """
        output = open(fileName, "w")
        #write binary martrix to file
        output.write("Key")
        for c in self.columnList:
            output.write("," + c)
        output.write("\n")
        #write each row of miRNA and matrix
        for i in range(len(self.rowList)):
            output.write(self.rowList[i])
            for value in self.matrix[i]:
                output.write("," + str(int(value)))
            output.write("\n")

    def createClusters(self):
        """
        Use R script to find clusters in matrix and create cluster objects.
        """
        matrixFileName = "tempMatrix.txt"
        clusterFileName = "tempClusters.txt"
        self.writeMatirx(matrixFileName)
        # Define command and arguments
        command = 'Rscript'
        path2script = 'simCluster.R'
        # Build subprocess command
        cmd = [command, path2script] + [matrixFileName, clusterFileName]
        # check_output will run the command and store to result
        x = subprocess.check_output(cmd, universal_newlines=True)
        print(x)

        # read in clusters from txt file created by R script
        clusterFile = open(clusterFileName, "r")
        i = 0
        numberOfClusters = 0
        for line in clusterFile:
            if i % 3 == 2:
                rows = line.strip().split(" ")
            if i % 3 == 0 and i > 0:
                columns = line.strip().split(" ")
                for j in range(len(columns)):
                    columns[j] = columns[j].replace(".", " ")#.replace("X", "")
                newCluster = Cluster(rows, columns)
                self.clusters.append(newCluster)
                numberOfClusters += 1
            i += 1
        clusterFile.close()
        os.remove(clusterFileName)
        os.remove(matrixFileName)
        return numberOfClusters

    def getClusterInteractions(self, clusterChoice="All"):
        """
        Return list of interactions from a given cluster or all clusters
        """
        clusterRows = set()
        clusterColumns = set()

        # Collect items from clusters
        # Defualt to crete from all clusters
        if clusterChoice == "All" or clusterChoice == "all":
            for clust in self.clusters:
                for row in clust.getRows():
                    clusterRows.add(row)
                for column in clust.getColumns():
                    clusterColumns.add(column)
        # Or user choose a particular cluster
        else:
            for row in self.clusters[clusterChoice].getRows():
                clusterRows.add(row)
            for column in self.clusters[clusterChoice].getColumns():
                clusterColumns.add(column)

        clusterInteractions = set()
        for items in self.interactions:
            if items[0] in clusterRows and items[1] in clusterColumns:
                clusterInteractions.add((items[0], items[1]))

        return list(clusterInteractions)

    def getClusterRows(self, clusterChoice="All"):
        """Create a txt file of the rows found in clusters"""
        rowItems = set()

        # Collect items from clusters
        # Defualt to crete from all clusters
        if clusterChoice == "All" or clusterChoice == "all":
            for clust in self.clusters:
                for row in clust.getRows():
                    rowItems.add(row)
        else:
            for row in self.clusters[clusterChoice].getRows():
                rowItems.add(row)

        return rowItems

    def getClusterColumns(self, clusterChoice="All"):
        """Create a txt file of the columns found in clusters"""
        columnItems = set()

        # Collect items from clusters
        # Defualt to crete from all clusters
        if clusterChoice == "All" or clusterChoice == "all":
            for clust in self.clusters:
                for column in clust.getColumns():
                    columnItems.add(column)
        else:
            for column in self.clusters[clusterChoice].getColumns():
                columnItems.add(column)

        return columnItems

    def shuffleMatrix(self):
        """Randomly shuffle matrix values. Does not preserve correct values"""
        np.random.shuffle(self.matrix)
        for row in self.matrix:
            np.random.shuffle(row)


class EdgeList():
    """
    Class to build and write an edgeList for network
    """
    def __init__(self, mirnas):
        super(EdgeList, self).__init__()
        self.mirnas = mirnas
        self.interactions = []

    def getInteractions(self):
        return self.interactions

    def readEdgeList(self, fileName):
        """
        Read miRNA-target edgelist.
        """
        db = open(fileName, "r")
        for line in db:
            elements = line.strip().split("\t")
            if elements[0] != "Source" and elements[1] != "Target":
                self.interactions.append([elements[0], elements[1]])
        db.close()

    def addInteractions(self, newInteractions):
        """
        Add a set of interactions to edge list
        """
        try:
            self.interactions.append(sorted(newInteractions))
        except TypeError:
            self.interactions.append(newInteractions)

    def removeIncompleteLinks(self):
        """
        Delete nodes in network that do not have a paths to all ends of the
        network.
        """
        # check if each mirna is in each interaction list and remove if not
        for interactionList in self.interactions:
            tempList = []
            for mirnaLink in self.mirnas:
                for interaction in interactionList:
                    if mirnaLink[0] in interaction or mirnaLink[1] in \
                    interaction:
                        tempList.append(mirnaLink)
                        break
            self.mirnas = tempList
        # check if each interaction is in the mirna list and remove if not
        reducedInteractions = []
        for interactionList in self.interactions:
            tempList = []
            for interaction in interactionList:
                for mirnaLink in self.mirnas:
                    if interaction[0] in mirnaLink or interaction[1] \
                    in mirnaLink:
                        tempList.append(interaction)
                        break
            reducedInteractions.append(tempList)
        self.interactions = reducedInteractions

    def simplifyMirnaNames(self):
        simplified = []
        for interactionList in self.interactions:
            tempList = []
            for interaction in interactionList:
                if len(interaction[0].split("-")) > 2:
                    tempList.append([interaction[0].split("-")[1] + "-" +
                    ''.join([n for n in interaction[0].split("-")[2]
                    if n in '1234567890']), interaction[1]])
                else:
                    tempList.append(interaction)
            simplified.append(tempList)
        self.interactions = simplified

    def writeEdgeList(self, fileName, writeMirnaInteractions=True):
        """
        Write the current edgeList to file
        """
        interactionText = ""
        if len([items for interactionList in self.interactions for items in
        interactionList]) > 0:
            edgeListFile = open(fileName, "w")
            edgeListFile.write("Source\tTarget\n")
            for interactionList in self.interactions:
                #print("\n\n\nInteraction list: ", interactionList)
                for items in interactionList:
                    interaction = str(items[0]) + "\t" + str(items[1]) + "\n"
                    edgeListFile.write(interaction)
                    interactionText += interaction
            if writeMirnaInteractions:
                for items in self.mirnas:
                    #print(items)
                    interaction = items[0] + "\t" + items[1] + "\n"
                    edgeListFile.write(interaction)
                    interactionText += interaction
            edgeListFile.close()
        return interactionText

    def circos(self, dirName):
        """
        create ciros images
        """
        confText = """karyotype = map.txt
chromosomes_display_default = yes
show_ticks = no
show_tick_labels = no
#chromosomes_scale   = /./=1rn
<links>
<link>
file = links.txt
radius = 0.99r
bezier_radius = 0r
color = black_a4
thickness = 2
</link>
</links>
<ideogram>
<spacing>
default = 0.005r
</spacing>
radius = 0.85r
thickness = 20p
fill = yes
stroke_color = dgrey
stroke_thickness = 2p
show_label = yes
label_font = default
label_radius = 1.01r
label_size = 15
label_parallel = no
</ideogram>
<image>
<<include /etc/circos/etc/image.conf>>
</image>
<<include /etc/circos/etc/colors_fonts_patterns.conf>>
<<include /etc/circos/etc/housekeeping.conf>>"""
        linkText = """ribbon           = yes
flat   = no
color            = black
stroke_color     = vdgrey
stroke_thickness = 1
thickness        = 1
radius           = 0.40r
bezier_radius    = 0r
crest                = 0.5
bezier_radius_purity = 0.75"""

        for interactionList in self.interactions:
            if not os.path.exists(dirName):
                os.makedirs(dirName)
            if not os.path.exists(dirName + "/etc/tracks"):
                os.makedirs(dirName + "/etc/tracks")
            linkConf = open(dirName + "/etc/tracks/link.conf", "w")
            linkConf.write(linkText)
            linkConf.close()
            circosConf = open(dirName + "/circos.conf", "w")
            circosConf.write(confText)
            circosConf.close()
            mapFile = open(dirName + "/map.txt", "w")
            linksFile = open(dirName + "/links.txt", "w")
            edgeListFile = open(dirName + "/edgeList.txt", "w")
            edgeListFile.write("Source\tTarget\n")
            nameDict1 = {}
            nameDict2 = {}
            for items in interactionList:
                if items[0] not in nameDict1:
                    nameDict1[items[0]] = items[2]
                else:
                    nameDict1[items[0]] += items[2]
                if items[1] not in nameDict2:
                    nameDict2[items[1]] = items[3]
                else:
                    nameDict2[items[1]] += items[3]
                colour = int((items[4] / 25) + 1)
                if colour > 9:
                    colour = 9
                linksFile.write(
                "{0} {1} {2} {3} {4} {5} color=ylorrd-9-seq-{6}\n".format(
                items[0], nameDict1[items[0]] - items[2], nameDict1[items[0]],
                items[1], nameDict2[items[1]] - items[3], nameDict2[items[1]],
                colour))
                edgeListFile.write("{0}\t{1}\n".format(items[0], items[1]))
            for n in sorted(list(nameDict1.keys())):
                mapFile.write("chr - {0} {0} 0 {1} red\n".format(
                    n, nameDict1[n]))
            for n in sorted(list(nameDict2.keys())):
                mapFile.write("chr - {0} {0} 0 {1} blue\n".format(
                    n, nameDict2[n]))
            linksFile.close()
            mapFile.close()
            edgeListFile.close()
            os.chdir(os.path.dirname(os.path.abspath(__file__))
                + "/" + dirName)
            os.system("circos")
            os.chdir(os.path.dirname(os.path.abspath(__file__)))

    def writeLists(self, fileNameStart, writeMirnaLists=True):
        """
        Wtite list of nodes used by node type for use in cytoscape when
        visualising network.
        """
        i = 1
        for interactionList in self.interactions:
            lstFile = open(fileNameStart + str(i) + ".txt", "w")
            for items in interactionList:
                lstFile.write(str(items[1]) + "\n")
            i += 1
            lstFile.close()

        if writeMirnaLists:
            mirna1File = open(fileNameStart + "mirna1.txt", "w")
            mirna2File = open(fileNameStart + "mirna2.txt", "w")
            for items in self.mirnas:
                mirna1File.write(items[0] + "\n")
                mirna2File.write(items[1] + "\n")
            mirna1File.close()
            mirna2File.close()

    def popInteractions(self):
        return(self.interactions.pop(-1))


class Gene():

    def __init__(self, geneName, entrezId):
        super(Gene, self).__init__()
        self.geneName = geneName
        self.entrezId = entrezId

    def getName(self):
        return self.geneName

    def getEntrez(self):
        return self.entrezId

    def __str__(self):
        return self.geneName


def mirnaFullNames(miRNAs):
    """
    serach MirBase for assesion numbers of mature miRNAs from input list
    """
    mirBase = open("miRNA.dat", "r")
    record = False
    interactions = []
    for line in mirBase:
        if(line[0] == "I" and line[1] == "D"):
            if("hsa" in line and "mir" in line):
                miRNA = line.split()[1]
                number = ''.join([n for n in miRNA.split("-")[2]
                if n in '1234567890'])
                for mi in miRNAs:
                    if number == mi.split("-")[1]:
                        simpleMirna = mi
                        record = True
        elif(line[0] == "F" and line[1] == "T" and record):
            if("/product" in line):
                product = line.split('"')[1]
                interactions.append([simpleMirna, product])
        if(line[0] == "/" and line[1] == "/"):
            record = False
    mirBase.close()
    return interactions


def simplifyMirnaNames(fullmiRNAs):
    """
    Strip miRNA names down to the number eg hsa-miR-1a-1 = miR-1
    For use for with interactions with diseases due to miriad data
    """
    interactions = []
    for mirna in fullmiRNAs:
        interactions.append([mirna.split("-")[1] + "-" + ''.join(
            [n for n in mirna.split("-")[2] if n in '1234567890']), mirna])
    return interactions


def readMiriadData():
    """
    Reads in data from Miriad Database to be processed
    """
    db = open("miRiaDDB.tsv", "r")
    interactions = []
    #read in file
    i = 0
    for line in db:
        if line != "\n":
            elements = line.strip().split("\t")
            if i != 0:
                interactions.append(["miR-" + elements[2], elements[4].replace(
                    ",", "").replace(" ", " ").replace("-", "-")])  # change?
            i += 1
    #finished with file
    db.close()
    return interactions


def readMirTarBase(mirnas="all", removeNonFunctional=True, removeWeak=True,
    useGeneClass=False):
    """
    Find targets for a given set of miRNAs. Select if non-functional and weak
    function interactions are included.

    Return miRNA-Target interactions
    """
    mirTarBase = open("hsa_MTI.csv", "r")
    interactions = []
    if mirnas == "all":
        for line in mirTarBase:
            elements = line.strip().split(",")
            if removeNonFunctional is True and elements[7] == \
            "Non-Functional MTI":
                continue
            if removeWeak is True and elements[7] == "Functional MTI (Weak)":
                continue
            interactions.append([elements[1], elements[3]]) # 3 = name, 4 = entrez
            #print(elements[4])
    elif mirnas == "allEntrez":
        for line in mirTarBase:
            elements = line.strip().split(",")
            if removeNonFunctional is True and elements[7] == \
            "Non-Functional MTI":
                continue
            if removeWeak is True and elements[7] == "Functional MTI (Weak)":
                continue
            interactions.append([elements[1], elements[4]])
    else:
        for line in mirTarBase:
            elements = line.strip().split(",")
            if elements[1] in mirnas:
                if removeNonFunctional is True and elements[7] == \
                "Non-Functional MTI":
                    continue
                if removeWeak is True and elements[7] == \
                "Functional MTI (Weak)":
                    continue
                if useGeneClass:
                    interactions.append([elements[1], Gene(elements[3],
                    elements[4])])
                else:
                    interactions.append([elements[1], elements[3]])

    mirTarBase.close()
    return interactions


def geneNameEntrezDict():
    """
    Create of dictionary of key - gene names, value - entrez number for
    all human genes in MirTarBase
    """
    mirTarBase = open("hsa_MTI.csv", "r")
    geneNameEntrezDict = {}
    for line in mirTarBase:
        elements = line.strip().split(",")
        geneNameEntrezDict[elements[3]] = elements[4]
    mirTarBase.close()
    return geneNameEntrezDict


def readMirnaTSI(lowest, highest=1.0):
    """
    read data for TSI values for miRNAs. Return miRNAs with a between a given
    TSI range
    """
    mirnaTSIFile = open("tsi_Values.csv")
    chosenMirnas = []
    i = 0
    for line in mirnaTSIFile:
        if i != 0:
            values = line.strip().split(",")
            tsiValue = (float(values[4]) + float(values[5])) / 2
            if tsiValue > lowest and tsiValue <= highest:
                #print(line.strip().split(",")[0], ":\t", tsiValue)
                chosenMirnas.append(line.strip().split(",")[0])
        i += 1
    return chosenMirnas


def readTissueMatrix(tsiMiRNAList):
    """
    Find which tissues highly express miRNAs in input list. Return miRNA-Tissue
    interactions
    """
    tissueMatrix = open("data_matrix_quantile.txt", "r")
    tissues = []
    miRNAs = []
    matrix = []
    interactions = []
    i = 0
    for line in tissueMatrix:
        if i == 0:
            #create column names
            for tissue in line.strip().split("\t")[1:]:
                tissues.append(tissue.replace(".", " ").replace(
                    "_", "").strip())  # space or _
        else:
            #create row names and matrix
            miRNAs.append(line.strip().split("\t")[0])
            matrix.append(line.strip().split("\t")[1:])
        i += 1
    i = 0
    tissueMatrix.close()
    for row in matrix:
        floats = [float(value.replace(",", ".")) for value in row]
        j = 0
        for value in floats:
            # check if highest and in high tsi List
            if value > st.median(floats) + st.stdev(floats) and \
            miRNAs[i] in tsiMiRNAList:
                interactions.append([miRNAs[i], tissues[j]])
            j += 1
        i += 1
    return interactions


def reactomePA(inputInteractions, pValue):
    """
    Call reactomePA script and return pathways
    """
    fileInput = "reactomePAInput.txt"
    fileOutput = "rectomePAOutput.txt"
    interactions = []
    interactionDict = {}
    for item in inputInteractions:
        interactionDict.setdefault(item[0], []).append(item[1].getEntrez())

    for key, value in list(interactionDict.items()):
        pathwayFile = open(fileInput, "w")
        for gene in value:
            pathwayFile.write(gene + "\n")
        pathwayFile.close()

        print("loading reactome pathway analysis....")
        command = 'Rscript'
        path2script = 'recatomePathwayScript.R'
        # Build subprocess command
        cmd = [command, path2script] + [fileInput, fileOutput, str(pValue)]
        # check_output will run the command and store to result
        x = subprocess.check_output(cmd, universal_newlines=True)
        print(x)

        rectomePAOutput = open(fileOutput, "r")

        i = 0
        for line in rectomePAOutput:
            if i != 0:
                miRNA = key
                pathway = line.strip().split("\t")[2].replace(" ", "_") \
                .replace("-", "_").replace("\"", "").replace(":", "_") \
                .replace("(", "_").replace(")", "_").replace("/", "_") \
                .replace("+", "").replace(",", "")
                interactions.append([miRNA, pathway])
            i += 1

        rectomePAOutput.close()

    os.remove(fileInput)
    os.remove(fileOutput)

    return interactions


def kegg(inputInteractions):
    from bioservices.kegg import KEGG
    k = KEGG()
    interactions = []
    for items in inputInteractions:
        print(items[1].getName())
        try:
            pathways = k.get_pathway_by_gene(items[1].getName(), "hsa")
            #print(pathways)
            if pathways:
                for key, value in list(pathways.items()):
                    interactions.append([items[0], value])
        except AttributeError:
            print("Gene name error!!!!!!!!!")
    return interactions


