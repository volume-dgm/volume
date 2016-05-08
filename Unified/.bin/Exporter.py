import os
import salome
import sys

import SMESH
import SALOMEDS
from salome.smesh import smeshBuilder
from salome.kernel import termcolor
# print dir(smeshBuilder)


def resize_list(l, size):
    l += [0] * (max(0, size - len(l)))


def find_selected_meshes():
    """

    :rtype : list
    """
    meshes = list()
    smesh = smeshBuilder.New(salome.myStudy)
    nrSelected = salome.sg.SelectedCount()  # total number of selected items

    foundMesh = False
    for i in range(nrSelected):
        selected = salome.sg.getSelected(i)
        selobjID = salome.myStudy.FindObjectID(selected)
        selobj = selobjID.GetObject()
        print "selected class: "
        print selobj.__class__
        if selobj.__class__ == SMESH._objref_SMESH_Mesh or selobj.__class__ == salome.smesh.smeshBuilder.meshProxy:
            # mName = selobjID.GetName().replace(" ","_")
            foundMesh = True
            mesh = smesh.Mesh(selobj)
            meshes.append(mesh)

    if not foundMesh:
        print "Script Error: you have to select a mesh"

    return meshes

outFolder = 'OutMeshes'
if not os.path.exists(outFolder):
    os.makedirs(outFolder)


meshes = find_selected_meshes()

if len(meshes) > 0:
    detectorMeshes = list()
    geomMeshes = list()

    for mesh in meshes:
        if "[D]" in mesh.GetName():
            detectorMeshes.append(mesh)
        else:
            geomMeshes.append(mesh)

    salomeCellType = None
    salomeBoundaryType = None

    for mesh in geomMeshes:
        outFileName = outFolder + '/Salome_%(meshName)s' % {'meshName': mesh.GetName()}
        outFile = open(outFileName + '.mesh', 'w')
        paramsFile = open(outFileName + '.params', 'wb')

        dimension = mesh.MeshDimension()

        if dimension == 2:
            salomeCellType = SMESH.FACE
            salomeBoundaryType = SMESH.EDGE
        elif dimension == 3:
            salomeCellType = SMESH.VOLUME
            salomeBoundaryType = SMESH.FACE

        print "Mesh: ", mesh.GetName()
        print "Dimension: ", dimension
        nodeIds = mesh.GetNodesId()
        outFile.write('%(nodesCount)d\n' % {'nodesCount': len(nodeIds)})

        for nodeId in nodeIds:
            xyz = mesh.GetNodeXYZ(nodeId)
            for coordIndex in range(dimension):
                outFile.write('%(coord)f ' % {'coord': xyz[coordIndex]})
            outFile.write('\n')

        cells = mesh.GetElementsByType(salomeCellType)
        outFile.write('\n')
        outFile.write('%(indicesCount)d\n' % {'indicesCount': len(cells) * (dimension + 1)})

        for cell in cells:
            cellNodes = mesh.GetElemNodes(cell)
            for nodeIndex in cellNodes:
                outFile.write('%(node)d ' % {'node': nodeIndex - 1})
            outFile.write('\n')

        maxContactTypeNumber = 0
        edgesOfContactTypeCount = []

        maxBoundaryTypeNumber = 0
        edgesOfBoundaryTypeCount = []

        groupTypeNumber = []
        groupTypes = []  # 0 is contact 1 is boundary

        totalContactEdges = 0
        totalBoundaryEdges = 0

        for groupIndex in range(len(mesh.GetGroups())):
            group = mesh.GetGroups()[groupIndex]
            groupTypeNumber += [-1]
            groupTypes += [-1]

            if group.GetType() == salomeBoundaryType:
                name = group.GetName()
                number = int(-1)
                try:
                    number = int(name.split('[')[1].split(']')[0])
                except:
                    pass

                if number >= 0:
                    # print "contact %d" % groupIndex
                    if number > maxContactTypeNumber:
                        maxContactTypeNumber = number
                    groupTypeNumber[groupIndex] = number
                    groupTypes[groupIndex] = 0  # contact
                    resize_list(edgesOfContactTypeCount, number + 1)
                    edgesOfContactTypeCount[number] += len(group.GetIDs())
                    totalContactEdges += len(group.GetIDs())

                number = int(-1)
                try:
                    number = int(name.split('{')[1].split('}')[0])
                except:
                    pass

                if number >= 0:
                    # print "boundary %d" % groupIndex
                    if number > maxBoundaryTypeNumber:
                        maxBoundaryTypeNumber = number
                    groupTypeNumber[groupIndex] = number
                    groupTypes[groupIndex] = 1  # boundary
                    resize_list(edgesOfBoundaryTypeCount, number + 1)
                    edgesOfBoundaryTypeCount[number] += len(group.GetIDs())
                    totalBoundaryEdges += len(group.GetIDs())

        outFile.write('\n%(totalContactEdges)d\n' % {'totalContactEdges': totalContactEdges})
        for typeIndex in range(len(edgesOfContactTypeCount)):
            for groupIndex in range(len(mesh.GetGroups())):
                group = mesh.GetGroups()[groupIndex]
                if groupTypeNumber[groupIndex] == typeIndex and groupTypes[groupIndex] == 0:
                    for edgeId in group.GetIDs():
                        edgeNodes = mesh.GetElemNodes(edgeId)
                        # double because left side == right side and will be duplicated in mesh builder
                        for num in range(2):
                            for nodeIndex in range(dimension):
                                outFile.write('%(node)d ' % {'node': edgeNodes[nodeIndex] - 1})

        print "Contact types: %d " % len(edgesOfContactTypeCount)
        outFile.write('\n')
        outFile.write('%d ' % len(edgesOfContactTypeCount))
        outFile.write('\n')
        for typeIndex in range(len(edgesOfContactTypeCount)):
            outFile.write('%d ' % edgesOfContactTypeCount[typeIndex])

        outFile.write('\n%(totalBoundaryEdges)d\n' % {'totalBoundaryEdges': totalBoundaryEdges})

        for typeIndex in range(len(edgesOfBoundaryTypeCount)):
            for groupIndex in range(len(mesh.GetGroups())):
                group = mesh.GetGroups()[groupIndex]
                if groupTypeNumber[groupIndex] == typeIndex and groupTypes[groupIndex] == 1:
                    for edgeId in group.GetIDs():
                        edgeNodes = mesh.GetElemNodes(edgeId)
                        for nodeIndex in range(dimension):
                            outFile.write('%(node)d ' % {'node': edgeNodes[nodeIndex] - 1})

        print "Boundary types: %d " % len(edgesOfBoundaryTypeCount)
        outFile.write('\n')
        outFile.write('%d ' % len(edgesOfBoundaryTypeCount))
        outFile.write('\n')
        for typeIndex in range(len(edgesOfBoundaryTypeCount)):
            outFile.write('%d ' % edgesOfBoundaryTypeCount[typeIndex])
            print "Boundary edges of type %(type)d: %(count)d " % {'type' : typeIndex, 'count' : edgesOfBoundaryTypeCount[typeIndex]}


        totalDetectorsCount = 0
        for detectorMesh in detectorMeshes:
            totalDetectorsCount += len(detectorMesh.GetNodesId())

        print "Detectors count: %d\n" % totalDetectorsCount
        outFile.write('%d ' % totalDetectorsCount)

        for detectorMesh in detectorMeshes:
            for nodeId in detectorMesh.GetNodesId():
                xyz = detectorMesh.GetNodeXYZ(nodeId)
                for coordIndex in range(dimension):
                    outFile.write('%(coord)f ' % {'coord': xyz[coordIndex]})
                outFile.write('\n')

        elements = mesh.GetElementsId()
        elementSubmeshIndices = len(elements) * [0]  # we'll have edge elements as dummies in this array as well but they won't be used

        print "Total cells count: %d" % len(cells)

        for groupIndex in range(len(mesh.GetGroups())):
            group = mesh.GetGroups()[groupIndex]
            groupType = -1

            if group.GetType() == salomeCellType:
                name = group.GetName()
                try:
                    groupType = int(name.split('<')[1].split('>')[0])
                except:
                    pass
                print "Cell group" + name + "is index %d" % groupType

            if groupType >= 0:
                for elementIndex in group.GetIDs():
                    if mesh.GetElementType(elementIndex, True) == salomeCellType:
                        elementSubmeshIndices[elementIndex - 1] = groupType

        cellSubmeshIndices = []
        for elementIndex in range(len(elementSubmeshIndices)):
            if mesh.GetElementType(elementIndex + 1, True) == salomeCellType:
                cellSubmeshIndices.append(elementSubmeshIndices[elementIndex])

        print "Cell params count %d" % len(cellSubmeshIndices)

        paramsFileByteArray = bytearray(cellSubmeshIndices)

        print "Cell params array size %d" % len(paramsFileByteArray)
        paramsFile.write(paramsFileByteArray)

        outFile.close()
        paramsFile.close()

        print "Exporting done"