
import sys
import maya.OpenMayaMPx as omMPx
import maya.OpenMaya as om
import math
import numpy as np
import random
import png
from scipy.spatial import Voronoi as vor
from collections import defaultdict
from region import Region

''' how to rounding
points = [(-7.66584,0,1.287069),(-7.067755,0,-0.577034),(-5.871584,0,-4.305242),(-0.435031,0,5.852464),(5.794102,0,-8.565758),(6.457351,0,5.83384),(9.227738,0,1.930198),(10.612932,0,-0.0216233)]
import maya.cmds as cmds

cmds.file( newFile=True, force=True )

cmds.unloadPlugin( 'terrainNode.py' )
cmds.loadPlugin( 'terrainNode.py' )

cmds.createNode( 'TerrainNode', name='terrainNode1' )

cmds.createNode("transform", name="terrain1")
cmds.createNode( 'mesh', name='terrainShape1', parent="terrain1" )
cmds.sets("terrainShape1", add="initialShadingGroup")

cmds.curve( p=points )

# Connect the attributes.
cmds.connectAttr( 'curveShape1.local', 'terrainNode1.inputCurve' )

cmds.connectAttr( 'terrainNode1.outputMesh', 'terrainShape1.inMesh' )
'''




# Plug-in information:
kPluginNodeName = 'TerrainNode'           # The name of the node.
kPluginNodeId = om.MTypeId( 0xBEEF4 ) # A unique ID associated to this node type.

##########################################################
# Plug-in 
##########################################################
class TerrainNode(omMPx.MPxNode):
    # Static variables which will later be replaced by the node's attributes.
    inputCurveAttribute = om.MObject()
    outputMeshAttribute = om.MObject()

    class Voronoi:
        def __init__(self, vorObjnp ):
            self.vertices = vorObjnp.vertices.tolist()

            #remove all infinites values
            self.regions = self.removeNegatives( vorObjnp.regions )

    
        def removeNegatives( self, listPoints ):
            rList = []
            for x in listPoints:
                nonNegative = True
                for y in x:
                    if y == -1:
                        nonNegative = False
                if nonNegative is True:
                    #not empty
                    if x:
                        rList.append(x)
            return rList

    def __init__(self):
        ''' Constructor. '''
        # (!) Make sure you call the base class's constructor.
        omMPx.MPxNode.__init__(self)
        self.regionsDictObj={}
        self.adjVertices = None

        #vertices to be used to generate terrain mesh
        self.vertices3d = None
        
    
    def compute(self, pPlug, pDataBlock):
        ''' Here, we will create a voxelized version of the input mesh. '''
        
        if( pPlug == TerrainNode.outputMeshAttribute ):
            
            # Get our custom input node attributes and values.
            # print TerrainNode.inputCurveAttribute
            inputMeshHandle = pDataBlock.inputValue( TerrainNode.inputCurveAttribute )
            curveObj = inputMeshHandle.asNurbsCurve()
            inCurve = om.MFnNurbsCurve(curveObj)

            self.cPoints = self.computePointsCurve( inCurve )
            # print self.cPoints

            #add limit to map to avoid infinity values
            rPoints = self.computeBoundary(self.cPoints)

            vorObj = self.computeVoronoi( rPoints )

            #create a dict with all the region Objects
            self.createRegionObjects( vorObj, self.cPoints )

            #vertex X is shared by N faces
            self.computeAdjacentVertices( )

            self.computeRegionBoundary( )

            self.computeRiverAltitude( )

            ########
            #Image generation
            ########
            #the maximum point in the terrain is
            #ALWAYS the last value in rPoints
            maxPoint = rPoints[-1][0]
            self.generateImage ( maxPoint, inCurve )


            ########
            #Node output
            ########
            outputHandle = pDataBlock.outputValue( TerrainNode.outputMeshAttribute )

            dataCreator = om.MFnMeshData()
            newOutputData = dataCreator.create()

            self.generateTerrain( vorObj, newOutputData )

            outputHandle.setMObject(newOutputData)
            pDataBlock.setClean(pPlug)

            
        else:
            return om.kUnknownParameter

    def generateTerrain( self, vorObj, outData ):

        self.vertices3d = []
        # print self.adjVertices

        #insert zeros to Y values
        for i in range( 0, len( vorObj.vertices ) ):
            tv = vorObj.vertices[i]
            tv.insert( 1, 0.0 )
            self.vertices3d.append( tv )

        rObj = self.regionsDictObj
        for index, regionsI in self.adjVertices.items():
            maxY = 0
            # print 'v: ' + str(index)
            for regionI in regionsI:
                rPointY = rObj[regionI].pointY
                # print 'r: ' + str(regionI) + ' alt: ' + str(rPointY)
                if rPointY > maxY:
                    maxY = rPointY
                    # print 'max: ' + str(maxY)
            self.vertices3d[index][1] = maxY


        # for region in self.regionsDictObj.values():
        #     for vertexIndex in region.regionI:
        #         # self.vertices3d[vertexIndex][1] = region.pointY
        #         if region.nextEdge is not None:
        #             # if vertexIndex in region.nextEdge:
        #                 #sum Y value to vertex
        #             yDist = region.pointY
        #             print yDist
        #             self.vertices3d[vertexIndex][1] = yDist

        vtx = []
        for v in self.vertices3d:
            vtx.append (om.MFloatPoint(v[0], v[1], v[2]))

        numVertices = len(vtx)

        points = om.MFloatPointArray()
        points.setLength(numVertices)
        for i in range(0, numVertices):
            points.set(vtx[i], i)


        vorFaces = []
        pointCount = 0
        faceCount = 0

        numPolConnects = 0

        #get the number of polygons on the terrain
        numPolCounts = len(self.regionsDictObj)
        polCounts = om.MIntArray()
        polCounts.setLength(numPolCounts)

        polCountI = 0
        for region in self.regionsDictObj.values():
            numVerticesRegion = len(region.regionI)
            polCounts.set(numVerticesRegion, polCountI)
            numPolConnects = numPolConnects + numVerticesRegion
            polCountI = polCountI + 1

        # print numPolConnects
        polConnects = om.MIntArray()
        polConnects.setLength(numPolConnects)


        faceConI = 0
        for region in self.regionsDictObj.values():
            for vertexIndex in region.regionI:
                polConnects.set(vertexIndex, faceConI)
                faceConI = faceConI + 1


        #checklist to create a mesh:
        ## numVertices - ok
        ## numPolygons - numPolCounts
        ## vertexArray - points
        ## polygonCounts - polCounts
        ## polygonConnects - polConnects
        ## 
        meshFS = om.MFnMesh()
        newMesh = meshFS.create(numVertices, numPolCounts, points, polCounts, polConnects, outData)

        return newMesh

    def computeVoronoi( self, points ):
        #transform array to numpy
        npPoints = np.array( points )
        # print npPoints
        #compute voronoi
        npVorObj = vor( npPoints )
        return self.Voronoi( npVorObj )

    def createRegionObjects ( self, vorObj, cPoints ):
        regionVertices = self.indexToVertex( vorObj.vertices, vorObj.regions )
        # print len(regionVertices)
        for region, regionIndex in zip(regionVertices, vorObj.regions):
            for i in range(0, len(cPoints)):
                if self.pointInRegion ( cPoints[i], region ):
                    self.regionsDictObj[i] = Region( i, cPoints[i], regionIndex )
                    continue

    def computeRegionBoundary( self ):

        for i in range( 0, len(self.regionsDictObj) ):
            rI = self.regionsDictObj[i].regionI
            for j in range ( i+1, len(self.regionsDictObj) ):
                rJ = self.regionsDictObj[j].regionI
                # print 'comparing ' + str(rI) + ' with ' + str(rJ)
                edge = self.intersect(rI, rJ)
                if  edge is not None:
                    if self.regionsDictObj[i].id == self.regionsDictObj[j].id-1:
                        self.regionsDictObj[i].setBoundary(edge)
                        continue
    

    def intersect( self, a, b):
        result = list(set(a) & set(b))
        if not result:
            return None
        else:
            return result


    def computeRiverAltitude( self, angle=0.1 ):
        #default angle is 1 degree

        #generate a random value from 0 to angle
        rAngle = []
        for i in range( 0, len(self.regionsDictObj) ):
            rAngle.append( random.uniform(0, angle) )
        # print 'angle ' + str(rAngle)
        tY = 0 #tY is the total sum of y 
        eY = 0
        #for from second point to last
        for i in range( 1, len(self.regionsDictObj) ):
            dist = self.calcDist(self.regionsDictObj[i].point, self.regionsDictObj[i-1].point)
            # print 'dist ' + str(dist)
            eY = (dist/2)*math.tan(rAngle[i]) + tY
            tY += dist*math.tan(rAngle[i])
            self.regionsDictObj[i].setEPAltitude( tY )
            self.regionsDictObj[i-1].setNextEdgeAltitude( eY )
            # print eY

    def calcDist( self, pA, pB ):
        result = ( pA[0] - pB[0] )*( pA[0] - pB[0] ) + ( pA[1] - pB[1] )*( pA[1] - pB[1] )
        result = math.sqrt(result)
        return result

    #receives Region indeces
    #returns Region with vertices values
    def indexToVertex( self, vorVertices, vorRegions ):
        regionVertices = []
        for region in vorRegions:
            # print region
            rVertices = []
            for vIndex in region:
                rVertices.append( vorVertices[vIndex] )
            regionVertices.append( rVertices )
        return regionVertices


        #check if point is inside region
    def pointInRegion( self, point, region ):
        poly = region
        x = point[0]
        y = point[1]

        n = len(poly)
        inside = False

        p1x,p1y = poly[0]
        for i in range(n+1):
            p2x,p2y = poly[i % n]
            if y > min(p1y,p2y):
                if y <= max(p1y,p2y):
                    if x <= max(p1x,p2x):
                        if p1y != p2y:
                            xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xints:
                            inside = not inside
            p1x,p1y = p2x,p2y

        return inside

    def computeAdjacentVertices( self ):
        self.adjVertices = defaultdict( list )
        for key, region in self.regionsDictObj.items():
            # print region.regionI
            for vertex in region.regionI:
                if key not in self.adjVertices[vertex]:
                    self.adjVertices[vertex].append (key)
                # print self.adjVertices[vertex]

    def computePointsCurve( self, curve, divisions = 1 ):
        cLen = curve.length()
        nPars =  curve.findParamFromLength(cLen)
        # print 'num of parameters' + str(nPars)
        #equal distance between each parameter
        #math.ceil solve problem with rounding
        lenPars = (nPars / math.ceil(nPars)) / divisions
        # print 'len btw par' + str(lenPars)
        sumPars = 0
        cPoints = []
        #p is a point object
        p = om.MPoint()

        while sumPars <= nPars:
            curve.getPointAtParam(sumPars, p)
            #ignoring y because voronoi is 2D and y SHOULD be always zero
            pointL = [p.x, p.z]
            #append each point point to the array
            cPoints.append(pointL)
            # print sumPars
            sumPars = sumPars + lenPars
        # print cPoints
        return cPoints

    def computeBoundary( self, cPoints ):
        maxPoint = self.computeMaxValue( cPoints )
        boundary = maxPoint * 2
        # print boundary

        #append maximum point to the point list
        cPoints.append( [-boundary, -boundary] )
        cPoints.append( [boundary, -boundary] )
        cPoints.append( [-boundary, boundary] )
        cPoints.append( [boundary, boundary] )
        return cPoints
    
    #points is a list of two tuples
    def computeMaxValue( self, points ):
        # unzip points to get the max value
        xP, zP = zip( *points )
        maxX = max( map( abs, xP ) )
        maxZ = max( map( abs, zP ) )
        maxPoint = max( maxX, maxZ )

        return maxPoint

##########################################################
# Image Generation.
##########################################################
    def generateImage( self, maxPoint, curve ):
        #width, heigth
        imgRes = [ 100, 100 ]
        points = self.computePointsCurve( curve, 100 )

        # image for PNG output, old implementation
        image = self.createImgArrayForPNG( points, maxPoint, imgRes )
        # print image
        imgObj = png.from_array(image, 'RGB')
        imgObj.save( 'test.png')
        # for i in imgObj.rows
        #     print i
        # for line in image:
        #     print line


        #TODO
        '''
        imgObj = om.MImage()
        imgObj.create(imgRes[0], imgRes[1], 3, OpenMaya.MImage.kFloat )

        #image for NODE output
        self.populateImgArrayForNode( imgObj.floatPixels(), points, maxPoint, imgRes )
        '''

    def populateImgArrayForNode( self, points, maxPoint, imgRes ):
        #this function should be used to output an image to the crowd sim
        # image = self.fillZeros ( imgRes )
        index = 0
        for point in points:
            pointInImage = self.newRange (maxPoint, point, imgRes)
            #set 1 where the curve pass
            pixel = (imgRes[0] * pointInImage[0]) + (pointInImage[1]*3)
            image[ pixel + 0 ] = 255
            image[ pixel + 1 ] = 255
            image[ pixel + 2 ] = 255

            # image[ pointInImage[0]   ][ pointInImage[1]*3 + 0  ] = 255
            # image[ pointInImage[0]   ][ pointInImage[1]*3 + 1  ] = 255
            # image[ pointInImage[0]   ][ pointInImage[1]*3 + 2  ] = 255

            # # #set 1 to the NEIGHBOURS
            # image[ pointInImage[0]+1 ][ pointInImage[1]*3 + 0  ] = 255
            # image[ pointInImage[0]+1 ][ pointInImage[1]*3 + 1  ] = 255
            # image[ pointInImage[0]+1 ][ pointInImage[1]*3 + 2  ] = 255

            # image[ pointInImage[0]   ][ (pointInImage[1]+1)*3 + 0 ] = 255
            # image[ pointInImage[0]   ][ (pointInImage[1]+1)*3 + 1 ] = 255
            # image[ pointInImage[0]   ][ (pointInImage[1]+1)*3 + 2 ] = 255

            # image[ pointInImage[0]+1 ][ pointInImage[1]+1 ] = white

            # image[ pointInImage[0]-1 ][ pointInImage[1]   ] = white

            # image[ pointInImage[0]   ][ pointInImage[1]-1 ] = white

            # image[ pointInImage[0]-1 ][ pointInImage[1]-1 ] = white

        # #initial point should be RED
        # pointInImage = self.newRange (maxPoint, points[0], imgRes)
        # image[ pointInImage[0]   ][ pointInImage[1]*3 + 0  ] = 255
        # image[ pointInImage[0]   ][ pointInImage[1]*3 + 1  ] = 0
        # image[ pointInImage[0]   ][ pointInImage[1]*3 + 2  ] = 0

        # ##goal point is GREEN
        # pointInImage = self.newRange (maxPoint, points[-1], imgRes)
        # image[ pointInImage[0]   ][ pointInImage[1]*3 + 0  ] = 0
        # image[ pointInImage[0]   ][ pointInImage[1]*3 + 1  ] = 255
        # image[ pointInImage[0]   ][ pointInImage[1]*3 + 2  ] = 0
        return image


    #create array to be used with PNG Library
    def createImgArrayForPNG( self, points, maxPoint, imgRes ):
        image = self.fillZeros ( imgRes )
        # print image
        for point in points:
            pointInImage = self.newRange (maxPoint, point, imgRes)
            #set 1 where the curve pass
            # print image[ pointInImage[0]   ][ pointInImage[1]*3 + 0  ]
            image[ pointInImage[0]   ][ pointInImage[1]*3 + 0  ] = 255
            image[ pointInImage[0]   ][ pointInImage[1]*3 + 1  ] = 255
            image[ pointInImage[0]   ][ pointInImage[1]*3 + 2  ] = 255

            # #set 1 to the NEIGHBOURS
            image[ pointInImage[0]+1 ][ pointInImage[1]*3 + 0  ] = 255
            image[ pointInImage[0]+1 ][ pointInImage[1]*3 + 1  ] = 255
            image[ pointInImage[0]+1 ][ pointInImage[1]*3 + 2  ] = 255

            image[ pointInImage[0]   ][ (pointInImage[1]+1)*3 + 0 ] = 255
            image[ pointInImage[0]   ][ (pointInImage[1]+1)*3 + 1 ] = 255
            image[ pointInImage[0]   ][ (pointInImage[1]+1)*3 + 2 ] = 255

        #initial point should be RED
        pointInImage = self.newRange (maxPoint, points[0], imgRes)
        image[ pointInImage[0]   ][ pointInImage[1]*3 + 0  ] = 255
        image[ pointInImage[0]   ][ pointInImage[1]*3 + 1  ] = 0
        image[ pointInImage[0]   ][ pointInImage[1]*3 + 2  ] = 0

        ##goal point is GREEN
        pointInImage = self.newRange (maxPoint, points[-1], imgRes)
        image[ pointInImage[0]   ][ pointInImage[1]*3 + 0  ] = 0
        image[ pointInImage[0]   ][ pointInImage[1]*3 + 1  ] = 255
        image[ pointInImage[0]   ][ pointInImage[1]*3 + 2  ] = 0

        return image

    def newRange ( self, maxValue, value, imgRes ):
        oldRange = 2*maxValue
        # oldRange = (OldMax - OldMin)
        newValue = [ 0, 0 ]
        newValue[0] = int((((value[1] + maxValue) * imgRes[0]) / oldRange))
        newValue[1] = int((((value[0] + maxValue) * imgRes[1]) / oldRange))
        return newValue


    def fillZeros( self, imgRes, PNG=True ):

        if PNG:
            #list row flat pixel construction
            image = []
            for i in range( 0, imgRes[0] ):
                line = []
                for j in range( 0, imgRes[1] ):
                    r = 0
                    line.append( r )
                    g = 0
                    line.append( g )
                    b = 0
                    line.append( b )

                image.append( line )

        else:
            #flat row flat pixel construction
            image = []
            for i in range( 0, imgRes[0] ):
                for j in range( 0, imgRes[1] ):
                    r = 0
                    image.append( r )
                    g = 0
                    image.append( g )
                    b = 0
                    image.append( b )

        return image
        
##########################################################
# Plug-in initialization.
##########################################################
def nodeCreator():
    ''' Creates an instance of our node class and delivers it to Maya as a pointer. '''
    return omMPx.asMPxPtr( TerrainNode() )

def nodeInitializer():
    ''' Defines the input and output attributes as static variables in our plug-in class. '''
    # The following MFnNumericAttribute function set will allow us to create our attributes.
    nAttr = om.MFnNumericAttribute()
    
    # This one allows us to create our input and output mesh attributes.
    tAttr = om.MFnTypedAttribute()
    
    #==================================
    # INPUT NODE ATTRIBUTE(S)
    #==================================
    # We will need an input mesh attribute.
    # TerrainNode.inputCurveAttribute = nAttr.create( 'inputCurve', 'ic', om.MFnNumericData.k3Double )
    # nAttr.setWritable( True )
    # nAttr.setReadable( False )
    # nAttr.setStorable( False )
    # nAttr.setHidden( False )
    # nAttr.setArray( True )
    # TerrainNode.addAttribute( TerrainNode.inputCurveAttribute )
    TerrainNode.inputCurveAttribute = tAttr.create( 'inputCurve', 'ic', om.MFnData.kNurbsCurve )
    tAttr.setWritable( True )
    tAttr.setReadable( False )
    tAttr.setStorable( False )
    tAttr.setHidden( False )
    TerrainNode.addAttribute( TerrainNode.inputCurveAttribute )

    #==================================
    # OUTPUT NODE ATTRIBUTE(S)
    #==================================
    TerrainNode.outputMeshAttribute = tAttr.create( 'outputMesh', 'om',
                                                                 om.MFnData.kMesh )
    tAttr.setWritable( False )
    tAttr.setReadable( True )
    tAttr.setStorable( False )
    tAttr.setHidden( False )
    TerrainNode.addAttribute( TerrainNode.outputMeshAttribute )
    
    #==================================
    # NODE ATTRIBUTE DEPENDENCIES
    #==================================
    # If any of the inputs change, the output mesh will be recomputed.
    TerrainNode.attributeAffects( TerrainNode.inputCurveAttribute, TerrainNode.outputMeshAttribute )
    
    
def initializePlugin( mobject ):
    ''' Initialize the plug-in '''
    mplugin = omMPx.MFnPlugin( mobject )
    try:
        mplugin.registerNode( kPluginNodeName, kPluginNodeId, nodeCreator, nodeInitializer )
    except:
        sys.stderr.write( 'Failed to register node: ' + kPluginNodeName )
        raise
    print 'loaded'
    
def uninitializePlugin( mobject ):
    ''' Uninitializes the plug-in '''
    mplugin = omMPx.MFnPlugin( mobject )
    try:
        mplugin.deregisterNode( kPluginNodeId )
    except:
        sys.stderr.write( 'Failed to deregister node: ' + kPluginNodeName )
        raise