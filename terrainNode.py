
import sys
import maya.OpenMayaMPx as omMPx
import maya.OpenMaya as om
import math
import numpy as np
import random
from scipy.spatial import Voronoi as vor
from collections import defaultdict
from region import Region

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

            self.vertices3d = None
    
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
        self.adjVertices = defaultdict( list )
        
    
    def compute(self, pPlug, pDataBlock):
        ''' Here, we will create a voxelized version of the input mesh. '''
        
        if( pPlug == TerrainNode.outputMeshAttribute ):
            
            # Get our custom input node attributes and values.
            # print TerrainNode.inputCurveAttribute
            inputMeshHandle = pDataBlock.inputValue( TerrainNode.inputCurveAttribute )
            curveObj = inputMeshHandle.asNurbsCurve()
            inCurve = om.MFnNurbsCurve(curveObj)

            cLen = inCurve.length()
            nPars =  inCurve.findParamFromLength(cLen)
            # print 'num of parameters' + str(nPars)

            self.cPoints = self.computePointsCurve( nPars, inCurve )
            print self.cPoints

            #add limit to map to avoid infinity values
            rPoints = self.computeBoundary(self.cPoints)

            vorObj = self.computeVoronoi( rPoints )

            #create a dict with all the region Objects
            self.createRegionObjects( vorObj, self.cPoints )

            #vertex X is shared by N faces
            self.computeAdjacentVertices( )

            self.computeRegionBoundary( )

            self.computeRiverAltitude( )

            
        else:
            return om.kUnknownParameter


    def computeVoronoi( self, points ):
        #transform array to numpy
        npPoints = np.array( points )
        print npPoints
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
            print eY

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
        for key, region in self.regionsDictObj.items():
            for vertex in region.regionI:
                self.adjVertices[vertex].append (key)

    def computePointsCurve( self, nPars, curve ):
        #equal distance between each parameter
        #math.ceil solve problem with rounding
        lenPars = nPars / math.ceil(nPars)
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

            sumPars = sumPars + lenPars

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