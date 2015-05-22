
import sys
import maya.OpenMayaMPx as omMPx
import maya.OpenMaya as om
import math

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
    
    
    def __init__(self):
        ''' Constructor. '''
        # (!) Make sure you call the base class's constructor.
        omMPx.MPxNode.__init__(self)
        
    
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
            print 'num of parameters' + str(nPars)

            #equal distance between each parameter
            #math.ceil solve problem with rounding
            lenPars = nPars / math.ceil(nPars)
            print 'len btw par' + str(lenPars)
            sumPars = 0
            while sumPars <= nPars:
                point = om.MPoint()
                inCurve.getPointAtParam(sumPars, point)
                print point.x
                print point.y
                print point.z
                sumPars = sumPars + lenPars
                print 'sum' + str(sumPars)
            
        else:
            return om.kUnknownParameter
    
        
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