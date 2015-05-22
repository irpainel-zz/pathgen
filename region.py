class Region:
	
	def __init__(self, idEP, point, regionI ):
		self.id = idEP
		#edit point coordinates
		self.point = point
		self.pointY = 0
		self.regionI = regionI
		self.nextEdge = None
		self.nextEdgeY = 0

	#set the edge that is boundary with the next EP ( next river Point )
	def setBoundary( self, edge ):
		self.nextEdge = edge

	def setEPAltitude( self, y ):
		self.pointY = y

	def setNextEdgeAltitude( self, y ):
		self.nextEdgeY = y
