import math

class slab():

	"""This class contains the methods and properties common to slabs, from it inherit the two types of loas, One-Way and Two-Way
	"""

	def __init__(self):

		"""This method is the constructor. Contains the common inputs on the slabs, loads and characteristics of the materials

		"""

		self._deadLoadWithoutSlabWeight=50
		self._liveLoad=100

		self._strengthReductionFactor=0.90

		self._specifiedYieldStrengthReinforcement=60000
		self._specifiedCompressiveStrengthConcrete=4000
		self._concreteDensity=150
		self._steelDensity=0.283564814814815

	def geometricParameters(self,lengthSpanX,columnDimensionX,lengthSpanY,columnDimensionY):

		"""This method contains the particular inputs on the slabs, characteristics of the materials
		Arguments:
			
			lengthSpanX = length of span parallel to the x-axis, measured center-to-center of supports, in
			lengthSpanY = length of span parallel to the y-axis, measured center-to-center of supports, in
			columnDimensionX = column dimensions parallel to the x-axis, in
			columnDimensionY= column dimensions parallel to the y-axis, in


		Returns:

			lengthSpanX = length of span parallel to the x-axis, measured center-to-center of supports, in
			lengthSpanY = length of span parallel to the y-axis, measured center-to-center of supports, in
			columnDimensionX = column dimensions parallel to the x-axis, in
			columnDimensionY= column dimensions parallel to the y-axis, in


			lengthClearSpanX = length of clear span measured face-to-face of supports parallel to the x-axis, in
			lengthClearSpanY = length of clear span measured face-to-face of supports parallel to the y-axis, in
			lengthSpanMinimum = length of span measured face-to-face of supports minimum, in
			lengthSpanMaximum = length of span measured face-to-face of supports maximum, in
			lengthClearSpanMinimum = length of clear span measured face-to-face of supports minimum, in
			lengthClearSpanMaximum = length of clear span measured face-to-face of supports maximum, in
			ratioLengthSpan = ratio of longer to shorter panel dimensions

		"""

		self._lengthSpanX=lengthSpanX
		self._columnDimensionX=columnDimensionX
		self._lengthSpanY=lengthSpanY
		self._columnDimensionY=columnDimensionY

		self._lengthClearSpanX=self._lengthSpanX-self._columnDimensionX
		self._lengthClearSpanY=self._lengthSpanY-self._columnDimensionY
		self._lengthSpanMinimum=min(self._lengthSpanX,self._lengthSpanY)
		self._lengthSpanMaximum=max(self._lengthSpanX,self._lengthSpanY)
		self._lengthClearSpanMinimum=min(self._lengthClearSpanX,self._lengthClearSpanY)
		self._lengthClearSpanMaximum=max(self._lengthClearSpanX,self._lengthClearSpanY)
		self.ratioLengthSpan=self._lengthSpanMaximum/self._lengthSpanMinimum

	def _linealInterpolate(self,x,filaX,filaY,matriz): 

		"""	This method interpolates the x value in an matrix, that we define the values ​​in x in rowX and y values ​​in row Y
		Arguments:

			x = value in X to interpolate
			filaX = row number where array x is located
			filaY = row number where array y is located
			matriz = matrix containing arrays x and y

		Returns:
			
			y = interpolated value

		"""

		vectorX=[]
		vectorY=[]
		for i in range(len(matriz[filaX])):
			if x>=matriz[filaX][i] and x<matriz[filaX][i+1]:
				vectorX.append(matriz[filaX][i])
				vectorX.append(matriz[filaX][i+1])
				vectorY.append(matriz[filaY][i])
				vectorY.append(matriz[filaY][i+1])
		y=(vectorY[0])+(((vectorY[1])-(vectorY[0]))/((vectorX[1])-(vectorX[0])))*(x-(vectorX[0]))
		return y

	def _loadForMinimumStiffnessThickness(self):

		"""This method calculates the factored load for minimum stiffness thickness

		Returns:

			deadLoadForMinimumStiffnessThickness = dead load for minimum stiffness thickness, psi
			factoredLoadForMinimumStiffnessThickness = factored load for minimum stiffness thickness, psi

		"""

		#In deadLoadSlab multiply by (1/12) to convert units in-->ft
		deadLoadSlab=self._concreteDensity*self._minimumThickness*(1/12)
		self._deadLoadForMinimumStiffnessThickness=deadLoadSlab+self._deadLoadWithoutSlabWeight

		loadCombination1=1.4*self._deadLoadForMinimumStiffnessThickness
		loadCombination2=1.2*self._deadLoadForMinimumStiffnessThickness+1.6*self._liveLoad

		#In factoredLoadForMinimumStiffnessThickness multiply by (1/144) to convert units psf-->psi
		self._factoredLoadForMinimumStiffnessThickness=max(loadCombination1,loadCombination2)*(1/144) 
			
	def _factoredLoadForDefinitiveThickness(self):

		"""This method calculates the factored load for definitive thickness

		Returns:

			factoredLoad = factored load for definitive thickness, psi

		"""

		#In deadLoadSlab multiply by (1/12) to convert units in-->ft
		deadLoadSlab=self._concreteDensity*self.slabThickness*(1/12)
		deadLoad=deadLoadSlab+self._deadLoadWithoutSlabWeight

		loadCombination1=1.4*deadLoad
		loadCombination2=1.2*deadLoad+1.6*self._liveLoad

		#In factoredLoad multiply by (1/144) to convert units psf-->psi
		self._factoredLoad=max(loadCombination1,loadCombination2)*(1/144)

	def _slabReinforcementArea(self):

		"""This method calculates the reinforcement area to verify with the minimum requirements and the factored moment

		Returns:
			
			reinforcementArea = reinforcement area for minimum requirements and the factored moment, in²
			strengthReductionFactor = strength reduction factor for tension controller
		"""

		self._strengthReductionFactor=0.9
		self._reinforcementArea=[]
		i=0
		for i in range(len(self._lengthAppliesMoments)):
			reinforcementForMomentRequirements=(0.85*self._specifiedCompressiveStrengthConcrete*self._lengthAppliesMoments[i]/self._specifiedYieldStrengthReinforcement*(self._distanceToTensionReinforcement-((self._distanceToTensionReinforcement**2-((2*self._factoredMoments[i])/(self._strengthReductionFactor*0.85*self._specifiedCompressiveStrengthConcrete*self._lengthAppliesMoments[i]))))**0.5))
			reinforcementForMinimumRequirements=(0.0018*self._lengthAppliesMoments[i]*self.slabThickness)
			self._reinforcementArea.append(max(reinforcementForMomentRequirements,reinforcementForMinimumRequirements))
		
	def _slabMomentDesign(self):

		"""This method calculates the moment design of the slab

		Returns:
			
			momentDesign = moment design of the slab, lbf*in

		"""

		self._momentDesign=[]
		i=0
		for i in range(len(self._lengthAppliesMoments)):
			self._momentDesign.append(self._strengthReductionFactor*self._reinforcementArea[i]*self._specifiedYieldStrengthReinforcement*(self._distanceToTensionReinforcement-0.59*(self._reinforcementArea[i]*self._specifiedYieldStrengthReinforcement/(self._lengthAppliesMoments[i]*self._specifiedCompressiveStrengthConcrete))))
		
class oneWay(slab):

	"""This class inherits the common methods from the Class slabs and adds the particular ones from One-Way slabs. 
		The class calculates the quantities of steel and concrete in a slab
	"""
	
	def design(self,lengthSpanX,columnDimensionX,lengthSpanY,columnDimensionY):

		"""This method specifies the order to execute the methods to calculate the quantities of steel and concrete of a slab
		Arguments:
			
			lengthSpanX = length of span parallel to the x-axis, measured center-to-center of supports, in
			lengthSpanY = length of span parallel to the y-axis, measured center-to-center of supports, in
			columnDimensionX = column dimensions parallel to the x-axis, in
			columnDimensionY= column dimensions parallel to the y-axis, in

		Returns:
			
			slabThickness = thickness of the slab, in
			steelQuantities = quantities of steel, psf
			concreteQuantities = quantities of concrete, ft³

		"""

		self.geometricParameters(lengthSpanX,columnDimensionX,lengthSpanY,columnDimensionY)
		self._minimumStiffnessThickness()
		self._loadForMinimumStiffnessThickness()
		self._verificationOfLoadCondition()
		self._factoredShear()
		self._tickness()
		self._factoredLoadForDefinitiveThickness()
		self._factoredMomentsForDefinitiveThickness()
		self._slabReinforcementArea()
		self._slabMomentDesign()
		self._quantities()

	def _minimumStiffnessThickness(self):

		"""This method calculates the minimum thickness of the slab for stiffness requirements

		Returns:

			minimumThicknes = minimum thickness of the slab for stiffness requirements, in

		"""

		self._minimumThickness=math.ceil((self._lengthClearSpanMinimum/20)*(0.4+self._specifiedYieldStrengthReinforcement/100000))


	def _verificationOfLoadCondition(self):

		"""This method checks that the load condition is verified for the chosen calculation procedure used

		Returns:
			loadCondition = boolean value indicating that the load condition is verified

		"""

		if self._liveLoad>3*self._deadLoadForMinimumStiffnessThickness:
			self._loadCondition=True
		else:
			self._loadCondition=False

	def _factoredShear(self):

		"""This method calculates the factored shear for minimum stiffness thickness

		Returns:
			shear = factored shear for minimum stiffness thickness, lbf/in

		"""

		self._shear=1.15*self._factoredLoadForMinimumStiffnessThickness*self._lengthClearSpanMaximum/2


	def _tickness(self):

		"""This method calculates the the thickness of the slab, which is the greatest of the shear and stiffness requirements.

		Returns:
			slabThickness = thickness of the slab, which is the greatest of the shear and stiffness requirements, in
			distanceToTensionReinforcement = distance from extreme compression fiber to centroid of longitudinal tension reinforcement, in

		"""

		shearThickness=math.ceil((self._shear)/(2*0.8*0.75*self._specifiedCompressiveStrengthConcrete**0.5))
		self.slabThickness=max(shearThickness,self._minimumThickness)
		self._distanceToTensionReinforcement=self.slabThickness*0.8

	def _factoredMomentsForDefinitiveThickness(self):

		"""This method calculates the factored moments for definitive thickness, and the length where it is applies 

		Returns:

			factoredMoments = factored moments for definitive thickness, lbf*in
			lengthAppliesMoments = length where it is applies the factored moments, in

		"""
	
		factoredLoadPerUnitLength=self._factoredLoad*self._lengthSpanMaximum
		self._factoredMoments=[factoredLoadPerUnitLength*self._lengthClearSpanMinimum**2/10,factoredLoadPerUnitLength*self._lengthClearSpanMinimum**2/14]
		self._lengthAppliesMoments=[self._lengthSpanMaximum,self._lengthSpanMaximum]

	def _quantities(self):

		"""This method calculates the quantities of concrete and steel for the slab

		Returns:
			
			steelQuantities = quantities of steel, psf
			concreteQuantities = quantities of concrete, ft³

		"""

		steelVolumen=(self._reinforcementArea[0]*self._lengthSpanMinimum+self._reinforcementArea[1]*self._lengthSpanMinimum)*2
		
		#In steelQuantities multiply by (144) to convert units ft²-->in²
		self.steelQuantities=self._steelDensity*steelVolumen/(self._lengthSpanX*self._lengthSpanY)*144

		#In concreteQuantities multiply by (1/728) to convert units in³-->ft³
		self.concreteQuantities=self._lengthSpanX*self._lengthSpanY*self.slabThickness*(1/1728)

class twoWay(slab):

	"""This class inherits the common methods from the Class slabs and adds the particular ones from Two-Way slabs. 
		The class calculates the quantities of steel and concrete in a slab
	"""

	def design(self,lengthSpanX,columnDimensionX,lengthSpanY,columnDimensionY):

		"""This method specifies the order to execute the methods to calculate the quantities of steel and concrete of a slab
		Arguments:
			
			lengthSpanX = length of span parallel to the x-axis, measured center-to-center of supports, in
			lengthSpanY = length of span parallel to the y-axis, measured center-to-center of supports, in
			columnDimensionX = column dimensions parallel to the x-axis, in
			columnDimensionY= column dimensions parallel to the y-axis, in

		Returns:
			
			slabThickness = thickness of the slab, in
			steelQuantities = quantities of steel, psf
			concreteQuantities = quantities of concrete, ft³

		"""

		self.geometricParameters(lengthSpanX,columnDimensionX,lengthSpanY,columnDimensionY)
		self._minimumStiffnessThickness()
		self._loadForMinimumStiffnessThickness()
		self._verificationOfLoadCondition()
		self._factoredShear()
		self._tickness()
		self._factoredLoadForDefinitiveThickness()
		self._factoredMomentsForDefinitiveThickness()
		self._slabReinforcementArea()
		self._slabMomentDesign()
		self._quantities()

	def _minimumStiffnessThickness(self):

		"""This method calculates the minimum thickness of the slab for stiffness requirements

		Returns:

			minimumThicknes = minimum thickness of the slab for stiffness requirements, in

		"""

		minimumThicknessBySpecifiedYieldStrengthReinforcement=[[40000,60000,80000],[self._lengthClearSpanMaximum/36,self._lengthClearSpanMaximum/33,self._lengthClearSpanMaximum/30]]
		self._minimumThickness=math.ceil(self._linealInterpolate(self._specifiedYieldStrengthReinforcement,0,1,minimumThicknessBySpecifiedYieldStrengthReinforcement))


	def _verificationOfLoadCondition(self):

		"""This method checks that the load condition is verified for the chosen calculation procedure used

		Returns:
			loadCondition = boolean value indicating that the load condition is verified

		"""

		if self._liveLoad>2*self._deadLoadForMinimumStiffnessThickness:
			self._loadCondition=True
		else:
			self._loadCondition=False

	def _factoredShear(self):

		"""This method calculates the factored shear for minimum stiffness thickness

		Returns:
			shear = factored shear for minimum stiffness thickness, lbf/in

		"""

		self._shear=self._factoredLoadForMinimumStiffnessThickness*self._lengthSpanMinimum/2

	def _tickness(self):

		"""This method calculates the the thickness of the slab, which is the greatest of the shear and stiffness requirements.

		Returns:
			slabThickness = thickness of the slab, which is the greatest of the shear and stiffness requirements, in
			distanceToTensionReinforcement = distance from extreme compression fiber to centroid of longitudinal tension reinforcement, in

		"""

		sizeEffectFactorForShear=min(1,(2/(1+(0.8*self._minimumThickness/10)))**0.5)
		shearThickness=math.ceil(self._shear/(2*0.8*0.75*sizeEffectFactorForShear*self._specifiedCompressiveStrengthConcrete**0.5))
		self.slabThickness=max(shearThickness,self._minimumThickness)
		self._distanceToTensionReinforcement=self.slabThickness*0.8


	def _factoredMomentsForDefinitiveThickness(self):

		"""This method calculates the factored moments for definitive thickness, and the length where it is applies 

		Returns:

			factoredMoments = factored moments for definitive thickness, lbf*in
			lengthAppliesMoments = length where it is applies the factored moments, in

		"""

		momentToDistributeInX=self._factoredLoad*self._lengthSpanX*self._lengthClearSpanY**2/8
		momentToDistributeInY=self._factoredLoad*self._lengthSpanY*self._lengthClearSpanX**2/8
		self._factoredMoments=[0.70*momentToDistributeInY,0.70*momentToDistributeInY,0.57*momentToDistributeInY,0.70*momentToDistributeInX,0.70*momentToDistributeInX,momentToDistributeInX*0.57]
		self._lengthAppliesMoments=[self._lengthSpanY,self._lengthSpanY,self._lengthSpanY,self._lengthSpanX,self._lengthSpanX,self._lengthSpanX]

	def _quantities(self):

		"""This method calculates the quantities of concrete and steel for the slab

		Returns:

			steelQuantities = quantities of steel, psf
			concreteQuantities = quantities of concrete, ft³

		"""

		steelVolumen=self._reinforcementArea[0]*self._lengthSpanX/3+self._reinforcementArea[1]*self._lengthSpanX/3+self._reinforcementArea[2]*self._lengthSpanX+self._reinforcementArea[3]*self._lengthSpanY/3+self._reinforcementArea[4]*self._lengthSpanY/3+self._reinforcementArea[5]*self._lengthSpanY
		
		#In steelQuantities multiply by (144) to convert units ft²-->in²
		self.steelQuantities=self._steelDensity*steelVolumen/(self._lengthSpanX*self._lengthSpanY)*144

		#In concreteQuantities multiply by (1/728) to convert units in³-->ft³
		self.concreteQuantities=self._lengthSpanX*self._lengthSpanY*self.slabThickness*(1/1728)


if __name__ == "__main__":

	"""This main function contains an example
	"""

	mySlab=slab()
	mySlab.geometricParameters(315,15,315,15)

	if mySlab.ratioLengthSpan > 2:
		mySlab=oneWay()
		print("One-way slab")
	else:
		mySlab=twoWay()
		print("Two-way slab")

	mySlab.design(315,15,315,15)


	print("The slab thickness is ",mySlab.slabThickness," in")
	print("The Steel Quantities are ",round(mySlab.steelQuantities,2)," psf")
	print("The Concrete Quantities are ",round(mySlab.concreteQuantities,2)," ft³")
	print("Design finished")
	input(" ")

