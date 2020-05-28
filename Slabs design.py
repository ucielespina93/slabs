
class slab():
	def __init__(self):
		self.__Dws=50 #psf
		self.__L=100 #psf
		self.__lx=240 #in
		self.__cx=15 #in
		self.__ly=240 #in
		self.__cy=15 #in
		self.__fy=60000 #psi
		self.__fc=4000
		self.__wc=150 #psi
		self.__steelDensity=0.283564814814815 #pci

	def __input(self,Dws,L,lx,cx,ly,cy,fy,fc,wc): #INPUT PARAMETER
		self.__Dws=Dws
		self.__L=L
		self.__lx=lx
		self.__ly=ly
		self.__cx=cx
		self.__cy=cy
		self.__fy=fy
		self.__fc=fc
		self.__wc=wc

	#Start Auxiliary functions

	def __linterpolate(self,x,filaX,filaY,matriz): #Function that interpolates the x value in an matrix, that we define the values ​​in x in rowX and y values ​​in row Y
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

	def __rundingup(self,x): #rounding up
		if (x%(x//1) ==0):
			y=x
		else:
			y=x//1+1
		return y

	#End Auxiliary functions

	def __ratioo(self): #1. Ratio of longer to shorter panel dimensions
		self.__lnx=self.__lx-self.__cx
		self.__lny=self.__ly-self.__cy
		self.__lmin=min(self.__lx,self.__ly)
		self.__lmax=max(self.__lx,self.__ly)
		self.__lnmin=min(self.__lnx,self.__lny)
		self.__lnmax=max(self.__lnx,self.__lny)
		self.__ratio=self.__lmax/self.__lmin

	def __minThickness(self): #2. Minimum slab thickness
		if self.__ratio>2:
			self.__hMin=self.__rundingup((self.__lnmin/20)*(0.4+self.__fy/100000))
		else:
			self.__minimumThickness=[[40000,60000,80000],[self.__lnmax/36,self.__lnmax/33,self.__lnmax/30]]
			self.__hMin=self.__rundingup(self.__linterpolate(self.__fy,0,1,self.__minimumThickness))

	def __Deadhm(self): #3. Factored load
		self.__DslabhMin=self.__wc*self.__hMin*(1/12) #1/12 units in-->ft
		self.__DeadhMin=self.__DslabhMin+self.__Dws

	def __quhm(self):
		self.__U1hMin=1.4*self.__DeadhMin
		self.__U2hMin=1.2*self.__DeadhMin+1.6*self.__L
		self.__quhMin=max(self.__U1hMin,self.__U2hMin)*(1/144) #1/144 units psf-->psi
			
	def __LoadCond(self): #Load condition
		if self.__ratio>2:
			if self.__L>3*self.__DeadhMin:
				print("Does not verify the load condition.")
			else:
				print("The load condition verifies.")
		if self.__ratio<=2:
			if self.__L>2*self.__DeadhMin:
				print("Does not verify the load condition.")
			else:
				print("The load condition verifies.")	

	def __Vu(self): #4. Shear Verification
		if self.__ratio>2:
			self.__Vuu=1.15*self.__quhMin*self.__lnmax/2
		else:
			self.__Vuu=self.__quhMin*self.__lmin/2

	def __hforShear(self):
		self.__dd=0.8*self.__hMin
		self.__lamdaS=min(1,(2/(1+(self.__dd/10)))**0.5)
		if self.__ratio>2:
			self.__hShear=self.__rundingup((self.__Vuu)/(2*0.8*0.75*self.__fc**0.5))
		else:
			self.__hShear=self.__rundingup(self.__Vuu/(2*0.8*0.75*self.__lamdaS*self.__fc**0.5))
		self.__h=max(self.__hShear,self.__hMin)
		self.__d=self.__h*0.8

	def __Deadd(self): #Factored load
		self.__Dslab=self.__wc*self.__h*(1/12) #1/12 units in-->ft
		self.__Dead=self.__Dslab+self.__Dws

	def __quu(self):
		self.__U1=1.4*self.__Dead
		self.__U2=1.2*self.__Dead+1.6*self.__L
		self.__qu=max(self.__U1,self.__U2)*(1/144) #1/144 units psf-->psi

	def __factoredMoment(self): #5. Factored Moment
		if self.__ratio > 2:
			self.__wu=self.__qu*self.__lmax
			self.__Mu=[self.__wu*self.__lnmin**2/10,self.__wu*self.__lnmin**2/14]
		else:
			self.__M0x=self.__qu*self.__lx*self.__lny**2/8
			self.__M0y=self.__qu*self.__ly*self.__lnx**2/8
			self.__Mu=[0.70*self.__M0y,0.70*self.__M0y,0.57*self.__M0y,0.70*self.__M0x,0.70*self.__M0x,self.__M0x*0.57]

	def __bb(self):	
		if self.__ratio>2:
			self.__b=[self.__lmax,self.__lmax]
		else:
			self.__b=[self.__ly,self.__ly,self.__ly,self.__lx,self.__lx,self.__lx]

	def __Asss(self): #6. Design moment
		self.__fi=0.9
		def Ass():
			As=[]
			i=0
			for i in range(len(self.__b)):
				Asoli=(0.85*self.__fc*self.__b[i]/self.__fy*(self.__d-((self.__d**2-((2*self.__Mu[i])/(self.__fi*0.85*self.__fc*self.__b[i]))))**0.5))
				AMini=(0.0018*self.__b[i]*self.__h)
				Ai=max(Asoli,AMini)
				As.append(Ai)
			return As
		self.__As=Ass()

	def __Mddd(self):
		def Mdd():
			Md=[]
			i=0
			for i in range(len(self.__b)):
				Md.append(self.__fi*self.__As[i]*self.__fy*(self.__d-0.59*(self.__As[i]*self.__fy/(self.__b[i]*self.__fc))))
			return Md
		self.__Md=Mdd()

	def __quantities(self):
		if self.__ratio >2:
			self.__steelVolumen=(self.__As[0]*self.__lmin+self.__As[1]*self.__lmin)*2
		else:
			self.__steelVolumen=self.__As[0]*self.__lx/3+self.__As[1]*self.__lx/3+self.__As[2]*self.__lx+self.__As[3]*self.__ly/3+self.__As[4]*self.__ly/3+self.__As[5]*self.__ly

		self.__steelQuantities=self.__steelDensity*self.__steelVolumen/(self.__lx*self.__ly)*144 #psf 1728 units ft²-->in²
		self.__concreteQuantities=self.__lx*self.__ly*self.__h*(1/1728)

	def design(self,Dws,L,lx,cx,ly,cy,fy,fc,wc):
		self.__input(Dws,L,lx,cx,ly,cy,fy,fc,wc)
		self.__ratioo()
		self.__minThickness()
		self.__Deadhm()
		self.__quhm()

		if self.__ratio>2:
			print("One-way slab")
		else:
			print("Two-way slab")

		self.__LoadCond()
		self.__Vu()
		self.__hforShear()
		self.__Deadd()
		self.__quu()
		self.__factoredMoment()
		self.__bb()
		self.__Asss()
		self.__Mddd()
		self.__quantities()

		print("The slab thickness is (in): ")
		print(self.__h)
		print("The Steel Quantities are (psf): ")
		print(self.__steelQuantities)
		print("The Concrete Quantities are (ft³): ")
		print(self.__concreteQuantities)
		print("Design finished")
		input(" ")

def main():
	mySlab=slab()
	mySlab.design(50,100,240,8,156,8,60000,4000,150)

if __name__ == "__main__":
    main()

