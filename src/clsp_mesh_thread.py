#!/usr/bin/python3.5
# -*-coding: utf-8 -*

from clsp_ga_library import *
from chromosome import *

#--------------------
# Class : MeshThreadsManager
# author : Tafsir GNA
# purpose : Describing a the structure of a manager of slave threads that compute chromosomes for main threads and return the results
#--------------------

class MeshThreadsManager(object):

	"""docstring for MeshThreadsManager"""
	NbMaxPopulation = 0
	crossOverRate = 0
	listMainThreads = 0
	NumberOfMigrants = 0

	def __init__(self, mainThread):

		super(MeshThreadsManager, self).__init__()
		
		# i create and initialize a mesh of thread 
		self.grid = []
		self.fitnessMean = 0
		self.meshRows = 0
		self.meshCols = 0
		self.mainThread = mainThread
		self.rangedThreadsList = []
		self.migrants = []
		self.immigrants = []
		self.locker = threading.Lock()
		self.memory = []

		roundSqrt = math.floor(math.sqrt(MeshThreadsManager.NbMaxPopulation))

		# i make calculations to get the number of rows and columns of my grid
		i = roundSqrt
		while i >= 1:
			if MeshThreadsManager.NbMaxPopulation % i == 0:
				self.meshRows = i
				self.meshCols = int(MeshThreadsManager.NbMaxPopulation / i)
				break
			i -= 1

		#print(self.meshRows, " : ", self.meshCols)

		# i fill the grid with the threads
		i = 0
		prevMeshThread = 0
		while i < self.meshRows:
			row = []
			j = 0
			while j < self.meshCols:

				meshThread = MeshThread(self.mainThread)
				meshThread.queue = copy.deepcopy(self.mainThread.queue)
				meshThread.row = i
				meshThread.column = j

				if i == 0 and j == 0:
					prevMeshThread = meshThread
				else:
					meshThread.flagToWait = prevMeshThread.doneEvent
					prevMeshThread = meshThread

				row.append(meshThread)

				j += 1

			self.grid.append(row)

			i += 1

		# i set the neighbor threads of each thread in the grid 
		i = 0
		while i < self.meshRows:

			j = 0
			while j < self.meshCols:

				# i set the neighbors of each thread
				if i - 1 >= 0:
					(self.grid[i][j]).neighbors.append((self.grid[i-1][j]))
					if j - 1 >= 0:
						(self.grid[i][j]).neighbors.append((self.grid[i-1][j-1]))
						(self.grid[i][j]).neighbors.append((self.grid[i][j-1]))
					if j + 1 <= self.meshCols-1:
						(self.grid[i][j]).neighbors.append((self.grid[i-1][j+1]))
				
				if i + 1 <= self.meshRows-1:
					(self.grid[i][j]).neighbors.append((self.grid[i+1][j]))
					if j + 1 <= self.meshCols-1:
						(self.grid[i][j]).neighbors.append((self.grid[i+1][j+1]))
						(self.grid[i][j]).neighbors.append((self.grid[i][j+1]))
					if j - 1 >= 0:
						(self.grid[i][j]).neighbors.append((self.grid[i+1][j-1]))

				j += 1

			i += 1

	def sendMigrants(self):

		migrants = []
		for i in range(MeshThreadsManager.NumberOfMigrants):
			migrants.append(copy.deepcopy(self.rangedThreadsList[i].chromosome))

		#print("yes " + str(MeshThreadsManager.NumberOfMigrants) + " " + str(migrants))
		for thread in MeshThreadsManager.listMainThreads:
			if thread.getName() != self.mainThread.name:
				thread.meshThreadsManager.receiveImmigrants(migrants)


	def isValidPop(self):
		
		# first condition for a population to be valid

		for row in self.grid:
			for thread in row:
				if thread.chromosome.solution == []:
					return False

		return True

	def receiveImmigrants(self, immigrants):

		self.locker.acquire()
		self.immigrants += copy.deepcopy(immigrants)
		self.locker.release()
	

	def contains(self, obj):

		for thread in self.rangedThreadsList:
			if thread.chromosome.solution == obj:
				return True
		return False


	def putInRank(self, thread):

		popSize = len(self.rangedThreadsList)
		# for inserting the new thread, i range it in the list from the fittest to the least fit
		if (self.rangedThreadsList == []):

			self.rangedThreadsList.append(thread)
			#print(threadQueue)
	
		elif popSize == 1 and (self.rangedThreadsList[0]).chromosome.fitnessValue == 0:

			self.rangedThreadsList.append(thread)
			#print(threadQueue) 

		else:
			# i sort the list of zeroperiods from the most convenient place to the least convenient one
			prevValue = 0
			j = 0
			found = False
			while j < popSize:

				if thread.chromosome.fitnessValue >= prevValue and thread.chromosome.fitnessValue <= (self.rangedThreadsList[j]).chromosome.fitnessValue:
					found = True
					self.rangedThreadsList = self.rangedThreadsList[:j] + [thread] + self.rangedThreadsList[j:]
					break

				prevValue = (self.rangedThreadsList[j]).chromosome.fitnessValue

				j += 1

			if found is False:
				self.rangedThreadsList.append(thread)
		#print("log add : ", (self.grid[len(self.grid)-1]).chromosome)


	def start(self):

		for row in self.grid:
			for thread in row:
				thread.action = 0 # i want to initialize the first population of the algorithm
				thread.start()

		(self.grid[self.meshRows-1][self.meshCols-1]).doneEvent.wait()

		if not self.isValidPop():
			return

		print(" log start : ", self.printPopulation())
		self.putPopInRank()
		#print("Yes : ", self.printRangedThreads())

		# TODO : make sure that the initial population is quite diverse

		'''
		# then, i send the migrants
		self.sendMigrants()

		i = 0
		while True:			

			# i put all the immigrants into the population before starting a new run
			i = len(self.rangedThreadsList) - 1
			for chromosome in self.immigrants:
				(self.rangedThreadsList[i]).chromosome = copy.deepcopy(chromosome)
				i -= 1

			for row in self.grid:
				for thread in row:
					thread.doneEvent.clear()

			for row in self.grid:
				for thread in row:
					thread.action = 2  # i launch the search to find out the best solution to the problem 
					thread.run()

			if self.lacksDiversity():
				#print("Solution is : ", (self.grid[0][0]).chromosome)
				self.memory.append(copy.deepcopy((self.grid[0][0]).chromosome))

				c = copy.deepcopy((self.grid[0][0]).chromosome)
				c.advmutate()

				if c == (self.grid[0][0]).chromosome:
					break
				else:
					(self.grid[0][0]).chromosome = copy.deepcopy(c)
				
				#break

			(self.grid[self.meshRows-1][self.meshCols-1]).doneEvent.wait()

			#print(self.printPopulation())
				
			i += 1

		if self.memory != []:
			c = self.memory[0]
			for chromosome in self.memory:
				if chromosome.fitnessValue < c.fitnessValue:
					c = chromosome

		print("Solution is : ", (self.grid[0][0]).chromosome)



		#(self.grid[self.meshRows-1][self.meshCols-1]).doneEvent.wait()

		#print("Yes : ", self.printPopulation())
		'''
	
	def putPopInRank(self):
		for row in self.grid:
			for thread in row:
				self.putInRank(thread)

	def lacksDiversity(self):

		fitnessValue = (self.grid[0][0]).chromosome.fitnessValue

		for row in self.grid:
			for thread in row:
				if thread.chromosome.fitnessValue != fitnessValue:
					return False
		return True

	def printPopulation(self):

		retVal = ""
		
		if (self.meshRows * self.meshCols) > 0:

			i = 0
			while i < self.meshRows:

				j = 0
				while j < self.meshCols:
					retVal += str(i) + ", " + str(j) + str((self.grid[i][j]).chromosome)
					j += 1

				i += 1

		return retVal

	def printRangedThreads(self):

		retVal = ""
		
		if (self.meshRows * self.meshCols) > 0:

			for thread in self.rangedThreadsList:
				retVal += str(thread.chromosome) + " "

		return retVal


class MeshThread(Thread):

	"""docstring for MeshThread"""

	def __init__(self, mainThread):

		super(MeshThread, self).__init__()

		self._chromosome = Chromosome()
		self.neighbors = []
		#self.threadLocker = threading.Lock()
		self.doneEvent = Event()
		self.flagToWait = 0
		self.row = 0
		self.column = 0
		self.action = -1
		self.mainThread = mainThread
		self.queue = []

	def run(self):

		if self.action == 0: # i want to get the first chromosome of this thread

			self.initSearch()
			#pass

		elif self.action == 2: # i want to mate the current chromosome with its neighbors
		
			if self._chromosome.solution != []:
				c = copy.deepcopy((self.neighbors[0]).chromosome)
				for thread in self.neighbors:
					if thread.chromosome.fitnessValue < c.fitnessValue:
						c = copy.deepcopy(thread.chromosome)
				#print("Log run : ", self._chromosome, " and ", c)

				if c.solution != []:
					#print(" log run : ", self.mainThread.name , self._chromosome, " : ", c)
					self.mateWith(c)

		if self.flagToWait != 0:
			self.flagToWait.wait()
		self.doneEvent.set()
		

	def mateWith(self, chromosome):

		solution3 = []
		solution4 = []

		#print(Population.crossOverRate)
		if (randint(0,100) < (MeshThreadsManager.crossOverRate*100)):

			#print("Log mateWith : ", self.chromosome, " and ", chromosome)
			chromosomes = []
			chromosomes.append(self._chromosome)
			chromosomes.append(chromosome)
			# i retrieve a table that stores the period each item has been manufactered for
			ranks1 = self._chromosome.itemsRank
			ranks2 = chromosome.itemsRank

			ranks3 = []
			ranks4 = []

			randomIndice = randint(1,Chromosome.problem.nbTimes-1)

			#print(" ")

			#print(" 1 - solution1 : ", chromosome1.solution, " ranks1 : ", chromosome1.itemsRank, " solution2 : ", chromosome2.solution, " ranks2 : ", chromosome2.itemsRank)
			#print(" randomIndice : ", randomIndice)
			#print(" ranks1 : ", ranks1, " ranks2 : ", ranks2)

			
			solution3 = self._chromosome.solution[:randomIndice]
			solution4 = chromosome.solution[:randomIndice]

			ranks3 = ranks1[:randomIndice]
			ranks4 = ranks2[:randomIndice]


			solution3 += chromosome.solution[randomIndice:]
			solution4 += self._chromosome.solution[randomIndice:]

			ranks3 += ranks2[randomIndice:]
			ranks4 += ranks1[randomIndice:]

			# Once, the two resulting chromosomes have been formed, i make each of them feasible with regards of the constraints

			#print("Log mateWith 1 : ", self.chromosome, " and "," randomIndice : ", randomIndice)
			#print("Log mateWith 2 : ", self.chromosome, " and "," 2 - solution3 : ", solution3, " ranks3 : ", ranks3, " solution4 : ", solution4, " ranks4 : ", ranks4)

			chromosome3 = Chromosome(solution3, ranks3)
			chromosome3.getFeasible()
			#chromosome3.advmutate()
			chromosomes.append(chromosome3)

			chromosome4 = Chromosome(solution4, ranks4)
			chromosome4.getFeasible()
			#chromosome4.advmutate()
			chromosomes.append(chromosome4)

			#print("Log mateWith 3 : ", self.chromosome, " and ", " 2 - solution3 : ", chromosome3.solution, " ranks3 : ", ranks3, " solution4 : ", chromosome4.solution, " ranks4 : ", ranks4)

			c = chromosomes[0]
			for obj in chromosomes:
				if obj.fitnessValue < c.fitnessValue:
					c = obj

			#print("Log mateWith 2: ", self._chromosome, " and ", chromosome, chromosome3, chromosome4, " and ", c)
			#tempc = copy.deepcopy(self.chromosome)
			#self.chromosom = int(0) #Chromosome(c.solution)
			self._chromosome = copy.deepcopy(c)
			#print("Log mateWith 3: ", " and ", self._chromosome, " and ", c)

			self._chromosome.mutate()
		

	def _set_chromosome(self, new_value):
		if isinstance(new_value, list):
			self._chromosome = Chromosome(new_value)
		elif isinstance(new_value, Chromosome):
			self._chromosome = copy.deepcopy(new_value)

	def _get_chromosome(self):
		return copy.deepcopy(self._chromosome)

	def initSearch(self):

		queueSize = len(self.queue)
		currentNode = copy.deepcopy(self.queue[queueSize-1])
		del self.queue[queueSize-1]
		#print("Queue : ", self.queue)

		i = 0
		while True:

			if currentNode.isLeaf():

				self._chromosome = Chromosome(list(currentNode.solution))
				#self._chromosome.advmutate()
				break

			else:

				print ("current Node : ", currentNode)
				l = currentNode.getChildren()
				print("Children : ", l)
				self.queue += l

			#if i == 10:
			#	break

			#print("inter : ", self.queue)
			queueSize = len(self.queue)
			if queueSize == 0:
				break
	
			currentNode = copy.copy(self.queue[queueSize-1])
			del self.queue[queueSize-1]

			i += 1
	

	def improvChromosome(self):
		self._chromosome.advmutate()

		if self.flagToWait != 0:
			self.flagToWait.wait()
		self.doneEvent.set()
		pass

	# Class' Properties
	chromosome = property(_get_chromosome, _set_chromosome)