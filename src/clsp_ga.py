# -*-coding: utf-8 -*

from clsp_main_thread import *
from clsp_mesh_thread import *

#--------------------
# Class : GeneticAlgorithm
# author : Tafsir GNA
# purpose : Describing the structure of the kind of genetic algorithm used in the program
#--------------------

class GeneticAlgorithm:

	#	Class' variables
	NbMaxPopulation = 25
	mutationRate = 0.05
	crossOverRate = 0.80
	FITNESS_PADDING = 1
	NumberOfMigrants = 1
	MigrationRate = 0 # this variable holds the number of generations needed before a migration occurs during the search
	nbMainThreads = 3
	NbGenToStop = 4

	# Builder
	def __init__(self, inst):

		self.instance = inst
		self.hashTable = {} #hashTable is a dictionnary
		self.ga_memory = []
		self.listMainThreads = [] # i initialize a list that's intended to contain all the main threads of this genetic algorithm program


		# i impart some parameters to the chromosome class and population class
		MeshThreadsManager.NbMaxPopulation = GeneticAlgorithm.NbMaxPopulation
		MeshThreadsManager.crossOverRate = GeneticAlgorithm.crossOverRate
		MeshThreadsManager.FITNESS_PADDING = GeneticAlgorithm.FITNESS_PADDING

		Chromosome.mutationRate = GeneticAlgorithm.mutationRate
		Chromosome.problem = self.instance
		Chromosome.hashTable = self.hashTable

		MeshThreadsManager.listMainThreads = self.listMainThreads
		ClspThread.NbGenToStop = GeneticAlgorithm.NbGenToStop
		MeshThreadsManager.NumberOfMigrants = GeneticAlgorithm.NumberOfMigrants


	#--------------------
	# function : start
	# Class : GeneticAlgorithm
	# purpose : Initializing le population of chromosomes to be processed during the first iteration
	#--------------------

	def start(self):

		# In order to create this new population, i use the deep first search(DFS) to create some potentially good chromosomes

		# i pick the item, i will start the scheduling with
		item = 1#randint(1, Chromosome.problem.nbItems)
		period = Chromosome.problem.deadlineDemandPeriods[item-1][0][0]
		nbperiods = (period + 1) * Chromosome.problem.nbCapacities
		rootPerThread = math.ceil(nbperiods / GeneticAlgorithm.nbMainThreads)

		print (" log start : ", item, " : ", period)

		# i initialize each thread and put it into the corresponding list 
		i = period
		threadQueue = []
		threadCounter = 0

		#Node.evaluate([[0,0], [1,1], [2,2], [1,2], [3,1], [0,0], [0,0], [3,2], [2,2], [3,2]])
		#Node.evaluate([[1, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]])
		
		while i >= 0:
		
			# i initialize the node from which the search of each thread will start

			for j in range(Chromosome.problem.nbCapacities):

				root = Node()
				#print("root 1 : ", root)
				root.currentItem = item
				root.currentPeriod = 1	

				solution = []
				for k in range(Chromosome.problem.nbTimes * 2):
					solution.append([0,0])

				solution[i + j * Chromosome.problem.nbTimes][0] = item
				solution[i + j * Chromosome.problem.nbTimes][1] += 1
				root.solution = solution
				root.currentQuantity -= 1

				print("root : ", root)
				self.putInThreadQueue(root, threadQueue)


			if  len(threadQueue) >= rootPerThread or i == 0:

				# i initialize the thread and put it into a list of threads created for this purpose
				clspThread = ClspThread(threadCounter, list(reversed(threadQueue)))
				self.listMainThreads.append(clspThread)

				threadQueue = []
				threadCounter += 1

				#if threadCounter == 1:
				#	break
			i -= 1
		
		# want to make sure the parent process will wait for the child thread before exiting
		#for thread in self.listMainThreads:
		#	thread.start()
		#	thread.join()

		self.printResults()
	

	def putInThreadQueue(self, root, threadQueue):

		if (threadQueue == []):

			threadQueue.append(copy.deepcopy(root))
			#print(threadQueue)
			
		elif len(threadQueue) == 1 and (threadQueue[0]).fitnessValue == 0:

			threadQueue.append(copy.deepcopy(root))
			#print(threadQueue)

		else:
			
			# i sort the list of zeroperiods from the most convenient place to the least convenient one
			size = len(threadQueue)
			prevValue = 0
			j = 0
			found = False
			while j < size:
				
				if root.fitnessValue >= prevValue and root.fitnessValue <= (threadQueue[j]).fitnessValue:
					threadQueue = threadQueue[:j] + [copy.deepcopy(root)] + threadQueue[j:]
					found = True
					break
				prevValue = (threadQueue[j]).fitnessValue
				j += 1

			if found is False:
				threadQueue.append(copy.deepcopy(root))

	#--------------------
	# function : printResults
	# Class : GeneticAlgorithm
	# purpose : This prints the results of the search once the search has ended
	#--------------------
	def printResults(self):

		#print(self.ga_memory)
		if len(self.ga_memory) != 0:

			self.bestSolution = self.ga_memory[0] # variable that holds the best solution found so far in the search

			for chromosome in self.ga_memory:
				if chromosome.fitnessValue < self.bestSolution.fitnessValue:
					self.bestSolution = chromosome

			print("The best solution found so far is : ", self.bestSolution.solution)
			#print(self.population)
			print("The fitness of this solution is : ", self.bestSolution.fitnessValue)

		else:
			
			print("No solutions have been found so far. Please Try again.")
	