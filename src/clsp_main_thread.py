#!/usr/bin/python3.5
# -*-coding: utf-8 -*

from clsp_mesh_thread import *

#--------------------
# Class : ClspThread
# author : Tafsir GNA
# purpose : Describing the structure of an instance to the algorithm
#--------------------

class ClspThread(Thread):

	NbGenToStop = 0

	# Builder 
	def __init__(self, threadId, queue):

		Thread.__init__(self)
		self.threadId = threadId
		self.name = "Thread - " + str(threadId)
		self.queue = queue
		self.thread_memory = []
		self.migrants = []
		self.meshThreadsManager = MeshThreadsManager(self)

	def run(self):

		self.meshThreadsManager.start()

		'''
		self.population.getImproved()

		print (self.name, " ", "Initial Population : ", self.population)
		
		self.population.thread_memory = self.thread_memory
		self.population.startPopData = []
		self.population.startPopData.append(copy.deepcopy(self.population.chromosomes))
		self.population.startPopData.append(copy.deepcopy(self.population.listFitnessData))

		# i send the best chromosomes of the population to its neighbors
		self.sendMigrants()

		nbGenB4Stop = 0
		storedFitnessMean = self.population.fitnessMean

		# After the initial population has been created, i launch the search process
		i = 1
		indiceMigration = 0

		while True:

			if self.migrants != []:
				for chromosome in self.migrants:
					self.population.replace(chromosome)
				self.population.listFitnessData = []
				self.population.getFitnessData()

			print ("Thread : ", self.name, "Population : ", i, self.population.chromosomes, " and ", self.population.listFitnessData, " Population's memory : ", Population.ga_memory, self.population.chromosomes[0] in Population.ga_memory)
			#print("population : ", self.population.fitnessMean)

			if len(self.population.chromosomes) <= 1:
				#print("the thread ", self.name, " has exited!")
				break


			#if i == 10:
			#	break 

			population = Population()
			retVal = population.initialize(indiceMigration, self.population)

			if population.fitnessMean >= self.population.fitnessMean:
				nbGenB4Stop += 1
				if nbGenB4Stop == ClspThread.NbGenToStop and self.thread_memory != []:
					#print("stop")
					break
			else:
				nbGenB4Stop == 0

			if retVal == 1: # this signals it's time for migration
				#print("Yo")
				self.sendMigrants()

			self.population = population

			# i wait for two conditions to be true before exiting
			# first is :
			if len(self.thread_memory) > 0:
				print("len")
				break

			
			population.renew()
			population.listFitnessData = []
			population.getFitnessData()

			#print ("Current population : ", self.population.chromosomes, " and ", self.population.listFitnessData)

			i += 1
		'''	

	'''
	def sendMigrants(self):

		if self.population.chromosomes != []:

			chromosomes = []
			i = 0
			while i < ClspThread.NumberOfMigrants:
				chromosomes.append(copy.deepcopy(self.population.chromosomes[i]))
				i += 1

			for thread in ClspThread.listMainThreads:
				if thread.getName() != self.name:
					thread.receiveMigrants(copy.deepcopy(chromosomes))

	def receiveMigrants(self, chromosomes):
		self.migrants += chromosomes
	'''