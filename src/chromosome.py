#!/usr/bin/python3.5
# -*-coding: utf-8 -*

from clsp_ga_library import *

class Chromosome(object):

	mutationRate = 0
	problem = 0
	hashTable = {}

	# Builder 
	def __init__(self, solution = [], itemsRank = []):

		# Variables
		self._solution = []
		self._itemsRank = []
		self._fitnessValue = 0
		self.itemsRankFlag = False
		self._hashSolution = ""
		self.manufactItemsPeriods = getManufactPeriodsGrid(Chromosome.problem.nbItems, Chromosome.problem.deadlineDemandPeriods) #Chromosome.problem.manufactItemsPeriods 

		if solution != []:

			self._solution = list(solution)
			#if itemsRank == []:
			#	self._get_itemsRanks()

			if self.isFeasible():
				self._get_hashSolution()
				self._get_fitnessValue()

		if itemsRank != []:
			self._itemsRank = list(itemsRank)

	# Getters

	def _get_fitnessValue(self):

		if self.isFeasible() is False:
			self._fitnessValue = 0
			return self._fitnessValue

		if self._fitnessValue == 0:

			if self.hashSolution not in Chromosome.hashTable:
				
				self._fitnessValue = Node.evaluate(self._solution)

				hashTableData = []
				hashTableData.append(self._fitnessValue)
				#hashTableData.append(self.itemsRank)
				Chromosome.hashTable[self.hashSolution] = list(hashTableData)

			else:
				hashTableData = Chromosome.hashTable[self.hashSolution]
				self._fitnessValue = hashTableData[0]

		#if self._fitnessValue < 375:
		#	print(self._solution, self._fitnessValue)

		return self._fitnessValue

	def _get_solution(self):
		return self._solution

	def _get_itemsRanks(self):

		#if self.itemsRankFlag is False:

		if self.hashSolution not in Chromosome.hashTable:

			if self.itemsRankFlag is False:

				self._itemsRank = []
				gridCounters = [1] * Chromosome.problem.nbItems
				#print("grid : ", gridCounters)

				i = 0 
				while i < Chromosome.problem.nbTimes:

					if self._solution[i] != 0:

						item = self._solution[i]
						counter = gridCounters[item-1]
						self._itemsRank.append(counter)

						# then, i increment the counter of this item
						gridCounters[item-1] = (counter+1)

					else:
						self._itemsRank.append(0)

					i+=1

				self.itemsRankFlag = True

				hashTableData = []
				hashTableData.append(self.fitnessValue)
				hashTableData.append(self._itemsRank)
				Chromosome.hashTable[self.hashSolution] = list(hashTableData)

		else:

			hashTableData = Chromosome.hashTable[self.hashSolution]
			self._itemsRank = hashTableData[1]

		return self._itemsRank

	def _get_hashSolution(self):

		if self._hashSolution == "":

			i = 0
			while i < Chromosome.problem.nbTimes:
				self._hashSolution += str(self._solution[i])
				i+=1
			
		return self._hashSolution


	# Setters

	def _set_hashSolution(self, new_value):
		self._hashSolution = new_value

	def _set_itemsRanks(self, new_value):
		self._itemsRank = new_value

	def _set_solution(self, new_solution):
		self._solution = new_solution

	def _set_fitnessValue(self, new_value):
		self._fitnessValue = new_value

	def __repr__(self):
		return " {} : {} ".format(self._solution,self.fitnessValue)

	def __eq__(self, chromosome):
		return self._solution == chromosome.solution

	def __ne__(self, chromosome):
		return self._solution != chromosome.solution
	
	# Genetic operators and other function

	def isFeasible(self):

		# i check first that there's not shortage or backlogging

		manufactGrid = copy.deepcopy(Chromosome.problem.deadlineDemandPeriods)

		for row in manufactGrid:
			for cell in row:
				cell[1] = 0

		#print(" yop ", manufactGrid)

		#listCounters = [1] * Chromosome.problem.nbItems

		for i in range(Chromosome.problem.nbTimes):

			for j in range(Chromosome.problem.nbCapacities):

				item = self._solution[i + j * Chromosome.problem.nbTimes][0]
				quantity = self._solution[i + j * Chromosome.problem.nbTimes][1]
				#print (" item : ", item)
				if item != 0:

					# i want to get the interval of periods to which the current item belongs to
					mini = 0
					inc = 0
					for cell in Chromosome.problem.deadlineDemandPeriods[item-1]:
						if i >= mini and i <= cell[0]: 
							break
						mini = cell[0]
						inc += 1

					#print (" item : ", item, inc)
					manufactGrid[item-1][inc][1] += quantity

		#print(manufactGrid, " : ", self._solution)
		#print (" hop : ", Chromosome.problem.deadlineDemandPeriods)

		for i in range(Chromosome.problem.nbItems):

			quantity = 0
			#sum1 = 0
			#sum2 = 0
			for j in range(len(Chromosome.problem.deadlineDemandPeriods[i])):
				expected = Chromosome.problem.deadlineDemandPeriods[i][j][1]
				quantity += manufactGrid[i][j][1]

				#sum1 += expected
				#sum2 += manufact

				#print(i+1, " : ", j+1, " : ", expected, " : ", quantity)
				if  quantity < expected:
					#print ("Falses returned ! ")
					return False
				else:
					quantity -= expected
					
		#print ("True returned ! ")
		return True

	#--------------------
	# function : mutate
	# Class : Chromosome
	# purpose : Applying mutation to a given chromosome and returning the resulting one
	#--------------------

	def mutate(self):

		#print("M Start : ", self._solution)
		saved_solution = self._solution
		saved_fitnessValue = self._fitnessValue

		if (randint(0,100) < (Chromosome.mutationRate*100)): # then the chromsome has been selected for mutation 

			mutated = False

			while mutated is False:

				randomIndice = randint(0,(Chromosome.problem.nbTimes-1))

				item1 = self._solution[randomIndice]

				# i make sure that the randomIndice variable never corresponds to a zero indice
				while item1 == 0:
					randomIndice = randint(0,(Chromosome.problem.nbTimes-1))
					# i get the item corresponding the gene to be flipped
					item1 = self._solution[randomIndice]

				#print(" randomIndice : ", randomIndice)

				visitedItems = []

				i = randomIndice-1
				itemsRank = self.itemsRank
				while i >= 0:
					if self._solution[i] != item1 and self._solution[i] != 0:

						item2 = self._solution[i]
						#print(" item2 : ", item2)

						if item2 not in visitedItems:

							visitedItems.append(item2)

							item2DemandPeriod = Chromosome.problem.deadlineDemandPeriods[item2-1][itemsRank[i]-1]

							if item2DemandPeriod >= randomIndice:
								#print(i, randomIndice)
								solution = switchGenes(self._solution, randomIndice, i)
								c = Chromosome(solution)
								self._solution = c.solution
								self._fitnessValue = c.fitnessValue
								self._hashSolution = c.hashSolution
								self._itemsRank = c.itemsRank
								self.manufactItemsPeriods = list(c.manufactItemsPeriods)
								mutated = True
								break
					i-=1

				# if the first approach doesn't work ,then i apply another one in order to leave the current chromosome actually mutated

				if mutated is False:

					visitedItems = []

					i = randomIndice + 1
					while i < Chromosome.problem.nbTimes:

						if self._solution[i] != item1 and self._solution[i] != 0:

							item2 = self._solution[i]

							if item2 not in visitedItems:

								item1DemandPeriod = Chromosome.problem.deadlineDemandPeriods[item1-1][itemsRank[randomIndice]-1]

								if item1DemandPeriod >= i:
									#print(i, randomIndice)
									solution = switchGenes(self._solution, randomIndice, i)
									c = Chromosome(solution)
									self._solution = c.solution
									self._fitnessValue = c.fitnessValue
									self._hashSolution = c.hashSolution
									self._itemsRank = c.itemsRank
									self.manufactItemsPeriods = list(c.manufactItemsPeriods)
									mutated = True
									break					

						i += 1
			

		#print("F Start : ", self._solution)

	# TODO Revamp advmutate function as function mutate is

	def advmutate(self):

		solution1 = list(self._solution)
		itemsRank1 = self.itemsRank

		i = 0
		while i < Chromosome.problem.nbTimes:

			if solution1[i] != 0:

				item1 = solution1[i]

				item2 = 1

				while item2 <= Chromosome.problem.nbItems :
					
					if item2 != item1:

						item2DemandPeriods = Chromosome.problem.deadlineDemandPeriods[item2-1]

						#print(" i : ", i," item2 : ", item2, " item2DemandPeriods : ", item2DemandPeriods)
						j = i
						solution2 = []
						while j >= 0:
							if solution1[j] == item2:
								#print(" item's rank value : ", itemsRank[j], " j : ", j)
								if item2DemandPeriods[itemsRank1[j]-1] >= i:
									solution2 = switchGenes(solution1, j, i)
									itemsRank2 = switchGenes(itemsRank1, j, i)									
									break
							j-=1

						# i check if the resulting chromosome would have a better fitness than the current's fitness
						if solution2 != []:

							c = Chromosome(solution2)
							if c.fitnessValue < self.fitnessValue:
								self._solution = c.solution
								self._itemsRank = c.itemsRank
								self._fitnessValue = c.fitnessValue
								self._hashSolution = c.hashSolution

					item2 += 1

			i+=1

	
	def getFeasible(self):

		#print(" In Chromosome 1 : ", self._solution)
		#print(self._solution)

		if self.isFeasible() is False:

			#print(" grid : ", grid)
			copy_solution = list(self._solution)

			# i make sure that the number of goods producted isn't superior to the number expected
			i = 0
			while i < Chromosome.problem.nbTimes:

				if self._solution[i] != 0:

					item = self._solution[i]
					#print(" ok : ", self._solution, self._itemsRank)
					#print(" item picked : ", item)
					rank = self._itemsRank[i]
					#print(i, item-1, rank-1, self.manufactItemsPeriods)
					value = self.manufactItemsPeriods[item-1][rank-1]

					if value == -1:
						itemDemandPeriods = self.manufactItemsPeriods[item-1]
						itemDemandPeriods[rank-1] = i
						#print(" It isn't yet in the tab")
						#print(" == -1 ", item, i, rank)

					else:

						#print(" != -1 ", item, i, rank)
						#print(" It is already in the tab")
						cost1 = Chromosome.getCostof(value, item, rank, copy_solution, i)
						cost2 = Chromosome.getCostof(i, item, rank, copy_solution, value)

						#print(" cost 1 : ", cost1, " cost2 : ", cost2)
						if cost2 < cost1 :
							itemDemandPeriods = self.manufactItemsPeriods[item-1]
							itemDemandPeriods[rank-1] = i

							#print(" cost2 < cost1 : ", value, item)
							self._solution[value] = 0

						else:
							self._solution[i] = 0
				i+=1

			#print(" in middle getFeasible : ", self._solution, ", ", self._itemsRank)
			#print()
			# i make sure that the number of items producted isn't inferior to the number expected
			i = 0
			while i < Chromosome.problem.nbItems:

				j = 0
				nbmanufactItemsPeriods = len(self.manufactItemsPeriods[i])
				while j < nbmanufactItemsPeriods:

					if self.manufactItemsPeriods[i][j] == -1:
						if j == 0:
							lbound = 0
						else:
							lbound = self.manufactItemsPeriods[i][j-1]
						
						zeroperiods = []
						k = lbound+1
						while k <= Chromosome.problem.deadlineDemandPeriods[i][j]:
							if self._solution[k] == 0:
								zeroperiods.append(k)
							k+=1

						#print("zeroperiods : ", zeroperiods)
						nbZeroPeriods = len(zeroperiods)

						if nbZeroPeriods > 0:

							cost1 = Chromosome.getCostof(zeroperiods[0], i+1, j+1, copy_solution)
							#print(" cost1 : ", cost1 )

							k = 1 
							indice = zeroperiods[0]
							while k < nbZeroPeriods:
								cost2 = Chromosome.getCostof(zeroperiods[k], i+1, j+1, copy_solution)
								#print(" cost2 : ", cost2 , zeroperiods[k])
								if cost2 < cost1:
									#print(" cost2 < cost1 : ", cost1 , cost2 )
									indice = zeroperiods[k]
								k+=1

							self._solution[indice] = i+1

							itemDemandPeriods = self.manufactItemsPeriods[i]
							itemDemandPeriods[j] = indice

						else:
							
							# experimental code 

							# if there's no place to put this item, then i check all the other times in order to put this item there
							
							lbound = 0
							p = 1
							for deadline in Chromosome.problem.deadlineDemandPeriods[i]:

								zeroperiods = []
								k = lbound
								while k <= deadline:
									if self._solution[k] == 0:
										zeroperiods.append(k)
									k += 1
								lbound = deadline + 1

								nbZeroPeriods = len(zeroperiods)
								if nbZeroPeriods > 0:

									cost1 = Chromosome.getCostof(zeroperiods[0], i+1, p, copy_solution)
									#print(" cost1 : ", cost1 )

									k = 1 
									indice = zeroperiods[0]
									while k < nbZeroPeriods:
										cost2 = Chromosome.getCostof(zeroperiods[k], i+1, p, copy_solution)
										#print(" cost2 : ", cost2 , zeroperiods[k])
										if cost2 < cost1:
											#print(" cost2 < cost1 : ", cost1 , cost2 )
											indice = zeroperiods[k]
										k+=1

									self._solution[indice] = i+1

									itemDemandPeriods = self.manufactItemsPeriods[i]
									itemDemandPeriods[j] = indice

									break
								p += 1


					j+=1
				i+=1

		#print("at the end of getFeasible : ", self._solution)
		self._get_fitnessValue()
		self._get_itemsRanks()

	def getCostof(cls, indice, item, rank,solution, secondIndice = -1):

		solution = list(solution)

		if secondIndice != -1:
			solution[secondIndice] = 0

		cost = 0
		# stocking cost 
		deadline = cls.problem.deadlineDemandPeriods[item-1][rank-1]
		cost += (deadline - indice)* int(Chromosome.problem.holdingGrid[item-1])

		#print(" cost 1 : ", cost)

		# change-over cost 
		nItem, nIndice = nextPeriodItemOf(indice, solution)
		#print(" nItem : ", nItem, " nIndice : ", nIndice)
		if nItem != 0:
			cost += int(Chromosome.problem.chanOverGrid[item-1][nItem-1])

		#print(" cost 2 : ", cost)

		pItem, pIndice = previousPeriodItemOf(indice, solution)
		#print(" pItem : ", pItem, " pIndice : ", pIndice, " solution : ", solution)
		if pItem != 0:
			cost += int(Chromosome.problem.chanOverGrid[pItem-1][item-1])

		#print(" cost 3 : ", cost)

		return cost

	# Class' methods
	getCostof = classmethod(getCostof)

	# Properties
	solution = property(_get_solution,_set_solution)
	fitnessValue = property(_get_fitnessValue,_set_fitnessValue)
	itemsRank = property(_get_itemsRanks, _set_itemsRanks)
	hashSolution = property(_get_hashSolution, _set_hashSolution) 

class Node(object):

	def __init__(self):

		self._solution = [] #[0] * (Chromosome.problem.nbTimes * 2)
		self._currentItem = 0
		self._currentPeriod = 0
		self._currentQuantity = 0
		self.fitnessValue = 0
		self.remPeriods = []

		self.remItems = []
		for i in range(Chromosome.problem.nbItems):
			self.remItems.append(i + 1)

	def __repr__(self):
		return "Chromosome : " + str(self._solution) + ", " + str(self.fitnessValue)
		#" Current Item : " + str(self.currentItem) + " Current Period : " + str(self.currentPeriod) + " Item Counter : " + str(self.itemCounter) + " Fitness value : " + 

	def isLeaf(self):
		
		#if self.itemCounter == Chromosome.problem.nbItems and self.currentPeriod == len(Chromosome.problem.deadlineDemandPeriods[self.currentItem-1]):
		#	return True
		#return False

		if self.remItems == [] and self.remPeriods == [] and self._currentQuantity == Chromosome.problem.deadlineDemandPeriods[self._currentItem-1][self._currentPeriod-1][1]:
			return True
		return False

	def getChildren(self):
		
		children = []
		#print(self.remItems, " : ",self.remPeriods)

		nextItem = 0
		nextPeriod = 0
		nextQuantity = 0

		#print("log getChildren 0 : ", self._currentItem, " : ", self._currentPeriod, " : ", self._currentQuantity, " : ", self.remPeriods)#Chromosome.problem.deadlineDemandPeriods, " : ", Chromosome.problem.deadlineDemandPeriods[self._currentItem-1][self._currentPeriod-1][1])
		# i produce the successors of this current node
		if self._currentQuantity != Chromosome.problem.deadlineDemandPeriods[self._currentItem-1][self._currentPeriod-1][1]:

			#print("yes")
			nextItem = self._currentItem
			nextPeriod = self._currentPeriod
			nextQuantity = self._currentQuantity + 1

		elif self.remPeriods != []:

			nextItem = self.currentItem
			nextPeriod = self.remPeriods[randint(0, len(self.remPeriods)-1)]
			nextQuantity = 1 #Chromosome.problem.deadlineDemandPeriods[nextItem-1][nextPeriod-1][1]

		elif self.remItems != []:

			#self.remPeriods.remove(self._currentPeriod)

			nextItem = self.remItems[randint(0, len(self.remItems)-1)]
			nextPeriod = randint(1, len(Chromosome.problem.deadlineDemandPeriods[nextItem-1]))
			nextQuantity = 1 #Chromosome.problem.deadlineDemandPeriods[nextItem-1][nextPeriod-1][1]

		if nextItem != 0:

			#print("log getChildren : ", nextItem, " : ", nextPeriod, " : ", nextQuantity)
			children = list(self.putNextItem(nextItem, nextPeriod, nextQuantity))
			#print(self.queue)

		return children

	def buildRemPeriod(self):

		self.remPeriods = []
		for i in range(len(Chromosome.problem.deadlineDemandPeriods[self._currentItem-1])):
			self.remPeriods.append(i + 1)

	def putNextItem(self, nextItem, nextPeriod, nextQuantity):

		period = Chromosome.problem.deadlineDemandPeriods[nextItem-1][nextPeriod-1][0]
		#print("period : ", period)
		childrenQueue = []

		for it in range(Chromosome.problem.nbCapacities):

			i = period + it * Chromosome.problem.nbTimes
			while i >= 0 + it * Chromosome.problem.nbTimes:

				if self.solution[i][0] == 0 or self.solution[i][0] == nextItem: 

					solution = copy.deepcopy(self.solution)
					if solution[i][0] == 0 :
						solution[i][0] = nextItem
						solution[i][1] = 1
					elif solution[i][0] == nextItem:
						solution[i][1] += 1

					nextNode = copy.deepcopy(self)
					if nextItem != self.currentItem:
						nextNode.currentItem = nextItem
						#nextNode.buildRemPeriod()
						nextNode.currentPeriod = nextPeriod
						#print("ok 1 ")
					elif nextPeriod != self.currentPeriod:
						nextNode.currentPeriod = nextPeriod
						#print("ok 2 ")
					else:
						nextNode.currentQuantity = nextQuantity
						#print("ok 3 ")#: ", nextNode.currentQuantity, nextQuantity)
						#print(self._currentItem, self._currentPeriod)
						if nextNode.currentQuantity == Chromosome.problem.deadlineDemandPeriods[self._currentItem-1][self._currentPeriod-1][1]:
							nextNode.remPeriods.remove(nextPeriod)

					nextNode.solution = solution
					#print(" i : ", i, " node : ", nextNode, " yep : ", nextNode.currentQuantity, nextQuantity)

					nbChildren = len(childrenQueue)

					#print("Child Node : ", nextNode)

					if (childrenQueue == []):

						childrenQueue.append(copy.deepcopy(nextNode))
						#print(threadQueue)
				
					elif nbChildren == 1 and (childrenQueue[0]).fitnessValue == 0:

						childrenQueue.append(copy.deepcopy(nextNode))
						#print(threadQueue)

					else:
						# i sort the list of zeroperiods from the most convenient place to the least convenient one
						prevValue = 0
						j = 0
						found = False
						while j < nbChildren:

							if nextNode.fitnessValue >= prevValue and nextNode.fitnessValue <= (childrenQueue[j]).fitnessValue:
								found = True
								childrenQueue = childrenQueue[:j] + [copy.deepcopy(nextNode)] + childrenQueue[j:]
								break

							prevValue = (childrenQueue[j]).fitnessValue

							j += 1

						if found is False:
							childrenQueue.append(copy.deepcopy(nextNode))

				i -= 1

		#print("childrenQueue : ", list(reversed(childrenQueue)), "---")
		return reversed(childrenQueue)
		#print(self.queue, "---")
	

	def _get_currentItem(self):
		return self._currentItem

	def _set_currentItem(self, new_value):
		self._currentItem = new_value
		self.remItems.remove(new_value)
		self.buildRemPeriod()
		#self._currentQuantity = 0 

	def _get_currentPeriod(self):
		return self._currentPeriod

	def _set_currentPeriod(self, new_value):
		self._currentPeriod = new_value
		self._currentQuantity = 1
		if self.currentQuantity == Chromosome.problem.deadlineDemandPeriods[self._currentItem-1][self._currentPeriod-1][1]:
			self.remPeriods.remove(new_value)
		#else:
		#	self._currentQuantity += 1

	def _get_solution(self):
		return self._solution

	def _set_solution(self, new_value):
		self._solution = list(new_value)
		self.fitnessValue = Node.evaluate(self._solution) 

	def _set_currentQuantity(self, new_value):
		self._currentQuantity = new_value

	def _get_currentQuantity(self):
		return self._currentQuantity

	def evaluate(cls, sol):
			
		solution = list(sol)
		#print(solution)

		fitnessValue = 0
		grid = Chromosome.problem.chanOverGrid

		# Calculation of all the change-over costs

		for i in range(Chromosome.problem.nbCapacities):

			prevItem = 0
			for j in range(Chromosome.problem.nbTimes):

				#print(" intermediary cost : ", fitnessValue)
				item = solution[j + i * Chromosome.problem.nbTimes][0]
				if j == 0:
					if item != 0:
						fitnessValue += int(grid[i][item-1])
						#print(1)
				else:
					if item != 0 and prevItem != item:
						fitnessValue += int(grid[i][item-1])
						#print(2, int(grid[i][item-1]))

				#print(" j : ", j, " item : ", item)
				prevItem = item

		#print("log evaluate : ", fitnessValue)
				
		#print(" intermediary cost : ", fitnessValue)
		# Calculation of the sum of holding costs

		for i in range(Chromosome.problem.nbCapacities):

			itemCounter = [1] * Chromosome.problem.nbItems

			for j in range(Chromosome.problem.nbTimes):

				item = solution[j + i * Chromosome.problem.nbTimes][0]
				quantity = solution[j + i * Chromosome.problem.nbTimes][0]
				#print ("log evaluate 1 : ", item, " : ", Chromosome.problem.deadlineDemandPeriods)

				if item != 0:

					#print("log evaluate fitness : ", fitnessValue, " : ", int(Chromosome.problem.holdingGrid[item-1]), " : ", (Chromosome.problem.deadlineDemandPeriods[item-1][(itemCounter[item-1])-1][0] - j), " : ", quantity, " : ")
					fitnessValue += int(Chromosome.problem.holdingGrid[item-1])	\
					* (Chromosome.problem.deadlineDemandPeriods[item-1][(itemCounter[item-1])-1][0] - j) \
					* quantity

				it = 0
				for itemDeadlines in Chromosome.problem.deadlineDemandPeriods:
					for deadline in itemDeadlines:
						if deadline[0] == j:
							#print ("log evaluate 2 : ", item, " : ", j)
							itemCounter[it] += 1 
					it += 1

				#print ("log evaluate 2 : ", itemCounter)

		#print(fitnessValue)
		
		return fitnessValue

	evaluate = classmethod(evaluate)

	# Properties
	currentItem = property(_get_currentItem, _set_currentItem)
	currentPeriod = property(_get_currentPeriod, _set_currentPeriod)
	solution = property(_get_solution, _set_solution)
	currentQuantity = property(_get_currentQuantity, _set_currentQuantity)
