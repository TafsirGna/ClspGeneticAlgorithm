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