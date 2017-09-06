

import logging
import numpy as np
import re
import itertools as it

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s')

class Stoichiometry(dict):

	def __init__( self, count=None):
		self.count = count

		

class Stoichiometries(list):

	def __init__(self, l=[]):
		for i in l:
			self.append(i)

	def slim(self, number):
		l=[ x for x in self if len(x) == number]
		stoichiometries_slim = Stoichiometries(l=l)
		return stoichiometries_slim

	def read_stoichiometries(self, file_handle):

		read_list = []
		for line in file_handle.readlines():
			print line
			s, cnt = line.split()
			s_list = re.split('([A-Z])',s)
			stoich = Stoichiometry(count=int(cnt))
			for i in range(1, len(s_list)):
				if i%2 == 0:
					continue
				#kdrew: TODO deal with weird pdb format that has '?' for chain numbers
				print i, s_list[i], s_list[i+1]
				if s_list[i+1] == '':
					s_list[i+1] = 1

				stoich[s_list[i]] = int(s_list[i+1])
				
			read_list.append(stoich)

		for stoich in read_list:
			print stoich
			permut_set = set()
			for permut in it.permutations(stoich.values(), len(stoich)):
				permut_set.add(permut)
			for s in permut_set:
				#kdrew: split counts by how many permutations there are
				#kdrew: this is trickier than expected, splitting counts might not be the best idea
				new_stoich = Stoichiometry(count=(1.0*stoich.count/len(permut_set)))
				for i, key in enumerate(stoich):
					new_stoich[key] = s[i]
				self.append(new_stoich)


	def total_count(self, numOfProteins=None):
		
		total = 0.0
		for stoich in self:
			if numOfProteins == None or len(stoich) == numOfProteins:
				total += stoich.count

		return total




