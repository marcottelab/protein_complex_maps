import math as m

class MonteCarlo(object):

	def __init__( self, scorefunction, temp, random_module=None ):
		self.__scorefunction = scorefunction
		self.__temp = temp

		self.__low_scoring_bicluster = None
		self.__low_score = None
		self.__current_bicluster = None
		self.__current_score = None

		self.__accepts = 0
		self.__rejects = 0
		self.__thermal_accepts = 0
		self.__score_accepts = 0
		self.__lowscore_accepts = 0

		self.__result_history = []
		self.__score_history = []
		self.__score_diff_history = []

		self.__iterations = 0 
		self.__iterations_since_last_accept = 0


		if random_module == None:
			try:
				import random 
			except ImportError:
				self.__random_module = None
			else:
				self.__random_module = random
		else:
			self.__random_module = random_module


	def reset(self, ):
		self.__low_scoring_bicluster = None
		self.__low_score = None
		self.__current_bicluster = None
		self.__current_score = None

		self.__accepts = 0
		self.__rejects = 0
		self.__thermal_accepts = 0
		self.__score_accepts = 0
		self.__lowscore_accepts = 0

		self.__result_history = []
		self.__score_history = []
		self.__score_diff_history = []

		self.__iterations = 0 
		self.__iterations_since_last_accept = 0

	def lowscore_bicluster(self,):
		return self.__low_scoring_bicluster

	def lowscore(self,):
		return self.__low_scoring
		
	def current_bicluster(self,):
		return self.__current_bicluster

	def current_score(self,):
		return self.__current_score

	def accept_rate(self,):
		return 1.0*self.__accepts/self.__iterations
	
	def iterations_since_last_accept(self,):
		return self.__iterations_since_last_accept
	
	def accepts(self,):
		return self.__accepts

	def rejects(self,):
		return self.__rejects

	def thermal_accepts(self,):
		return self.__thermal_accepts

	def score_accepts(self,):
		return self.__score_accepts

	def lowscore_accepts(self,):
		return self.__lowscore_accepts

	def result_history(self,):
		return self.__result_history

	def score_history(self,):
		return self.__score_history

	def score_diff_history(self,):
		return self.__score_diff_history

	def boltzmann(self, data_matrix, trial_bicluster):
		trial_score = self.__scorefunction(data_matrix, trial_bicluster)

		if self.__current_score == None:
			score_diff = float('inf')
		else:
			score_diff = trial_score - self.__current_score

		self.__score_diff_history.append( score_diff )


		if( trial_score >= self.__current_score and self.__current_score != None ):
			print "bicluster does not decrease total score, test "
			#kdrew: if no, use logistic regression to predict membership value: x
			#kdrew: p(drop|x) = math.e**(-(1-x)/T)

			if self.__temp == 0.0:
				prob = 0.0
			else:
				prob = m.exp( (-1.0*score_diff)/self.__temp )

			random_prob = self.__random_module.random()

			print "score_diff: %s, prob: %s, random_prob: %s" % (score_diff, prob,random_prob,)

			if prob > self.__random_module.random():
				#kdrew: thermal accept
				self.__current_bicluster = trial_bicluster
				self.__current_score = trial_score
				self.__accepts += 1
				self.__thermal_accepts += 1
				self.__iterations_since_last_accept = 0

				self.__result_history.append( "thermal" )
				self.__score_history.append( trial_score ) 

			else:
				#kdrew: reject
				self.__rejects += 1
				self.__iterations_since_last_accept += 1
				self.__result_history.append( "reject" )
				self.__score_history.append( trial_score )

		else:
			print "bicluster decreases total score, automatically accept"

			if self.__low_score > trial_score or self.__low_scoring_bicluster == None:
				self.__low_scoring_bicluster = trial_bicluster
				self.__low_score = trial_score
				self.__lowscore_accepts += 1

				self.__result_history.append( "lowscore" )
				self.__score_history.append( trial_score )

			else:
				self.__result_history.append( "score" )
				self.__score_history.append( trial_score )

			self.__current_bicluster = trial_bicluster
			self.__current_score = trial_score
			self.__accepts += 1
			self.__score_accepts += 1

			self.__iterations_since_last_accept = 0


		self.__iterations += 1

		return self.__current_bicluster

