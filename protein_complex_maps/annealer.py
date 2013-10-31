
#kdrew: base class to schedule alterations to montecarlo temperature 
class Annealer(object):
	
	#kdrew: scale value should be < 1 if decreasing temperature
	def __init__( self, montecarlo, scale = 0.5 ):
		self.montecarlo = montecarlo
		self.scale = scale

	def anneal( self ):
		current_temp = self.montecarlo.temp()
		new_temp = current_temp*self.scale
		self.montecarlo.temp(new_temp)
		
		

#kdrew: annealer class to only accept lower scoring moves
class QuenchAnnealer(Annealer):
	def __init__( self, montecarlo, quench_iteration ):
		super(QuenchAnnealer, self).__init__(montecarlo, scale = 0.0)
		self.quench_iteration = quench_iteration

	def anneal( self ):
		if self.montecarlo.iterations() >= self.quench_iteration:
			super(QuenchAnnealer, self).anneal()


#kdrew: annealer class to keep acceptance rate fixed(-ish)
class RateAnnealer(Annealer):
	def __init__( self, montecarlo, scale = 0.5, rate = 0.33  ):
		super(RateAnnealer, self).__init__( montecarlo, scale )
		self.rate = rate

	def anneal( self ):
		curr_rate = self.montecarlo.accept_rate()
		current_temp = self.montecarlo.temp()
		#kdrew: if the current acceptance rate is equal to target rate, keep temp same
		if curr_rate == self.rate:
			return
		#kdrew: if the current acceptance rate is too high, lower the temperature
		elif curr_rate > self.rate:
			new_temp = current_temp*self.scale
		#kdrew: if the current acceptance rate is too low, raise the temperature (divide by scale which should be between 0 and 1)
		else:
			new_temp = (1.0*current_temp)/self.scale

		self.montecarlo.temp(new_temp)
		


#kdrew: annealer class to raise temperature after set number of iterations
class RampAnnealer(Annealer):
	#kdrew: ramp is the number of iterations to ramp
	def __init__( self, montecarlo, scale = 0.5, ramp = 500 ):
		super(RampAnnealer, self).__init__( montecarlo, scale )
		self.ramp = ramp

	def anneal( self ):
		if self.montecarlo.iterations() % self.ramp == 0 :
			self.montecarlo.reset_temp()
		

class RateQuenchAnnealer(Annealer):
	def __init__( self, montecarlo, quench_iteration, scale = 0.5, rate = 0.33 ):
		self.rate_annealer = RateAnnealer(montecarlo, scale, rate)
		self.quench_annealer = QuenchAnnealer(montecarlo, quench_iteration)

	def anneal( self, ):
		#kdrew: perform rate anneal first
		self.rate_annealer.anneal()
		#kdrew: perform quench last because once quench should stay quenched
		self.quench_annealer.anneal()



