import __main__
__main__.pymol_argv = [ 'pymol', '-qei' ]

import pymol
pymol.finish_launching()
cmd = pymol.cmd

import symmetry_axis as sa
import center_of_mass as cm
from pymol.cgo import *

class SubComplex(object):
	def __init__(self, prot_ids, translate_vector, color, connection_length=20.0):
		self.prot_ids = prot_ids
		self.translate_vector = translate_vector
		self.color = color
		self.init_com = self.get_com()
		self.conn_length = connection_length

	def get_com(self, ):
		sele = " + ".join(self.prot_ids) 
		cmd.select("sele", "%s" % sele)
		center1 = cm.get_com("sele")
		return center1

class Counter:
	def __init__(self):
		self.state = 1

counter = Counter()

def main():
	#cmd.fetch("4cr2", "4cr2")

	#cmd.png("~/data/proteasome/pymol_scripts/4cr2.png")
	#cmd.rotate('x', angle=100, selection='all')
	#cmd.png("~/data/proteasome/pymol_scripts/4cr2_rotate.png")


	cmd.load("../4cr2_objects.pse")
	cmd.disable()
	cmd.set_view (" -0.183798626,   -0.982478261,    0.030202165,\
		-0.941125035,    0.184772134,    0.283043444,\
		-0.283658504,    0.023604879,   -0.958611667,\
		0.000000000,    0.000000000, -949.717102051,\
		174.657592773,  310.533782959,  285.625671387,\
		774.580139160, 1124.854980469,  -20.000000000" )


	PNAS=0
	LINKAGE=1

	sc_list = []
	c_list = []

	if LINKAGE:
		sc1 = SubComplex(['a1','a2','a3','a4','a5','a6','a7','b2','b3','b5','b6'], [0,-75,0], "blue", 75.0)
		sc2 = SubComplex(['b1','b4','b7'], [0,-150,0], "red", 150.0)
		sc3 = SubComplex(['rpn9', 'rpn6', 'rpn7', 'rpn8'], [-200,0,0], "green", 200.0)
		sc4 = SubComplex(['rpn3', 'rpn5', 'rpn11'], [-100,0,0], "white", 100.0)
		sc5 = SubComplex(['rpn12'], [-300,0,0], "brown", 300.0)
		sc6 = SubComplex(['rpt4', 'rpt5'], [0,0,-150], "cyan", 150.0)
		sc7 = SubComplex(['rpt6', 'rpt3'], [0,0,-50], "magenta", 50.0)
		sc8 = SubComplex(['rpt1', 'rpt2', 'rpn1'], [50,0,0], "orange", 50.0)
		sc9 = SubComplex(['rpn2'], [0,50,0], "yellow", 50.0)

		sc_list.append(sc1)
		sc_list.append(sc2)
		sc_list.append(sc3)
		sc_list.append(sc4)
		sc_list.append(sc5)
		sc_list.append(sc6)
		sc_list.append(sc7)
		sc_list.append(sc8)
		sc_list.append(sc9)

		c_list.append(sc1)
		c_list.append(sc2)
		c_list.append(sc3)
		c_list.append(sc4)
		c_list.append(sc5)
		c_list.append(sc6)
		#c_list.append(sc7)
		c_list.append(sc8)
		c_list.append(sc9)

	if PNAS:
		sc_list.append(SubComplex(['a1','a2','a3','a4','a5','a6','a7'], [0,-50,0], "red"))
		sc_list.append(SubComplex(['b1','b2','b3','b4','b5','b6','b7'], [0,-100,0], "blue"))
		sc_list.append(SubComplex(['rpt1','rpt2','rpt3','rpt4','rpt5','rpt6'], [0,0,0], "cyan"))
		sc_list.append(SubComplex(['rpn9', 'rpn5', 'rpn6', 'rpn7', 'rpn3', 'rpn12'], [-150,0,0], "green"))
		sc_list.append(SubComplex(['rpn1'], [50,0,0], "orange"))
		sc_list.append(SubComplex(['rpn2'], [0,50,0], "yellow"))
		sc_list.append(SubComplex(['rpn8', 'rpn11'], [-75,0,0], "magenta"))

	
	explode(sc_list)
	draw_connections(c_list)
	

def draw_connections(c_list):

	print "drawing connections"

	for conn in c_list:       
		#sele = " + ".join(conn[0].prot_ids) 
		#print sele
		#cmd.select("sele", "%s" % sele)
		#center1 = cm.get_com("sele")
		#print center1

		#center2 = conn[1]

		#sele = " + ".join(conn[1].prot_ids)	
		#print sele
		#cmd.select("sele","%s" % sele)
		#center2 = cm.get_com("sele")
		#print center2

		center1 = conn.get_com()
		center2 = conn.init_com

		#i = center1[0] - center2[0]
		#j = center1[1] - center2[1]
		#k = center1[2] - center2[2]
		#sa.draw_axis(center1[0], center1[1], center1[2], i, j, k, conn.conn_length )

		width = 1.0
		r = 1.0
		g = 1.0
		b = 1.0

		x1 = center1[0]
		y1 = center1[1]
		z1 = center1[2]
		x2 = center2[0]
		y2 = center2[1]
		z2 = center2[2]

		obj = [
			LINEWIDTH, width,
			BEGIN, LINES,

			COLOR,   r,  g,  b,
			VERTEX, x1, y1, z1,
			VERTEX, x2, y2, z2,

			END
		]

		cmd.load_cgo(obj,'axis'+str(counter.state))
		counter.state += 1




def explode(sc_list):

	for name in cmd.get_names("objects"):
		for sc in sc_list:
			if name in sc.prot_ids:
				cmd.color(sc.color, name)
				cmd.enable(name)
				cmd.translate(sc.translate_vector, name)

def explode_PNAS():

	for name in cmd.get_names("objects"):
		print name


		if name.startswith('a'):
			cmd.color("red", name)
			print "red"
			cmd.translate([0,-50,0], name)
			cmd.enable(name)

		elif name.startswith('b'):
			cmd.color("blue", name)
			print "blue"
			cmd.translate([0,-100,0], name)
			cmd.enable(name)

		elif name.startswith('rpt'):
			cmd.color("cyan", name)
			print "cyan"
			cmd.enable(name)

		elif name in ['rpn9', 'rpn5', 'rpn6', 'rpn7', 'rpn3', 'rpn12']:
			cmd.color("green", name)
			print "green"
			cmd.enable(name)
			cmd.translate([-150,0,0], name)

		elif name in ['rpn1']:
			cmd.color("orange", name)
			cmd.enable(name)
			cmd.translate([50,0,0], name)

		elif name in ['rpn2']:
			cmd.color("yellow", name)
			cmd.enable(name)
			cmd.translate([0,50,0], name)

		elif name in ['rpn8', 'rpn11']:
			cmd.color("magenta", name)
			cmd.enable(name)
			cmd.translate([-75,0,0], name)
	




	#pymol.cmd.quit()

if __name__ == "__main__":
	main()

