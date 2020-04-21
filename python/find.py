#! /usr/bin/env python



def from_parser():
	import pyfplo.fploio as fploio 
	p=fploio.INParser()
	p.parseFile("=.in")
	d=p()
	WP = d('wyckoff_positions')

	my_list = []

	for i in range(48):
   		index = i+1
   		my_list.append([index,WP[i]('element').S,WP[i]('tau').listD])

	my_list.sort(key = lambda t: t[2][2])

	return my_list

def from_out(max_coor,vacuum):
	L=max_coor*2+vacuum
	file = open("out",'r')
	data = file.readlines()
	file.close()
	i = 0
	block = []
	for fila in data:
		if "Atom sites" in fila:
			block = data[i+6:i+6+48]
			break
		i+=1
	my_list = []
	for b in block:
		sp = b.split()
		z=eval(sp[6])
	        if(z<0):
		  	z += L	
		my_list.append([eval(sp[0]),sp[1],[eval(sp[4]),eval(sp[5]),z]])

	my_list.sort(key = lambda t: t[2][2])
	return my_list

def identify_blocks():
	max_coor=113.163
	vac=10
	my_list = from_out(max_coor,vac)
	#my_list2 = from_parser() #this one needs to be changed because of c shift

	Q1 = my_list[0:5] 
	S1 = my_list[5:12] 
	Q2 = my_list[12:17] 
	S2 = my_list[17:24] 
	Q3 = my_list[24:29] 
	S3 = my_list[29:36] 
	Q4 = my_list[36:41] 
	S4 = my_list[41:48]

	my_blocks = [Q1,S1,Q2,S2,Q3,S3,Q4,S4]

	for block in my_blocks:
		if(len(block)==5):
			print "Q-block"
		else:
			print "S-block"
		for i in range(len(block)):
			print block[i][0]

	return my_blocks


def process_weights():
	import pyfplo.common as com
	wds=com.WeightDefinitions()
	"example:"
	w=wds.add(name='all')
	w.addAtoms(atom='Bi',sites=[1],orbitals=['6p'],fac=1)
	print "w",w

def make_string_def(site,element):

	if(element == "Te" or element=="Bi"):

		if(element == "Te"):
			n=5
		else:
			n=6

		pstring="""
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)sp-1
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)sp0
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)sp+1 
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
		"""%locals()
		return pstring
	else:
		dstring="""
bwdef simple on
    contrib
        site %(site)s
        orbital 3d
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
		"""%locals()
		return dstring

def print_bwdef(blocks):

	file = open("=.mydef",'w')
	print >> file, "coefficients +coeff\n"
	print >> file, "bweights +bweights_mydef\n"

	for b in blocks:
		for atom in b:
			print>> file,make_string_def(atom[0],atom[1])
	file.close()

blocks = identify_blocks()
print_bwdef(blocks)
