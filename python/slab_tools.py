#! /usr/bin/env python


to_angs = 0.529177

class my_slab():

    def __init__(self,file_name,rel=False):

        self.rel = rel
        self.N_atoms = -1

	#orbital names and labels
	self.orbs = [['p+1up','p+1dn'], ['p-1up','p-1dn'], ['p+0up','p+0dn'] ]
        self.orbs_names = ['px','py','pz']
	if(self.rel):
		self.orbs = [['p3/2+3/2'],['p3/2-3/2'],['p3/2+1/2'],['p3/2-1/2'],['p1/2+1/2'],['p1/2-1/2']]
		self.orbs_names = ['32_32','32_-32','32_12','32_-12','12_12','12_-12']

	self.labd = ['d+2up','d+2dn','d+1up','d+1dn','d+0up','d+0dn','d-1up','d-1dn','d-2up','d-2dn']

	file = open(file_name,'r')
	data = file.readlines()
	file.close()

	i = 0
	block = []
	for fila in data:
	    if "lattice constants [aB]" in fila:
		vals = fila.split()
		self.c = eval(vals[-1])
	    if "Atom sites" in fila:
                sp = data[i+2]
                vals = sp.split()
                self.N_atoms = eval(vals[-1])
		block = data[i+6:i+6+self.N_atoms]
		break
	    i+=1

        print "c: ", self.c
        print "Number of atoms in the slab: ", self.N_atoms

	self.my_list = []
	for b in block:
		sp = b.split()
		z=eval(sp[6])
	        if(z<=0):
		  	z += self.c
		self.my_list.append([eval(sp[0]),sp[1],[eval(sp[4]),eval(sp[5]),z]])

	self.my_list.sort(key = lambda t: t[2][2])


    def identify_blocks(self):
        """
        This method for the moment is only for eight-layer MnBi4Te7.
        """
	Q1 = self.my_list[0:5] 
	S1 = self.my_list[5:12] 
	Q2 = self.my_list[12:17] 
	S2 = self.my_list[17:24] 
	Q3 = self.my_list[24:29] 
	S3 = self.my_list[29:36] 
	Q4 = self.my_list[36:41] 
	S4 = self.my_list[41:48]

	self.my_blocks = [Q1,S1,Q2,S2,Q3,S3,Q4,S4]

	for block in self.my_blocks:
		if(len(block)==5):
			print "Q-block"
		else:
			print "S-block"
		for i in range(len(block)):
			print "Atom: ",block[i][0],block[i][2][2]

	return self.my_blocks

    def identify_blocks_124(self):
        """
        This method for the moment is only for six-layer MnBi2Te4.
        """
	S1 = self.my_list[0:7] 
	S2 = self.my_list[7:14] 
	S3 = self.my_list[14:21] 
	S4 = self.my_list[21:28] 
	S5 = self.my_list[28:35] 
	S6 = self.my_list[35:42] 

	self.my_blocks = [S1,S2,S3,S4,S5,S6]

	for block in self.my_blocks:
		print "S-block"
		for i in range(len(block)):
			print "Atom: ",block[i][0],block[i][2][2]

	return self.my_blocks


    def make_string_def(self,site,element):

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

    def make_string_def_for_dos(self,site,element):

	if(element == "Te" or element=="Bi"):

		if(element == "Te"):
			n=5
		else:
			n=6

		pstring="""
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)ss
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)sp
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
bwdef simple on
    contrib
        site %(site)s
        orbital 4s
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.

		"""%locals()
		return dstring

    def make_string_def_relativistic(self,site,element):

	if(element == "Te" or element=="Bi"):

		if(element == "Te"):
			n=5
		else:
			n=6

		pstring="""
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)ss
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)sp3/2+3/2
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)sp3/2-3/2
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)sp3/2+1/2
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)sp3/2-1/2
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)sp1/2+1/2
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.
bwdef simple on
    contrib
        site %(site)s
        orbital %(n)sp1/2-1/2
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
bwdef simple on
    contrib
        site %(site)s
        orbital 4s
        xaxis 1. 0. 0.
        zaxis 0. 0. 1.

		"""%locals()
		return dstring

    def print_bwdef(self,for_dos=False):

	name_file = "=.mydef"
	if(self.rel):
		name_file = "=.mydef_rel"
	if(for_dos):
		name_file = "=.mydef_dos"
	file = open(name_file,'w')
	print >> file, "coefficients +coeff\n"
	print >> file, "bweights +bweights_mydef\n"

	for b in self.my_blocks:
		for atom in b:
                    if(not for_dos):
			if(self.rel):
				print>> file,self.make_string_def_relativistic(atom[0],atom[1])
			else:
				print>> file,self.make_string_def(atom[0],atom[1])
		    else:
			print>> file,self.make_string_def_for_dos(atom[0],atom[1])
	file.close()

    def layer_weights(self,only_Bi=False,only_Te=False):
	import pyfplo.common as com
	wds=com.WeightDefinitions()

	labp = ['p+1up','p+1dn','p+0up','p+0dn','p-1up','p-1dn']
	if(self.rel):
		labp = ['p3/2+3/2','p3/2-3/2','p3/2+1/2','p3/2-1/2','p1/2+1/2','p1/2-1/2']
	b_ind=0
	for block in self.my_blocks:
		length = len(block)
		print length
                if(length==5):
			name = """Q_%(b_ind)s"""%locals()
		else:
			name = """S_%(b_ind)s"""%locals()
		w=wds.add(name=name)
		labels_block = []
		for atom in block:
			site_index = atom[0]
			site = "{:03d}".format(site_index)
			if(atom[1]=='Bi' and not only_Te):
				for suffix in labp:
					label = """Bi(%(site)s)6%(suffix)s"""%locals()
					labels_block.append(label)
			if(atom[1]=='Te' and not only_Bi):
				for suffix in labp:
					label = """Te(%(site)s)5%(suffix)s"""%locals()
					labels_block.append(label)
			if(atom[1]=='Mn' and self.rel==False and not only_Te and not only_Bi):
				for suffix in self.labd:
					label = """Mn(%(site)s)3%(suffix)s"""%locals()
					labels_block.append(label)
        	w.addLabels(labels=labels_block,fac=1)
		print "w",w
		b_ind+=1
        bw=com.BandWeights('+bweights_mydef')
	output = '+bwsum_layer_resolved'
	if(only_Bi):
		output+="_onlyBi"
	if(only_Te):
		output+="_onlyTe"
        bw.addWeights(wds,output)

    def orbital_layer_weights(self,ewindow=None):
	import pyfplo.common as com

	b_ind = 0
	for i in range(len(self.orbs)):
		orb_name = self.orbs_names[i]
		labp = self.orbs[i]
		wds=com.WeightDefinitions()
		for block in self.my_blocks:
			length = len(block)
                	if(length==5):
				name = """Q_%(b_ind)s"""%locals()
			else:
				name = """S_%(b_ind)s"""%locals()
			w=wds.add(name=name)
			labels_block = []
			for atom in block:
				site_index = atom[0]
				site = "{:03d}".format(site_index)
				if(atom[1]=='Bi'):
					for suffix in labp:
						label = """Bi(%(site)s)6%(suffix)s"""%locals()
						labels_block.append(label)
				if(atom[1]=='Te'):
					for suffix in labp:
						label = """Te(%(site)s)5%(suffix)s"""%locals()
						labels_block.append(label)
        		w.addLabels(labels=labels_block,fac=1)
			b_ind += 1
        	bw=com.BandWeights('+bweights_mydef')
		if(ewindow==None):
        		bw.addWeights(wds,"""+bwsum_orb_%(orb_name)s_layer_resolved"""%locals())
		else:
			e_min = ewindow[0]
			e_max = ewindow[1]
        		bw.addWeights(wds,"""+bwsum_orb_%(orb_name)s_emin_%(e_min)s_emax_%(e_max)s_layer_resolved"""%locals(),ewindow=ewindow)

    def build_bulk_weight(self):
        self.weight = {}
        file_name = """weight_bulk.dat"""%locals()
        file = open(file_name,'w')

        for i in range(49):
            self.weight[i] = 0

        ind_atom = 0
        blocks_bulk = [3,4]
        for i_block in blocks_bulk:
                block = self.my_blocks[i_block]
		for atom in block:
			site_index = atom[0]
			z_coord = atom[2][2]
			self.weight[site_index] = 1
			print >> file, site_index,z_coord,1.0
			ind_atom +=1
	file.close()
	return self.weight


    def build_exp_decay_weight(self,lambda_0,from_septuple=False):
	"""
	Function that builds the exponential decay weight. For the moment, for the dafult is from the quintuple layer.

	Returns: a dict from atom index to weight
	"""
	import numpy as np
	z_coord_0 = -1000
	self.weight = {}
        file_name = """weight_lambda_%(lambda_0)s_from_QL.dat"""%locals()
        if(from_septuple):
            file_name = """weight_lambda_%(lambda_0)s_from_SL.dat"""%locals()

        file = open(file_name,'w')

        ind_atom = 0
        if(from_septuple==False):        
            for block in self.my_blocks:
		for atom in block:
			site_index = atom[0]
			z_coord = atom[2][2]
			print "WARNING: this only works for quintuple layer termination"
			if(ind_atom == 0):
				z_coord_0 = z_coord
			assert z_coord_0 != -1000

			rel_z = (z_coord - z_coord_0)/lambda_0
			weight = np.exp(-rel_z)
			self.weight[site_index] = weight
			print >> file, site_index,z_coord,rel_z,np.exp(-rel_z)
			ind_atom +=1
        else:
            for block in reversed(self.my_blocks):
		for atom in reversed(block):
			site_index = atom[0]
			z_coord = atom[2][2]
			print "WARNING: this only works for quintuple layer termination"
			if(ind_atom == 0):
				z_coord_0 = z_coord
			assert z_coord_0 != -1000

			rel_z = abs(z_coord - z_coord_0)/lambda_0
			weight = np.exp(-rel_z)
			self.weight[site_index] = weight
			print >> file, site_index,z_coord,rel_z,np.exp(-rel_z)
			ind_atom +=1
	file.close()
	return self.weight

    def orbital_exp_decay(self,lambda_0=20,only_Bi=False,only_Te=False,from_septuple=False):
	"""
	lamba_0: characteristic lenght for the exponential decay in Bohr
	"""
        import pyfplo.common as com

	self.build_exp_decay_weight(lambda_0,from_septuple)

        for i in range(len(self.orbs)):
                orb_name = self.orbs_names[i]
                labp = self.orbs[i]
                wds=com.WeightDefinitions()
		ind_atom = 0
                for block in self.my_blocks:
                        for atom in block:
			    if(atom[1]!='Mn'):
                        	labels_atom = []
                                site_index = atom[0]
                		w=wds.add(name="""exp_%(orb_name)s_%(site_index)s"""%locals())
				
				weight = self.weight[site_index]

				#determine label
                                site = "{:03d}".format(site_index)
                                if(atom[1]=='Bi' and not only_Te):
                                        for suffix in labp:
                                                label = """Bi(%(site)s)6%(suffix)s"""%locals()
                                                labels_atom.append(label)
                                if(atom[1]=='Te' and not only_Bi):
                                        for suffix in labp:
                                                label = """Te(%(site)s)5%(suffix)s"""%locals()
                                                labels_atom.append(label)
				ind_atom += 1
				if(len(labels_atom)>0):
	                        	w.addLabels(labels=labels_atom,fac=weight)
				print "w",w
                bw=com.BandWeights('+bweights_mydef')
		output = """+bwsum_orb_%(orb_name)s_lambda_%(lambda_0)s"""%locals()
		if(only_Bi):
			output += "_onlyBi"
		if(only_Te):
			output += "_onlyTe"
                bw.addWeights(wds,output)

	
    def dos_exp_decay(self,atoms,dos_folder=".",lambda_0=20):
	
	weights = self.build_exp_decay_weight(lambda_0)

	import dos_sum as DS
	my_sum = DS.dos_sum(dos_folder,atoms)
	my_sum.weighted_sum(atoms,weights)

    def orbital_dos_exp_decay(self,atoms,orbitals,dos_folder=".",lambda_0=20,type_file="lj",from_septuple=False):
	
	weights = self.build_exp_decay_weight(lambda_0,from_septuple)

	import dos_sum as DS
	my_sum = DS.dos_sum(dos_folder,atoms,orbitals,type_file=type_file)
	my_sum.weighted_orb_sum(atoms,orbitals,weights)

    def orbital_dos_bulk(self,atoms,orbitals,dos_folder=".",type_file="lj"):
	
	weights = self.build_bulk_weight()

	import dos_sum as DS
	my_sum = DS.dos_sum(dos_folder,atoms,orbitals,type_file=type_file)
	my_sum.weighted_orb_sum(atoms,orbitals,weights)
