#!/usr/bin/python

import sys,os


def list_files(directory = "./", site = 0,type_file="lj"):
   all_files_in_diretory = os.listdir(directory)

   if site < 10:
       stringsite = """site00%(site)s"""%locals()
   else:
       stringsite = """site0%(site)s"""%locals()

   files_of_site = ["""%(directory)s/%(a)s"""%locals() for a in all_files_in_diretory if (stringsite in a and "dos" in a and type_file in a)]
   if(type_file == "lj"):
	return files_of_site
   else:
	new_file_set = []
	for e in files_of_site:
		if("lj" not in e):
			new_file_set.append(e)
	print new_file_set
	return new_file_set


def filter_files_of_site(files_of_site,orbital_name,bdos=False):
   filtered_set = []
   print "orbital name", orbital_name
   for file_ in files_of_site:
	file = open(file_,'r')
        data = file.readlines()
   	file.close()
	if orbital_name in data[0]:
		if(bdos==True):
			if "bdos" in file_:
				filtered_set.append(file_)
		else:
			filtered_set.append(file_)
   print filtered_set	
   return filtered_set

def sum_of_file(file_name,previous_dict=None):
   file = open(file_name,'r')
   data = file.readlines()
   file.close()

   new_block = False
   dict_en_dos = {}
   if(previous_dict != None):
       dict_en_dos = previous_dict

   for fila in data[1:]:
       if("#" in fila):
           new_block = True
       vals = fila.split()
       if(len(vals) == 2):
           if(new_block == False and previous_dict==None):
              dict_en_dos[vals[0]] = eval(vals[1])
           else:
              dict_en_dos[vals[0]] += eval(vals[1])

   return dict_en_dos

def print_dict(dictionary,output):
    if('/' in output):
	new = list(output)
	i=0
	for p in new:
		if(p=='/'):
			new[i]=''
		i+=1
		 
	output=''.join(new)
    file = open(output,'w')
    keys = dictionary.keys()
    sorted_keys = sorted(keys,key=lambda x: eval(x))
    for key  in sorted_keys:
       print >> file, eval(key),dictionary[key]
    file.close()

class dos_sum:

    def sum_site_orbital_bdos(self,site,orbital_name,type_="lj"):
	
        files_site = list_files(directory = self.directory,site = site,type_file=self.type_file)
	files_site_orbital = filter_files_of_site(files_site,orbital_name,bdos=True)
	if(len(files_site_orbital) > 0):
	        print "\nFiles to sum: ", files_site_orbital,"\n"

        	dic_site_orbital = {}
	        i = 0
        	for one_file in files_site_orbital:
	            print " Adding --> ", one_file,"\n"

        	    if(i==0):
                	dic_site_orbital = sum_of_file(one_file)
	            else:
        	        dic_site_orbital = sum_of_file(one_file,previous_dict=dic_site)

	            i+=1

		return dic_site_orbital
	else:
		return None

    def sum_site(self,site):
        files_site = list_files(directory = self.directory,site = site,type_file=self.type_file)
        print "\nFiles to sum: ", files_site,"\n"

        dic_site = {}
        i = 0
        for one_file in files_site:
            print " Adding --> ", one_file,"\n"

            if(i==0):
                dic_site = sum_of_file(one_file)
            else:
                dic_site = sum_of_file(one_file,previous_dict=dic_site)

            i+=1

        print_dict(dic_site,"""total_site_%(site)s.dat"""%locals())

        return dic_site


    def __init__(self,directory,sites,orbitals=None,type_file="lj"):
        self.dicts = {}
        self.directory = directory
        self.sites = sites
	self.type_file = type_file
	print "Orbitals: ",orbitals
        for site in self.sites:
		if(orbitals==None):
	            self.dicts["""%(site)s"""%locals()] = self.sum_site(site)
		else:
		    for orbital in orbitals:
			ad = self.sum_site_orbital_bdos(site,orbital)
			if(ad != None):
		            self.dicts["""%(site)s_%(orbital)s"""%locals()] = ad


    def sum(self,list_of_sites):

        file_name = "sum"
        
        i = 0
        dict_sum = {}
        for site in list_of_sites:
            label = """%(site)s"""%locals()
            if(i==0):
                dict_sum = self.dicts[label]
            else:
                for key  in dict_sum:
                    dict_sum[key] += self.dicts[label][key]
            i+=1
            file_name += "_" + label
        file_name += ".dat"
        print_dict(dict_sum,file_name)

        return dict_sum

    def weighted_sum(self,list_of_sites,weight):

        i = 0
        dict_sum = {}
        for site in list_of_sites:
            label = """%(site)s"""%locals()
	    w = weight[site]
            if(i==0):
                dict_sum = self.dicts[label]
		for key  in dict_sum:
			dict_sum[key] *= w
            else:
                for key  in dict_sum:
                    dict_sum[key] += self.dicts[label][key] * w
            i+=1

        file_name = "weighted_dos.dat"

        print_dict(dict_sum,file_name)

        return dict_sum

    def weighted_orb_sum(self,list_of_sites,orbitals,weight):

	my_keys = self.dicts.keys()
	print "my_keys: ", my_keys
        for orbital in orbitals:
            i = 0
            dict_sum = {}
            for site in list_of_sites:
	            label = """%(site)s_%(orbital)s"""%locals()
		    if label in my_keys:
			    w = weight[site]
		            if(i==0):
        		        dict_sum = self.dicts[label]
				for key  in dict_sum:
					dict_sum[key] *= w
		            else:
        		        for key  in dict_sum:
                		    dict_sum[key] += self.dicts[label][key] * w
		            i+=1

            file_name = """weighted_%(orbital)s_dos.dat"""%locals()
            print_dict(dict_sum,file_name)

