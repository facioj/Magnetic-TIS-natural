#!/usr/bin/python

import sys,os


def list_files(directory = "./", site = 0):
   all_files_in_diretory = os.listdir(directory)

   if site < 10:
       stringsite = """site00%(site)s"""%locals()
   else:
       stringsite = """site0%(site)s"""%locals()

   files_of_site = ["""%(directory)s/%(a)s"""%locals() for a in all_files_in_diretory if (stringsite in a and "ldos" and "lj" in a)]

   return files_of_site 

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
   file = open(output,'w')
   for key  in dictionary:
       print >> file, eval(key),dictionary[key]
   file.close()

class dos_sum:

    def sum_site(self,site):
        files_site = list_files(directory = self.directory,site = site)
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


    def __init__(self,directory,sites):
        self.dicts = {}
        self.directory = directory
        self.sites = sites

        for site in self.sites:
            self.dicts["""%(site)s"""%locals()] = self.sum_site(site)


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


