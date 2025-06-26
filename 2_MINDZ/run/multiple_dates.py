#!/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/python/3.10.2/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 11:49:20 2024

@author: andy
"""

# For special dates :  
# Faire un vecteur special_dates avec toutes les dates
# calculer nb_of_steps = le nombre de dates
# dans la boucle, cur_date = special_dates[id_day]
# mettre l'une ou l'autre des options avec des conditions if

import subprocess
import os
import datetime

os.chdir("/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/MINDZ/run/")

#---- TO CHANGE DEPENDING ON THE RUN WE WANT : ----#
#_________________________________________________________________________________________________


#-- Choose density of particles : 
part_fac = 4


#-- Choose taxon :
taxon = "LH" 
if taxon == "LH" : num_tax = "5"
if taxon == "LF" : num_tax = "7"
if taxon == "YF" : num_tax = "8"


#-- Choose diapause : 
diapause = "0" # 0 = diapause, 1 = active

#-- Choose forward or backward :
bwd = "F" #backward run, T = true
# if bwd = T, put the last date at dstart (ie the most recent). 

time_resol = "day" #hour or day if we want hour or daily output resolution

#-- Choose dates :
lgt_output    = 10  # lenght of each output, in days
date_type     = "cont" #cont for continuous, of special for just certain dates
dstart        = datetime.datetime(2018, 1, 1) #7 15
dend          = datetime.datetime(2018, 12, 31) #8 16
special_dates =  [datetime.datetime(2024, 11, 4)]
#special_dates =  [datetime.datetime(2017, 2, 11), datetime.datetime(2017, 2, 14), datetime.datetime(2017, 3, 8), datetime.datetime(2017, 4, 2), datetime.datetime(2017, 4, 4), datetime.datetime(2017, 4, 13), datetime.datetime(2017, 4, 18), datetime.datetime(2017, 4, 20), datetime.datetime(2017, 4, 29), datetime.datetime(2017, 4, 30), datetime.datetime(2017, 5, 3), datetime.datetime(2017, 5, 6), datetime.datetime(2017, 5, 7), datetime.datetime(2017, 6, 12), datetime.datetime(2017, 6, 13), datetime.datetime(2017, 6, 15), datetime.datetime(2017, 6, 16), datetime.datetime(2017, 6, 17), datetime.datetime(2017, 6, 19), datetime.datetime(2017, 6, 23), datetime.datetime(2017, 6, 24), datetime.datetime(2017, 6, 25), datetime.datetime(2017, 6, 28), datetime.datetime(2017, 6, 30), datetime.datetime(2017, 7, 9), datetime.datetime(2017, 7, 10), datetime.datetime(2017, 7, 11), datetime.datetime(2017, 7, 14), datetime.datetime(2017, 7, 15), datetime.datetime(2017, 7, 16), datetime.datetime(2017, 7, 19), datetime.datetime(2017, 7, 21), datetime.datetime(2017, 7, 22), datetime.datetime(2017, 7, 23), datetime.datetime(2017, 7, 31), datetime.datetime(2017, 8, 1), datetime.datetime(2017, 9, 10), datetime.datetime(2017, 9, 11), datetime.datetime(2017, 9, 18), datetime.datetime(2018, 2, 7), datetime.datetime(2018, 2, 17), datetime.datetime(2018, 3, 31), datetime.datetime(2018, 4, 1), datetime.datetime(2018, 4, 2), datetime.datetime(2018, 4, 4), datetime.datetime(2018, 4, 11), datetime.datetime(2018, 4, 12), datetime.datetime(2018, 4, 17), datetime.datetime(2018, 4, 21), datetime.datetime(2018, 4, 23), datetime.datetime(2018, 4, 24), datetime.datetime(2018, 4, 25), datetime.datetime(2018, 4, 26), datetime.datetime(2018, 5, 1), datetime.datetime(2018, 5, 3), datetime.datetime(2018, 5, 11), datetime.datetime(2018, 5, 20), datetime.datetime(2018, 5, 22), datetime.datetime(2018, 5, 25), datetime.datetime(2018, 5, 27), datetime.datetime(2018, 5, 28), datetime.datetime(2018, 5, 30), datetime.datetime(2018, 6, 1), datetime.datetime(2018, 6, 5), datetime.datetime(2018, 6, 6), datetime.datetime(2018, 6, 7), datetime.datetime(2018, 6, 10), datetime.datetime(2018, 6, 12), datetime.datetime(2018, 6, 16), datetime.datetime(2018, 6, 19), datetime.datetime(2018, 6, 20), datetime.datetime(2018, 6, 25), datetime.datetime(2018, 6, 27), datetime.datetime(2018, 7, 1), datetime.datetime(2018, 7, 3), datetime.datetime(2018, 7, 4), datetime.datetime(2018, 7, 6), datetime.datetime(2018, 7, 9), datetime.datetime(2018, 7, 10), datetime.datetime(2018, 7, 11), datetime.datetime(2018, 7, 12), datetime.datetime(2018, 7, 14), datetime.datetime(2018, 7, 17), datetime.datetime(2018, 7, 18), datetime.datetime(2018, 7, 19), datetime.datetime(2018, 7, 20), datetime.datetime(2018, 7, 21), datetime.datetime(2018, 7, 22), datetime.datetime(2018, 7, 24), datetime.datetime(2018, 7, 25), datetime.datetime(2018, 7, 27), datetime.datetime(2018, 8, 1), datetime.datetime(2018, 8, 2), datetime.datetime(2018, 8, 3), datetime.datetime(2018, 8, 4), datetime.datetime(2018, 8, 6), datetime.datetime(2018, 8, 7), datetime.datetime(2018, 8, 8), datetime.datetime(2018, 8, 9), datetime.datetime(2018, 8, 10), datetime.datetime(2018, 8, 11), datetime.datetime(2018, 8, 12), datetime.datetime(2019, 2, 10), datetime.datetime(2019, 2, 21), datetime.datetime(2019, 2, 27), datetime.datetime(2019, 3, 9), datetime.datetime(2019, 3, 10), datetime.datetime(2019, 3, 17), datetime.datetime(2019, 4, 1), datetime.datetime(2019, 4, 14), datetime.datetime(2019, 4, 15), datetime.datetime(2019, 4, 21), datetime.datetime(2019, 5, 3), datetime.datetime(2019, 5, 9), datetime.datetime(2019, 5, 22), datetime.datetime(2019, 5, 25), datetime.datetime(2019, 5, 26), datetime.datetime(2019, 5, 28), datetime.datetime(2019, 5, 30), datetime.datetime(2019, 5, 31), datetime.datetime(2019, 6, 3), datetime.datetime(2019, 6, 5), datetime.datetime(2019, 6, 23), datetime.datetime(2019, 6, 24), datetime.datetime(2019, 6, 25), datetime.datetime(2019, 6, 27), datetime.datetime(2019, 6, 28), datetime.datetime(2019, 6, 29), datetime.datetime(2019, 6, 30), datetime.datetime(2019, 7, 1), datetime.datetime(2019, 7, 4), datetime.datetime(2019, 7, 5), datetime.datetime(2019, 7, 6), datetime.datetime(2019, 7, 7), datetime.datetime(2019, 7, 8), datetime.datetime(2019, 7, 9), datetime.datetime(2019, 7, 10), datetime.datetime(2019, 7, 11), datetime.datetime(2019, 7, 12), datetime.datetime(2019, 7, 23), datetime.datetime(2019, 7, 24), datetime.datetime(2019, 7, 27), datetime.datetime(2019, 7, 28), datetime.datetime(2019, 7, 29), datetime.datetime(2019, 7, 30), datetime.datetime(2019, 7, 31), datetime.datetime(2019, 8, 1), datetime.datetime(2019, 8, 2), datetime.datetime(2019, 8, 3), datetime.datetime(2019, 8, 4), datetime.datetime(2019, 8, 5), datetime.datetime(2019, 8, 6), datetime.datetime(2019, 8, 7), datetime.datetime(2019, 8, 8), datetime.datetime(2019, 8, 9), datetime.datetime(2019, 8, 10), datetime.datetime(2019, 8, 11), datetime.datetime(2019, 8, 12), datetime.datetime(2019, 8, 14), datetime.datetime(2019, 8, 15), datetime.datetime(2019, 8, 16), datetime.datetime(2019, 8, 20), datetime.datetime(2019, 8, 27), datetime.datetime(2019, 9, 3), datetime.datetime(2019, 9, 8), datetime.datetime(2019, 9, 30), datetime.datetime(2019, 10, 1), datetime.datetime(2019, 10, 4), datetime.datetime(2019, 10, 6), datetime.datetime(2019, 10, 11), datetime.datetime(2019, 10, 12), datetime.datetime(2019, 10, 17), datetime.datetime(2019, 10, 18), datetime.datetime(2019, 10, 19), datetime.datetime(2020, 2, 25), datetime.datetime(2020, 4, 23), datetime.datetime(2020, 4, 28), datetime.datetime(2020, 5, 10), datetime.datetime(2020, 5, 14), datetime.datetime(2020, 5, 22), datetime.datetime(2020, 5, 25), datetime.datetime(2020, 5, 27), datetime.datetime(2020, 5, 29), datetime.datetime(2020, 5, 31), datetime.datetime(2020, 6, 3), datetime.datetime(2020, 6, 4), datetime.datetime(2020, 6, 6), datetime.datetime(2020, 6, 10), datetime.datetime(2020, 6, 12), datetime.datetime(2020, 6, 13), datetime.datetime(2020, 7, 1), datetime.datetime(2020, 7, 6), datetime.datetime(2020, 7, 7), datetime.datetime(2020, 7, 8), datetime.datetime(2020, 7, 9), datetime.datetime(2020, 7, 11), datetime.datetime(2020, 7, 12), datetime.datetime(2020, 7, 13), datetime.datetime(2020, 7, 16), datetime.datetime(2020, 7, 22), datetime.datetime(2020, 7, 23), datetime.datetime(2020, 7, 24), datetime.datetime(2020, 7, 26), datetime.datetime(2020, 7, 28), datetime.datetime(2020, 7, 29), datetime.datetime(2020, 7, 30), datetime.datetime(2020, 7, 31), datetime.datetime(2020, 8, 1), datetime.datetime(2020, 8, 2), datetime.datetime(2020, 8, 3), datetime.datetime(2020, 8, 4), datetime.datetime(2020, 8, 6), datetime.datetime(2020, 8, 7), datetime.datetime(2020, 8, 11), datetime.datetime(2020, 8, 12), datetime.datetime(2020, 8, 13), datetime.datetime(2020, 8, 28), datetime.datetime(2020, 9, 1), datetime.datetime(2020, 9, 2), datetime.datetime(2020, 9, 3), datetime.datetime(2020, 9, 5), datetime.datetime(2020, 9, 6), datetime.datetime(2020, 9, 8), datetime.datetime(2020, 9, 9), datetime.datetime(2020, 9, 11), datetime.datetime(2020, 9, 14), datetime.datetime(2020, 9, 18), datetime.datetime(2020, 9, 21), datetime.datetime(2020, 9, 25), datetime.datetime(2020, 10, 2), datetime.datetime(2020, 10, 3), datetime.datetime(2020, 10, 5), datetime.datetime(2020, 10, 16), datetime.datetime(2020, 10, 18), datetime.datetime(2020, 10, 29), datetime.datetime(2020, 11, 5), datetime.datetime(2020, 12, 1), datetime.datetime(2021, 1, 16), datetime.datetime(2021, 1, 21), datetime.datetime(2021, 2, 3), datetime.datetime(2021, 2, 13), datetime.datetime(2021, 2, 18), datetime.datetime(2021, 2, 21), datetime.datetime(2021, 2, 25), datetime.datetime(2021, 3, 2), datetime.datetime(2021, 3, 6), datetime.datetime(2021, 3, 20), datetime.datetime(2021, 3, 29), datetime.datetime(2021, 4, 15), datetime.datetime(2021, 4, 18), datetime.datetime(2021, 4, 24), datetime.datetime(2021, 4, 27), datetime.datetime(2021, 4, 30), datetime.datetime(2021, 5, 1), datetime.datetime(2021, 5, 2), datetime.datetime(2021, 5, 4), datetime.datetime(2021, 5, 5), datetime.datetime(2021, 5, 7), datetime.datetime(2021, 5, 9), datetime.datetime(2021, 5, 14), datetime.datetime(2021, 5, 19), datetime.datetime(2021, 5, 20), datetime.datetime(2021, 5, 22), datetime.datetime(2021, 5, 23), datetime.datetime(2021, 5, 24), datetime.datetime(2021, 5, 29), datetime.datetime(2021, 5, 30), datetime.datetime(2021, 6, 1), datetime.datetime(2021, 6, 2), datetime.datetime(2021, 6, 7), datetime.datetime(2021, 6, 8), datetime.datetime(2021, 6, 11), datetime.datetime(2021, 6, 13), datetime.datetime(2021, 6, 14), datetime.datetime(2021, 6, 25), datetime.datetime(2021, 6, 27), datetime.datetime(2021, 6, 28), datetime.datetime(2021, 6, 30), datetime.datetime(2021, 7, 1), datetime.datetime(2021, 7, 2), datetime.datetime(2021, 7, 3), datetime.datetime(2021, 7, 7), datetime.datetime(2021, 7, 8), datetime.datetime(2021, 7, 9), datetime.datetime(2021, 7, 10), datetime.datetime(2021, 7, 18), datetime.datetime(2021, 7, 23), datetime.datetime(2021, 7, 25), datetime.datetime(2021, 7, 26), datetime.datetime(2021, 7, 29), datetime.datetime(2021, 7, 30), datetime.datetime(2021, 8, 1), datetime.datetime(2021, 8, 2), datetime.datetime(2021, 8, 3), datetime.datetime(2021, 8, 5), datetime.datetime(2021, 8, 6), datetime.datetime(2021, 8, 8), datetime.datetime(2021, 8, 10), datetime.datetime(2021, 8, 11), datetime.datetime(2021, 8, 15), datetime.datetime(2021, 8, 18), datetime.datetime(2021, 8, 19), datetime.datetime(2021, 8, 21), datetime.datetime(2021, 9, 13), datetime.datetime(2021, 9, 23), datetime.datetime(2021, 9, 25), datetime.datetime(2021, 9, 27), datetime.datetime(2021, 9, 30), datetime.datetime(2021, 10, 2), datetime.datetime(2021, 10, 9), datetime.datetime(2021, 10, 10), datetime.datetime(2021, 10, 16), datetime.datetime(2021, 10, 26), datetime.datetime(2021, 10, 30), datetime.datetime(2021, 10, 31), datetime.datetime(2021, 11, 2), datetime.datetime(2021, 11, 5)]



#special_dates = [datetime.datetime(2021, 1, 16), datetime.datetime(2021, 1, 21), datetime.datetime(2021, 2, 3), datetime.datetime(2021, 2, 13), datetime.datetime(2021, 2, 18), datetime.datetime(2021, 2, 21), datetime.datetime(2021, 2, 25), datetime.datetime(2021, 3, 2), datetime.datetime(2021, 3, 6), datetime.datetime(2021, 3, 20), datetime.datetime(2021, 3, 29), datetime.datetime(2021, 4, 15), datetime.datetime(2021, 4, 18), datetime.datetime(2021, 4, 24), datetime.datetime(2021, 4, 27), datetime.datetime(2021, 4, 30), datetime.datetime(2021, 5, 1), datetime.datetime(2021, 5, 2), datetime.datetime(2021, 5, 4), datetime.datetime(2021, 5, 5), datetime.datetime(2021, 5, 7), datetime.datetime(2021, 5, 9), datetime.datetime(2021, 5, 14), datetime.datetime(2021, 5, 19), datetime.datetime(2021, 5, 20), datetime.datetime(2021, 5, 22), datetime.datetime(2021, 5, 23), datetime.datetime(2021, 5, 24), datetime.datetime(2021, 5, 29), datetime.datetime(2021, 5, 30), datetime.datetime(2021, 6, 1), datetime.datetime(2021, 6, 2), datetime.datetime(2021, 6, 7), datetime.datetime(2021, 6, 8), datetime.datetime(2021, 6, 11), datetime.datetime(2021, 6, 13), datetime.datetime(2021, 6, 14), datetime.datetime(2021, 6, 25), datetime.datetime(2021, 6, 27), datetime.datetime(2021, 6, 28), datetime.datetime(2021, 6, 30), datetime.datetime(2021, 7, 1), datetime.datetime(2021, 7, 2), datetime.datetime(2021, 7, 3), datetime.datetime(2021, 7, 7), datetime.datetime(2021, 7, 8), datetime.datetime(2021, 7, 9), datetime.datetime(2021, 7, 10), datetime.datetime(2021, 7, 18), datetime.datetime(2021, 7, 23), datetime.datetime(2021, 7, 25), datetime.datetime(2021, 7, 26), datetime.datetime(2021, 7, 29), datetime.datetime(2021, 7, 30), datetime.datetime(2021, 8, 1), datetime.datetime(2021, 8, 2), datetime.datetime(2021, 8, 3), datetime.datetime(2021, 8, 5), datetime.datetime(2021, 8, 6), datetime.datetime(2021, 8, 8), datetime.datetime(2021, 8, 10), datetime.datetime(2021, 8, 11), datetime.datetime(2021, 8, 15), datetime.datetime(2021, 8, 18), datetime.datetime(2021, 8, 19), datetime.datetime(2021, 8, 21), datetime.datetime(2021, 9, 13), datetime.datetime(2021, 9, 23), datetime.datetime(2021, 9, 25), datetime.datetime(2021, 9, 27), datetime.datetime(2021, 9, 30), datetime.datetime(2021, 10, 2), datetime.datetime(2021, 10, 9), datetime.datetime(2021, 10, 10), datetime.datetime(2021, 10, 16), datetime.datetime(2021, 10, 26), datetime.datetime(2021, 10, 30), datetime.datetime(2021, 10, 31), datetime.datetime(2021, 11, 2), datetime.datetime(2021, 11, 6)]

#special_dates = [datetime.datetime(2017, 1, 1), datetime.datetime(2017, 1, 1), datetime.datetime(2017, 1, 1), datetime.datetime(2017, 1, 1), datetime.datetime(2017, 1, 1)]

#_________________________________________________________________________________________________



#-- Loop for each date :
if date_type == "cont" : 
    nb_of_steps = dend - dstart 
    #nb_of_steps = int(nb_of_steps.days / lgt_output)
    nb_of_steps = nb_of_steps.days
    cur_date    = dstart
    if bwd == "T" : cur_date = dend
else :
    nb_of_steps = len(special_dates)

for id_day in range(nb_of_steps):
    print("iteration ", id_day+1, " on ", nb_of_steps)
    
    part_fac_temp = part_fac 
    # in case we specify several part fac, to compute several ones in the loop


    # update current day : 
    if date_type == "special" : 
        cur_date = special_dates[id_day]
    
    # Transform the date into years, month, days :
    cur_year      = str(cur_date.year)
    
    cur_month     = cur_date.month
    if cur_month < 10 :
        cur_month = "0" + str(cur_month)
    else:
        cur_month = str(cur_month)
        
    cur_day       = cur_date.day
    if cur_day < 10 :
        cur_day = "0" + str(cur_day)
    else:
        cur_day = str(cur_day)
        
    cur_date1     = cur_date + datetime.timedelta(1)
    if bwd == "T" : cur_date1 = cur_date - datetime.timedelta(1)
    
    cur_year1     = str(cur_date1.year)

    cur_month1    = cur_date1.month
    if cur_month1 < 10 :
        cur_month1 = "0" + str(cur_month1)
    else:
        cur_month1 = str(cur_month1)
   
    cur_day1      = cur_date1.day
    if cur_day1 < 10 :
        cur_day1 = "0" + str(cur_day1)
    else:
        cur_day1 = str(cur_day1)
   
    cur_date_end  = cur_date + datetime.timedelta(lgt_output+1)
    if bwd == "T" : cur_date_end = cur_date - datetime.timedelta(lgt_output)
        
    cur_year_end  = str(cur_date_end.year)
    
    cur_month_end = cur_date_end.month
    if cur_month_end < 10 :
        cur_month_end = "0" + str(cur_month_end)
    else:
        cur_month_end = str(cur_month_end)

    cur_day_end   = cur_date_end.day
    if cur_day_end < 10 :
        cur_day_end = "0" + str(cur_day_end)
    else:
        cur_day_end = str(cur_day_end)
  
   
   # end date to be written on the file name : 
    cur_datefile     = cur_date_end - datetime.timedelta(2)
    cur_year_file  = str(cur_datefile.year)

    cur_month_file = cur_datefile.month
    if cur_month_file < 10 :
        cur_month_file = "0" + str(cur_month_file)
    else:
        cur_month_file = str(cur_month_file)

    cur_day_file   = cur_datefile.day
    if cur_day_file < 10 :
        cur_day_file = "0" + str(cur_day_file)
    else:
        cur_day_file = str(cur_day_file)




    #-- Change the informations in run.list : 
        
     # Read the file and store its content
    with open("run.list.CIOPS-E", 'r') as file:   lines = file.readlines()
    
     # Modify the lines
    lines[17] = " ystart         = "  + cur_year       + "," + "\n"
    lines[18] = " mstart         = "  + cur_month      + "," + "\n"
    lines[19] = " dstart         = "  + cur_day        + "," + "\n"
    lines[21] = " yend           = "  + cur_year_end   + "," + "\n"
    lines[22] = " mend           = "  + cur_month_end  + "," + "\n"
    lines[23] = " dend           = "  + cur_day_end    + "," + "\n"
    
    if time_resol == "hour" : 
        lines[25] = " infile         = '/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/to_transfer/hour_" + cur_year  + "/" + cur_year  + cur_month  + cur_day  + "00_000.nc'" + "," + "\n"
        lines[26] = " infile2        = '/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/to_transfer/hour_" + cur_year  + "/" + cur_year  + cur_month  + cur_day  + "00_001.nc'" + "," + "\n"
        lines[27] = " infile_loc     = '/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/to_transfer/hour_" + cur_year  + "/'"  + "," + "\n"
    else:
        #lines[25] = " infile        = '/home/bourgoa/Disque_20TO_DATA/DATA/CIOP-E_simba/netcdf_pa/" + cur_year + "/" + cur_year + cur_month + cur_day + "00_000.nc'" + "," + "\n"
        lines[25] = " infile         = '/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/to_transfer/" + cur_year + "/" + cur_year + cur_month + cur_day + "00_000.nc'" + "," + "\n"
        #lines[26] = " infile2       = '/home/bourgoa/Disque_20TO_DATA/DATA/CIOP-E_simba/netcdf_pa/" + cur_year1 + "/" + cur_year1 + cur_month1 + cur_day1 + "00_000.nc'" + "," + "\n"
        lines[26] = " infile2        = '/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/to_transfer/" + cur_year1 + "/" + cur_year1 + cur_month1 + cur_day1 + "00_000.nc'" + "," + "\n"
        lines[27] = " infile_loc     = '/home/ulaval.ca/anbou863/projects/def-frmap5/MINDZ_andeol/to_transfer/" + cur_year1 + "/'"  + "," + "\n"

    lines[62] = " part_spec     = "   + num_tax             + "," + "\n"    
    lines[65] = " in_diapause   = "   + diapause            + "," + "\n"
    lines[70] = " part_fac      = "   + str(part_fac_temp)  + "," + "\n"
    lines[76] = " output_freq   = "   + str(lgt_output*12)  + "," + "\n"
    lines[80] = " y_Astart      = "   + cur_year            + "," + "\n"
    lines[81] = " m_Astart      = "   + cur_month           + "," + "\n"
    lines[82] = " d_Astart      = "   + cur_day             + "," + "\n"
    lines[84] = " y_Aend        = "   + cur_year_end        + "," + "\n"
    lines[85] = " m_Aend        = "   + cur_month_end       + "," + "\n"
    lines[86] = " d_Aend        = "   + cur_day_end         + "," + "\n"
    
    if time_resol == "hour" : 
        lines[91] = " outfile       =  '../outputs/" + taxon + diapause + "_h_" + cur_year_file + "_" + cur_month_file + "_" + cur_day_file + "'," + "\n" 
    
    else :
        lines[91] = " outfile       =  '../outputs/" + taxon + diapause + "_" + cur_year_file + "_" + cur_month_file + "_" + cur_day_file + "'," + "\n" 
        #lines[91] = " outfile       =  '../outputs/" + taxon + diapause + "_5j" +  "'," + "\n"
    if bwd == "T" : lines[91] = " outfile       =  '../outputs/CIOPSE_bwd_" + taxon + diapause + "_"  + cur_year + "_" + cur_month + "_" + cur_day + "'," + "\n"
    
     # Write the modified content back to the file
    with open("run.list.CIOPS-E", 'w') as file:  file.writelines(lines)

    # run the model : 
    subprocess.run("./model_execution.x", shell=True)
    
    # update current day : 
    if date_type == "cont" : 
        cur_date = cur_date + datetime.timedelta(1) #cur_date_end


    
# run the model : 
#subprocess.run("Rscript ../output_analyses/eularian_grids.R", shell=True)
  
    
    
    
    





