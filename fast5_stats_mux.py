# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 08:27:31 2017

@author: Mark
"""
###############################################################################
                    # Fast5 Analysis Program - 12/4 #
###############################################################################                                     
"""
This program accepts a set of fast5 files of arbitrarily deep folder structure
and produces a set of analyses for MUX 1, MUX 2, MUX 3, MUX 4, and the overall
run. These analyses include skips per base and event, stays per base and event,
average current (in pA) per event, average read length (in base pairs), total 
reads, total bases, and total events. 
"""                                        
###############################################################################
                            # Imports #
###############################################################################

import tqdm
import time
import csv
import operator
import h5py
from glob import glob

###############################################################################
                        # Fast5 Sorting Module #
###############################################################################

t0 = time.clock()

print()
print("Searching for fast5 files...")
files = glob("./**/*.fast5", recursive=True)

total_reads = len(files)

files_not_processed = []
start_success = 0
start_times = {}

for fast5 in files:
    try:
        with h5py.File(fast5,'r') as hdf:
            Analyses = hdf.get('Analyses')
            Basecall_1D_000 = Analyses.get('Basecall_1D_000')
            BaseCalled_template = Basecall_1D_000.get('BaseCalled_template')
            Events = BaseCalled_template.get('Events')
            start_times[fast5] = Events['start'][0]
            start_success += 1
    except Exception as e:
        files_not_processed.append(fast5)

print()        
print("Start times found for ", start_success, "files.")
print(len(files_not_processed), "files not processed.")

print()
print("Sorting reads by start times...")    
sorted_start_times = sorted(start_times.items(), key=operator.itemgetter(1))

print()
print("Sorting files by MUX...")

mux_1_reads = []
mux_2_reads = []
mux_3_reads = []
mux_4_reads = []

total_reads_1 = 0
total_reads_2 = 0
total_reads_3 = 0
total_reads_4 = 0

index = 0

for fast5 in tqdm.tqdm(sorted_start_times):

    if fast5[1] < 28800:
        index += 1 
        mux_1_reads.append(fast5[0])
        total_reads_1 = index
    elif fast5[1] >= 28800 and fast5[1] < 28800*2:
        index += 1
        mux_2_reads.append(fast5[0])
        total_reads_2 = index
    elif fast5[1] >= 28800*2 and fast5[1] < 28800*3:
        index += 1
        mux_3_reads.append(fast5[0])
        total_reads_3 = index
    elif fast5[1] >= 28800*3:
        index += 1
        mux_4_reads.append(fast5[0])
        total_reads_4 = index
    else:
        print("fatal error")
        
print() 

if total_reads_1 != 0:
    print(total_reads_1, "reads found in MUX 1.")
else:
    print("0 reads found in MUX 4.")
if total_reads_2 != 0:
    print(total_reads_2 - total_reads_1, "reads found in MUX 2.")
else:
    print("0 reads found in MUX 2.")
if total_reads_3 != 0:
    print(total_reads_3 - total_reads_2, "reads found in MUX 3.")
else:
    print("0 reads found in MUX 3.")
if total_reads_4 != 0:
    print(total_reads_4 - total_reads_3, "reads found in MUX 4.")
else:
    print("0 reads found in MUX 4.")

t1 = time.clock()
print()
print("Module time:", round(t1 - t0, 2), "seconds")

###############################################################################
                    # Variable Definitions #
############################################################################### 

files_not_processed_1 = []
total_current_1 = 0
total_events_1 = 0
bases_called_1 = 0
total_skips_1 = 0
total_stays_1 = 0
success_1 = 0
num_reads_1 = 0

files_not_processed_2 = []
total_current_2 = 0
total_events_2 = 0
bases_called_2 = 0
total_skips_2 = 0
total_stays_2 = 0
success_2 = 0
num_reads_2 = 0

files_not_processed_3 = []
total_current_3 = 0
total_events_3 = 0
bases_called_3 = 0
total_skips_3 = 0
total_stays_3 = 0
success_3 = 0
num_reads_3 = 0
   
files_not_processed_4 = []
total_current_4 = 0
total_events_4 = 0
bases_called_4 = 0
total_skips_4 = 0
total_stays_4 = 0
success_4 = 0
num_reads_4 = 0

###############################################################################
skips_event_1 = 0      
skips_base_1 = 0
stays_event_1 = 0        
stays_base_1 = 0
ave_current_1 = 0
ave_read_1 = 0
num_reads_1 = 0

skips_event_2 = 0      
skips_base_2 = 0
stays_event_2 = 0        
stays_base_2 = 0
ave_current_2 = 0
ave_read_2 = 0
num_reads_2 = 0

skips_event_3 = 0      
skips_base_3 = 0
stays_event_3 = 0        
stays_base_3 = 0
ave_current_3 = 0
ave_read_3 = 0
num_reads_3 = 0

skips_event_4 = 0      
skips_base_4 = 0
stays_event_4 = 0        
stays_base_4 = 0
ave_current_4 = 0
ave_read_4 = 0
num_reads_4 = 0

###############################################################################
                    # MUX Processing Modules #
###############################################################################

if total_reads_1 != 0:
    
    print()
    print("Processing MUX 1...")
    
    for fast5 in tqdm.tqdm(mux_1_reads):      
        try:
            with h5py.File(fast5,'r') as hdf:
                Analyses = hdf.get('Analyses')
                Basecall_1D_000 = Analyses.get('Basecall_1D_000')
                BaseCalled_template = Basecall_1D_000.get('BaseCalled_template')
                Events = BaseCalled_template.get('Events')
                
                current = list(Events['mean'])
                skipstay = list(Events['move'])
                skips = len([1 for x in skipstay if x > 1])
                stays = skipstay.count(0)
                
                total_current_1 += sum(current)
                total_events_1 += len(current)
                
                total_skips_1 += skips
                total_stays_1 += stays
                
                bases_called_1 += 4 + sum(skipstay)
                
                success_1 += 1
                
        except Exception as e:
            files_not_processed_1.append(fast5) 
    
    skips_event_1 = float(total_skips_1) / float (total_events_1)       
    skips_base_1 = float(total_skips_1) / float (bases_called_1)
    stays_event_1 = float(total_stays_1) / float (total_events_1)        
    stays_base_1 = float(total_stays_1) / float(bases_called_1)
    ave_current_1 = float(total_current_1) / float(total_events_1)
    ave_read_1 = float(bases_called_1) / float(success_1)
    num_reads_1 = len(mux_1_reads)
    
    print()
    print("Skips/Event: ", skips_event_1)
    print("Skips/Base: ", skips_base_1)
    print("Stays/Event: ", stays_event_1)
    print("Stays/Base: ", stays_base_1)
    print("Average Current: ", ave_current_1, " pA")
    print("Average Read Length: ", ave_read_1, " bp")
    print("Total Reads: ", num_reads_1)
    print("Total Bases: ", bases_called_1)
    print("Total Events: ", total_events_1)
    print("Files not processed: ", len(files_not_processed_1))

else:
    pass

t2 = time.clock()
print()
print("Module time:", round(t2 - t1, 2), "seconds")

###############################################################################

if total_reads_2 != 0:

    print()
    print("Processing MUX 2...")
    
    for fast5 in tqdm.tqdm(mux_2_reads):
        try:
            with h5py.File(fast5,'r') as hdf:
                Analyses = hdf.get('Analyses')
                Basecall_1D_000 = Analyses.get('Basecall_1D_000')
                BaseCalled_template = Basecall_1D_000.get('BaseCalled_template')
                Events = BaseCalled_template.get('Events')
                
                current = list(Events['mean'])
                skipstay = list(Events['move'])
                skips = len([1 for x in skipstay if x > 1])
                stays = skipstay.count(0)
                
                total_current_2 += sum(current)
                total_events_2 += len(current)
                
                total_skips_2 += skips
                total_stays_2 += stays
                
                bases_called_2 += 4 + sum(skipstay)
                
                success_2 += 1
                
        except Exception as e:
            files_not_processed_2.append(fast5) 
    
    skips_event_2 = float(total_skips_2) / float (total_events_2)       
    skips_base_2 = float(total_skips_2) / float (bases_called_2)
    stays_event_2 = float(total_stays_2) / float (total_events_2)        
    stays_base_2 = float(total_stays_2) / float(bases_called_2)
    ave_current_2 = float(total_current_2) / float(total_events_2)
    ave_read_2 = float(bases_called_2) / float(success_2)
    num_reads_2 = len(mux_2_reads)
    
    print()
    print("Skips/Event: ", skips_event_2)
    print("Skips/Base: ", skips_base_2)
    print("Stays/Event: ", stays_event_2)
    print("Stays/Base: ", stays_base_2)
    print("Average Current: ", ave_current_2, " pA")
    print("Average Read Length: ", ave_read_2, " bp")
    print("Total Reads: ", num_reads_2)
    print("Total Bases: ", bases_called_2)
    print("Total Events: ", total_events_2)
    print("Files not processed: ", len(files_not_processed_2))
    
else:
    pass

t3 = time.clock()
print()
print("Module time:", round(t3 - t2, 2), "seconds")

###############################################################################

if total_reads_3 != 0:

    print()
    print("Processing MUX 3...")

    for fast5 in tqdm.tqdm(mux_3_reads):
        try:
            with h5py.File(fast5,'r') as hdf:
                Analyses = hdf.get('Analyses')
                Basecall_1D_000 = Analyses.get('Basecall_1D_000')
                BaseCalled_template = Basecall_1D_000.get('BaseCalled_template')
                Events = BaseCalled_template.get('Events')
                
                current = list(Events['mean'])
                skipstay = list(Events['move'])
                skips = len([1 for x in skipstay if x > 1])
                stays = skipstay.count(0)
                
                total_current_3 += sum(current)
                total_events_3 += len(current)
                
                total_skips_3 += skips
                total_stays_3 += stays
                
                bases_called_3 += 4 + sum(skipstay)
                
                success_3 += 1
                
        except Exception as e:
            files_not_processed_3.append(fast5) 
    
    skips_event_3 = float(total_skips_3) / float (total_events_3)       
    skips_base_3 = float(total_skips_3) / float (bases_called_3)
    stays_event_3 = float(total_stays_3) / float (total_events_3)        
    stays_base_3 = float(total_stays_3) / float(bases_called_3)
    ave_current_3 = float(total_current_3) / float(total_events_3)
    ave_read_3 = float(bases_called_3) / float(success_3)
    num_reads_3 = len(mux_3_reads)
    
    print()
    print("Skips/Event: ", skips_event_3)
    print("Skips/Base: ", skips_base_3)
    print("Stays/Event: ", stays_event_3)
    print("Stays/Base: ", stays_base_3)
    print("Average Current: ", ave_current_3, " pA")
    print("Average Read Length: ", ave_read_3, " bp")
    print("Total Reads: ", num_reads_3)
    print("Total Bases: ", bases_called_3)
    print("Total Events: ", total_events_3)
    print("Files not processed: ", len(files_not_processed_3))
    
else:
    pass

t4 = time.clock()
print()
print("Module time:", round(t4 - t3, 2), "seconds")

###############################################################################

if total_reads_4 != 0:

    print()
    print("Processing MUX 4...")
    
    for fast5 in tqdm.tqdm(mux_4_reads):
        try:
            with h5py.File(fast5,'r') as hdf:
                Analyses = hdf.get('Analyses')
                Basecall_1D_000 = Analyses.get('Basecall_1D_000')
                BaseCalled_template = Basecall_1D_000.get('BaseCalled_template')
                Events = BaseCalled_template.get('Events')
                
                current = list(Events['mean'])
                skipstay = list(Events['move'])
                skips = len([1 for x in skipstay if x > 1])
                stays = skipstay.count(0)
                
                total_current_4 += sum(current)
                total_events_4 += len(current)
                
                total_skips_4 += skips
                total_stays_4 += stays
                
                bases_called_4 += 4 + sum(skipstay)
                
                success_4 += 1
                
        except Exception as e:
            files_not_processed_1.append(fast5) 
    
    skips_event_4 = float(total_skips_4) / float (total_events_4)       
    skips_base_4 = float(total_skips_4) / float (bases_called_4)
    stays_event_4 = float(total_stays_4) / float (total_events_4)        
    stays_base_4 = float(total_stays_4) / float(bases_called_4)
    ave_current_4 = float(total_current_4) / float(total_events_4)
    ave_read_4 = float(bases_called_4) / float(success_4)
    num_reads_4 = len(mux_4_reads)
    
    print()
    print("Skips/Event: ", skips_event_4)
    print("Skips/Base: ", skips_base_4)
    print("Stays/Event: ", stays_event_4)
    print("Stays/Base: ", stays_base_4)
    print("Average Current: ", ave_current_4, " pA")
    print("Average Read Length: ", ave_read_4, " bp")
    print("Total Reads: ", num_reads_4)
    print("Total Bases: ", bases_called_4)
    print("Total Events: ", total_events_4)
    print("Files not processed: ", len(files_not_processed_4))
    
else:
    pass

t5 = time.clock()
print()
print("Module time:", round(t5 - t4, 2), "seconds")

###############################################################################
                # Combined MUX Module #
###############################################################################

overall_total_events = total_events_1 + total_events_2 + total_events_3 + total_events_4
overall_total_bases = bases_called_1 + bases_called_2 + bases_called_3 + bases_called_4
overall_num_reads = num_reads_1 + num_reads_2 + num_reads_3 + num_reads_4

overall_skips_event = float(total_skips_1 + total_skips_2 + total_skips_3 + total_skips_4) / float(overall_total_events)
overall_skips_base = float(total_skips_1 + total_skips_2 + total_skips_3 + total_skips_4) / float(overall_total_bases)
overall_stays_event = float(total_stays_1 + total_stays_2 + total_stays_3 + total_stays_4) / float(overall_total_events)
overall_stays_base = float(total_stays_1 + total_stays_2 + total_stays_3 + total_stays_4) / float(overall_total_bases)
overall_ave_current = float(total_current_1 + total_current_2 + total_current_3 + total_current_4) / float(overall_total_events)
overall_ave_read = float(bases_called_1 + bases_called_2 + bases_called_3 + bases_called_4) / float(overall_num_reads)
overall_files_not_processed = len(files_not_processed_1) + len(files_not_processed_2) + len(files_not_processed_3) + len(files_not_processed_4)
overall_reads_processed = success_1 + success_2 + success_3 + success_4

print()
print("Overall Skips/Event: ", overall_skips_event)
print("Overall Skips/Base: ", overall_skips_base)
print("Overall Stays/Event: ", overall_stays_event)
print("Overall Stays/Base: ", overall_stays_base)
print("Overall Average Current: ", overall_ave_current)
print("Overall Average Read Length: ", overall_ave_read)
print("Overall Total Reads: ", overall_num_reads)
print("Overall Total Bases: ", overall_total_bases)
print("Overall Total Events: ", overall_total_events)
print("Overall Reads Processed: ", overall_reads_processed)
print("Overall files not processed: ", overall_files_not_processed)

print()
print()

###############################################################################
                    # .csv output module #
###############################################################################

print("Writing .csv file...")

header1 = ['mux 1 skips/event','mux 1 skips/base','mux 1 stays/event', \
          'mux 1 stays/base','mux 1 ave current','mux 1 ave read length', \
          'mux 1 total reads','mux 1 total bases', 'mux 1 total events', \
          'mux 1 reads processed','mux 1 files not processed']

header2 = ['mux 2 skips/event','mux 2 skips/base','mux 2 stays/event', \
          'mux 2 stays/base','mux 2 ave current','mux 2 ave read length', \
          'mux 2 total reads','mux 2 total bases', 'mux 2 total events', \
          'mux 2 reads processed','mux 2 files not processed']

header3 = ['mux 3 skips/event','mux 3 skips/base','mux 3 stays/event', \
          'mux 3 stays/base','mux 3 ave current','mux 3 ave read length', \
          'mux 3 total reads','mux 3 total bases', 'mux 3 total events', \
          'mux 3 reads processed','mux 3 files not processed']

header4 = ['mux 4 skips/event','mux 4 skips/base','mux 4 stays/event', \
          'mux 4 stays/base','mux 4 ave current','mux 4 ave read length', \
          'mux 4 total reads','mux 4 total bases', 'mux 4 total events', \
          'mux 4 reads processed','mux 4 files not processed']

header = ['overall skips/event','overall skips/base','overall stays/event', \
          'overall stays/base','overall ave current','overall ave read length', \
          'overall total reads','overall total bases', 'overall total events', \
          'overall reads processed','overall files not processed']


values1 = [skips_event_1, skips_base_1, stays_event_1, stays_base_1, \
          ave_current_1, ave_read_1, num_reads_1, bases_called_1, \
          total_events_1, success_1, len(files_not_processed_1)]

values2 = [skips_event_2, skips_base_2, stays_event_2, stays_base_2, \
          ave_current_2, ave_read_2, num_reads_2, bases_called_2, \
          total_events_2, success_2, len(files_not_processed_2)]

values3 = [skips_event_3, skips_base_3, stays_event_3, stays_base_3, \
          ave_current_3, ave_read_3, num_reads_3, bases_called_3, \
          total_events_3, success_3, len(files_not_processed_3)]

values4 = [skips_event_4, skips_base_4, stays_event_4, stays_base_4, \
          ave_current_4, ave_read_4, num_reads_4, bases_called_4, \
          total_events_4, success_4, len(files_not_processed_4)]

values = [overall_skips_event, overall_skips_base, overall_stays_event, overall_stays_base, \
          overall_ave_current, overall_ave_read, overall_num_reads, overall_total_bases, \
          overall_total_events, overall_reads_processed, overall_files_not_processed]

with open('analyses_output.csv', 'w') as csvfile:
    
    writer = csv.writer(csvfile)
    
    try:
        writer.writerow(header1)
        writer.writerow(values1)
    except:
        pass
    writer.writerow([])
    
    try:
        writer.writerow(header2)
        writer.writerow(values2)
    except:
        pass
    writer.writerow([])
    
    try:
        writer.writerow(header3)
        writer.writerow(values3)
    except:
        pass
    writer.writerow([])
    
    try:
        writer.writerow(header4)
        writer.writerow(values4)
    except:
        pass
    writer.writerow([])
    
    try:
        writer.writerow(header)
        writer.writerow(values)
    except:
        pass
    writer.writerow([])
    
    try:
        writer.writerow(['total files found', 'total files not processed'])
        writer.writerow([len(files), len(files_not_processed)])
    except:
        pass
    
print()
print("File saved.")

print()
print("Total processing time:", round(time.clock()/60, 2), "minutes")

print()
print()
print("Program complete.")


    

        
        