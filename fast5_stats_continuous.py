# -*- coding: utf-8 -*-
"""
Created on Sat Jun 24 15:46:02 2017

@author: Mark
"""
##################################################################################################################################
"""
This program takes basecalled fast5 files as an input and produces skips and stays per base and event, total events, average read length, and 
average current (in picoamps) for 30-second time chunks and outputs to a .csv file.
"""
##################################################################################################################################

##################################################################################################################################
                                                             
##################################################################################################################################

from glob import glob
import numpy as np
import h5py
import math 
import time
import csv
import operator
##################################################################################################################################
                                                            
##################################################################################################################################

CHUNK_SIZE = 30

def display_percent(chunk_size, chunk_percent, last_percent, progress):
    
    """
    Used to monitor progress of a process. Example useage:
        
        Progress = 0
        chunk_percent = 10.0
        chunk_size = int(math.ceil(all_files*(chunk_percent/100)))
        
        for x in all_files:
            Progress += 1 
            last_percent = display_percent(chunk_size, chunk_percent, last_percent, Progress)
    """
    
    percent = int(progress / chunk_size)
    
    if percent > last_percent:
        print("{0}%".format(percent * chunk_percent))
        
    return percent

def colease(n):
    return 1 if n == 0 else n


def process_chunk(files, start_time):

    TotalCurrent = 0
    TotalEvents = 0
    TotalStays = 0
    TotalSkips = 0
    bases_called = 0            
    event_details = {}

    for f in files:
        with h5py.File(f[0], 'r') as hdf:
            Analyses = hdf.get('Analyses')
            Basecall_1D_000 = Analyses.get('Basecall_1D_000')
            BaseCalled_template = Basecall_1D_000.get('BaseCalled_template')
            Events = BaseCalled_template.get('Events')
            
            Current = list(Events['mean'])
            SkipStay = list(Events['move'])

            TotalCurrent += sum(Current)
            TotalEvents += len(Current)          

            Skips = len([1 for x in SkipStay if x > 1])
            Stays = SkipStay.count(0)

            TotalStays += Stays
            TotalSkips += Skips

            bases_called += 4 + sum(SkipStay)



    event_details["Average Current"] = TotalCurrent/colease(TotalEvents)
    event_details["Average Bases Called"] = bases_called/colease(len(files))
    event_details["Skips Per Base"] = float(TotalSkips)/float(colease(bases_called))
    event_details["Stays Per Base"] = float(TotalStays)/float(colease(bases_called))
    event_details["Skips Per Event"] = float(TotalSkips)/float(colease(TotalEvents))
    event_details["Stays Per Event"] = float(TotalStays)/float(colease(TotalEvents))
    event_details["Total Events"] = TotalEvents
    event_details["Total Files"] = len(files)
    event_details["Start Time"] = start_time
    return event_details

def get_file_chunks(sorted_dict, chunk_size_in_sec):
    start_chunk = 0
    end_chunk = 0
    chunked_files = []
    total_events = len(sorted_dict)
    time1 = sorted_dict[start_chunk][1]
    search_time = int(time1)

    while start_chunk < total_events - 1:
        search_time += chunk_size_in_sec

        if sorted_dict[start_chunk][1] > search_time:
            continue

        end_chunk = start_chunk + 1

        while end_chunk != total_events - 1 and sorted_dict[end_chunk][1] < search_time:
            end_chunk += 1

        curr_chunk = sorted_dict[start_chunk:(end_chunk + 1)]

        chunked_files.append([curr_chunk, search_time - chunk_size_in_sec])

        start_chunk = end_chunk + 1

    if start_chunk == total_events - 1:
        chunked_files.append([[sorted_dict[start_chunk]], search_time - chunk_size_in_sec])

    return chunked_files

def main():
    print("Searching for fast5 files...")
    Files = glob("./**/*.fast5", recursive=True)
    TotalFiles = len(Files)
    print("Constructing dictionary of file start times...")
    start_times = {}
    files_not_processed = []

    for fast5 in Files:
        try:
            with h5py.File(fast5,'r') as hdf:
                Analyses = hdf.get('Analyses')
                Basecall_1D_000 = Analyses.get('Basecall_1D_000')
                BaseCalled_template = Basecall_1D_000.get('BaseCalled_template')
                Events = BaseCalled_template.get('Events')
                start_times[fast5] = Events["start"][0]

        except Exception as e:
            print("Error processing {}.".format(fast5))
            files_not_processed.append(fast5)

    sorted_start_times = sorted(start_times.items(), key=operator.itemgetter(1))
    print("Finished")
    print("Splitting up files into corrosponding chunks of {} seconds".format(CHUNK_SIZE))
    chunks = get_file_chunks(sorted_start_times, CHUNK_SIZE)
    
    print("Found {} chunks.".format(len(chunks)))

    total_events = 0

    details = []
    i = 0
    with open('output.csv', 'w') as csvfile:
        field_names = ["Average Current", "Average Bases Called", "Skips Per Base", "Stays Per Base", "Skips Per Event", "Stays Per Event", "Total Events", "Total Files", "Start Time"]
        writer = csv.DictWriter(csvfile, fieldnames=field_names)

        writer.writeheader()
        for chunk in chunks:
            print("Processing chunk {0} of {1} with {2} files...".format(i + 1, len(chunks), len(chunk[0])))
            try:
                curr_detail = process_chunk(chunk[0], chunk[1])
                total_events += curr_detail["Total Events"]
                curr_detail["Total Events"] = total_events
                writer.writerow(curr_detail)
            except Exception as e:
                print("Error processing chunk {}".format(i))
            else:
                print('processed')
            i += 1
        print("Finished processing chunks.")

###########################################################################################################################################

###########################################################################################################################################
		
main()