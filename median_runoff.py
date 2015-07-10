#!/usr/bin/env python
# For a given set of files (given in command line), read them all line by line,
# extract the grid coordinates, state, and total excess runoff, and output
# the coordinates, state, and median runoff (across files) for each point.
# Example call: median_runoff.py tmp_*1991.csv > median_excess_runoff.1991.csv
# Assumption: all inputs have the same coordinates in the same order!

import csv, sys

def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[((len(lst) + 1) / 2) - 1]
    else:
            return float(sum(lst[(len(lst) / 2) - 1:(len(lst) / 2) + 1])) / 2.0


files = [csv.reader(open(fn, 'r')) for fn in sys.argv[1:]]
lines = [reader.next() for reader in files]  # Dump header lines
print "LON,LAT,State,median.excess.runoff"

while True:
    try:
        lines = [reader.next() for reader in files]
    except csv.Error:
        print "Error reading CSV inputs"
    except StopIteration:
        break

    runoffs = [line[5] for line in lines]
    print lines[0][0] + "," + \
            lines[0][1] + "," + \
            lines[0][6] + "," + \
            median(runoffs)

