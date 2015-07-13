#!/usr/bin/env python
# For a given set of files (given in command line), read them all line by line,
# extract all fields and output them unmodified, except for one field for
# which only the median across the files is reported.
# So for example, if we have 5 files with fields X, Y, Z, and field Z is
# selected (by its index starting in zero, as the first argument),
# Then each row of output will consist of a value of X (must be identical
# across all files in the same row), a value of Y (ditto), and the median value
# of Z across the 5 files.
# Example call: median_runoff.py 5 tmp_*1991.csv > median_excess_runoff.1991.csv

import csv, sys

def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[((len(lst) + 1) / 2) - 1]
    else:
            return float(sum(lst[(len(lst) / 2) - 1:(len(lst) / 2) + 1])) / 2.0


index = int(sys.argv[1])
files = [csv.reader(open(fn, 'r')) for fn in sys.argv[2:]]
lines = [reader.next() for reader in files]  # Dump header lines
print ",".join(lines[0])  # Output CSV header identical to the original files

while True:
    try:
        lines = [reader.next() for reader in files]
    except csv.Error:
        print "Error reading CSV inputs"
    except StopIteration:
        break

    # Output fields that precede index:
    assert all(lines[x][0] == lines[0][0] for x in range(1, index))
    for i in range(0, index):
        sys.stdout.write('"' + lines[0][i] + '"')
        sys.stdout.write(',')

    data = [line[index] for line in lines]
    sys.stdout.write(median(data))

    # Output fields that follow index:
    assert all(lines[x][0] == lines[0][index + 1] for x in range(index + 2, len(lines[0])))
    for i in range(index + 1, len(lines[0])):
        sys.stdout.write(',')
        sys.stdout.write(lines[0][i])

    sys.stdout.write('\n')
