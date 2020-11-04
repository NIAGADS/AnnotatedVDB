#!/usr/bin/env python
# pylint: disable=invalid-name,unused-import
"""
Generates and load the index bins
gets db connection info from $GUS_HOME/config/gus.config
"""
from __future__ import print_function

import argparse

from psycopg2.extras import NumericRange

from GenomicsDBData.Util.utils import warning, xstr, die
from GenomicsDBData.Util.postgres_dbi import Database
from collections import OrderedDict

def read_chr_map():
    ''' read chr map file and store as dictionary '''
    result = OrderedDict()
    with open(args.chromosomeMap, 'r') as fh:
        for line in fh:
            chrom, chrLength = line.rstrip().split('\t')
            result[chrom] = int(chrLength)
            
    return result


def load_bins():
    ''' generate and load bins '''
    for chrom, seqLength in chrMap.iteritems():
        warning("Processing", chrom)
        binRoot = chrom
        level = 0
        increments[0] = seqLength
        generate_bins(binRoot, 0, seqLength, level, seqLength)
        
        warning("Done with", chrom)
        if args.commit:
            database.commit()
            warning("Committed")
        else:
            database.rollback()
            warning("Rolled back")


def generate_bins(binRoot, locStart, locEnd, level, seqLength):
    ''' recursive function for generating bin index '''
    global binCount
    if level >= numLevels: return
    
    lowerBound = locStart
    upperBound = locStart + increments[level]

    currentBin= 0

    if locEnd > seqLength: locEnd = seqLength
    
    while lowerBound < locEnd:
        binCount = binCount + 1
        currentBin = currentBin + 1
        binLabel = binRoot + ".B" + xstr(currentBin) if level != 0 else xstr(binRoot)
        if upperBound > seqLength:
            upperBound = seqLength
        if upperBound > locEnd: 
            upperBound = locEnd

        insert_bin(level, binCount, binLabel, lowerBound, upperBound)
        
        nextLevel = level + 1
        if nextLevel <= numLevels:
            if lowerBound == 0:
                warning("New Level:", level)
            
            generate_bins(binLabel + ".L" + xstr(nextLevel), lowerBound, upperBound, nextLevel, seqLength)

        lowerBound = upperBound
        upperBound = upperBound + increments[level]

def insert_bin(level, binId, binPath, locStart, locEnd):
    ''' inserts bin in to db '''
    values = binPath.split('.')
    chrom = values[0]
    cursor.execute(insertSql,( chrom, level, binId, binPath, NumericRange(locStart, locEnd, '(]')))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate and load the BinIndex reference table')
    parser.add_argument('-m', '--chromosomeMap', help="full path file containing mapping of chr names to length; tab-delim, no header", required=True)
    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile', '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    args = parser.parse_args()

    increments = [-1, 64000000, 32000000, 16000000, 8000000, 4000000, 2000000, 1000000, 500000, 250000, 125000, 62500, 31250, 15625]
    binCount = 0
    numLevels = len(increments)

    chrMap = read_chr_map()
    insertSql = "INSERT INTO BinIndexRef (chromosome, level, global_bin, global_bin_path, location) VALUES (%s, %s, %s, %s, %s)"

    database = Database(args.gusConfigFile)
    database.connect()

    cursor = database.cursor()
    load_bins()

    database.close()
