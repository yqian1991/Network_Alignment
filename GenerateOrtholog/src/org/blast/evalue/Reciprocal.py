'''
Created on May 21, 2014

@author: yqian33
'''
#! /usr/bin/env python
import sys, csv

# number below which to discard results
E_CUTOFF = 1e-3

def load_csv_to_dict(filename):
    """
    Load the CSV file into a dictionary, trying query sequences to subject
    sequences.
    """
    d = {}

    for (query_name, subject_name, score,expect) in csv.reader(open(filename)):
        # fix the names that start with a 'gi|XXXXXXX|'
        query_name = demangle_name(query_name)
        subject_name = demangle_name(subject_name)

        # convert the e-value into a floating point number
        expect = float(expect)
        if query_name not in d and expect < E_CUTOFF:
            # if we don't already have an entry for this, put it in.
            d[query_name] = subject_name

    return d

def demangle_name(name):
    """
    This functions strips off the 'gi|XXXXX|' name encoding that NCBI uses.

    Note that NCBI does this automatically for subject sequences.
    """
    if name.startswith('gi|'):
        name = name.split('|')
        name = name[2:]
        name = "|".join(name)

    return name

def load_blast_to_dict(infile1):
    d= {}
    
    ins = open(infile1, 'r')
    
    for line in ins:
        items = line.split('\t')
        d[items[0]] = items[1]
    
    return d
        
###

# This is the code that's run when you actually type 'find-reciprocal.py'
# at the command line; the above are function definitions, that define
# reusable blocks or chunks of code.

#in_file_1 = sys.argv[1]
#in_file_2 = sys.argv[2]
def run():
    in_file_1 = '/home/yqian33/APPI/dataset/E-value/clean_flyyeast_uniprot_blast1e-7.txt'
    in_file_2 = '/home/yqian33/APPI/dataset/E-value/clean_yeastfly_uniprot_blast1e-7.txt'
    outfile='/home/yqian33/APPI/dataset/E-value/fly_yeast_reciprocal_test.txt'
    
    d1 = load_blast_to_dict(in_file_1)
    d2 = load_blast_to_dict(in_file_2)
    print 'load done'
    
    #output = csv.writer(sys.stdout)
    output=open(outfile, 'w')
    
    count=0
    for seqname in d1:
        count+=1
        print d1[seqname]
    
    print count
    
    count=0
    for seqname in d1:
        seqmatch1 = d1[seqname]
        seqmatch2 = d2.get(seqmatch1)
        #print seqmatch2
    
        if seqmatch2 == seqname:
            #output.writerow([seqname, seqmatch1])
            count +=1
            print count
            output.write(seqname+'\t'+seqmatch1)
        else:
            print '0'
    
    output.close()

if __name__ == '__main__':
    run()
    pass