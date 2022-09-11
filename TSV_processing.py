#!/usr/bin/env python

## TSV processing
## output format blast6: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def TSV_proc():

    data=pd.read_csv('alignment.b6',sep = '\t', header = None)

    ## set list with a column names for blast output type 6 (-outfmt 6)
    Columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    data.columns = Columns ## add column names

    ## prepare dataset for ploting a histogram
    ## if there is more than one aligment (NC_[number]) for a given query (qseqid, read) choose one aligment - better one
    ## in the final table will be one aligment per each read (I hope I understood the task correctly ...)
    ## I will choose aligment by evalue which is a measure of getting a result by chance - lower evalue means that getting
    ## this result by chance is less likely
    
        
    DF = data.sort_values('evalue', ascending=True) ## sort by descending evalue
    DF_1 = DF.drop_duplicates(subset=['qseqid'],keep='first') ## remove when duplicate values in a query (qseqid) and keep first row (the one with lower evalue)
    DF_2 = DF_1.sort_index(axis = 0) ## sort dataframe by index
    
    
    
    ## check if there are all reads in a filtered dataframe
    x = data.qseqid.unique() ## create an numpy array with unique values for each read from an input dataframe
    y = DF_2.qseqid.unique() ## create an numpy array with unique values for each read from a filtered dataframe

    ## compare arrays to check if each read from the input dataframe is on filtered dataframe
    comparison = x == y
    equal_arrays = comparison.all()
    print(equal_arrays)

    ## plot histogram

    bins = DF_2['length'].nunique() ## set number of bins equal to number of unique length values

    plt.hist(DF_2['length'], bins=bins) ## I setted bins to 100 because histogram doesn't really show exact values otherwise 
    ## (ex. it shows that for 91 value frequency is more that 300 which is not true..)
    plt.ylabel('Frequency count')
    plt.xlabel('Aligment length');
    plt.title('A histogram showing a frequency of aligment lenght values')
    #plt.show()
    plt.savefig('Histogram_aligment_length.pdf')  

    ## create a frequency table to save a csv

    DF_freq = DF_2.length.value_counts()
    #type(DF_freq) ## check type of obtained data
    DF_freq = DF_freq.to_frame().reset_index() ## change pandas_core_series to a dataframe
    column_names = ['aligment_length', 'abundance'] ## create a variable with new names for a table
    DF_freq.columns = column_names ## set new names in a frequency table
    DF_freq = DF_freq.sort_values('aligment_length',ascending=False)
    DF_freq = DF_freq.reset_index(drop = True) ## reset index beacue it looks confusing somehow

    #print(DF_freq)
    DF_freq.to_csv('Aligment_lenght_abundance.csv', index=False, sep = ',') ## save the dataframe to csv file 

if __name__ == "__main__":
    TSV_proc()


