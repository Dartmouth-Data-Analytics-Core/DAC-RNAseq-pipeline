#! /usr/bin/env python

### Method for finding RPKM and TPM courtesy of:
### https://btep.ccr.cancer.gov/question/faq/what-is-the-difference-between-rpkm-fpkm-and-tpm/
### https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/


####### 
## Method for RPKM
#######

### Let million_scaled_read_count = total_read_count (for sample) / 10**6

### For every individual read (like a "fragment", especially for single end),
###     indiv_read_per_million_reads = indiv_read_count / million_scaled_read_count

### Finally, RPKM takes into account the size of the individual read:
###     indiv_read_per_kilobase_million = indiv_read_per_million_reads / (num_bases_in_indiv_read_gene / 10**3)

### This is technically equivalent to
###     indiv_read_per_kilobase_million = indiv_read_count * 10**9 / (total_read_count * num_bases_in_indiv_read_gene)

#### This assumes that the fragments that are good enough quality to be mapped are proportional to the reads per gene
#### Therefore, TPKM follows the same formula

####### 
## Method for TPM
#######

### For each individual read, let:
###     indiv_read_per_kilobase = indiv_read_count / (read_bases_count / 10**3)

### Then, we sum every individual read_per_kilobase value and divide by 10**6:
###     sum_read_per_base_rate_per_million = sum_over_reads{indiv_read_per_kilobase} / 10**6

### Finally, we divide an individual read's per kilobase value by the sum over one million:
###     indiv_transcript_per_million = indiv_read_per_kilobase / sum_read_per_base_rate_per_million

### It turns out that the 10**3's cancel out, so it is technically equivalent to:
###     indiv_transcript_per_million = (indiv_read_count / read_bases_count * 10**6) / (sum_over_reads{indiv_read_count / read_bases_count})

##### Input and Output
##### sys[1] is the tsv file name to take in and normalize, should be formatted in same way as featurecount's 


##### sys[2] is whether the data is single or paired end. Uses RPKM or TPKM and stores to correct accordingl
import sys
import numpy as np
import pandas as pd



def to_rpkm(arr):
    tot_rc = np.sum(arr[1:,6:], axis = 0)
    return 10**9 * arr[1:,6:] * np.outer((1 / arr[1:,5]), (1 / tot_rc))

def to_tpm(arr):
    ind_read_per_base = arr[1:,6:] / arr[1:,5].reshape(arr[1:,6:].shape()[0], 1)

    return (10**6 * ind_read_per_base / np.sum(ind_read_per_base, axis = 0).reshape(ind_read_per_base.shape()[0], 1))

## Basically, where all of the labels are gone (so arr is actually a number array)
def to_rpkm_pre_sliced(arr):
    tot_rc = np.sum(arr[:,1:], axis = 0)
    return arr[:,1:] / np.outer(arr[:,0], tot_rc) * 10**9

def to_tpm_pre_sliced(arr):
    ind_read_per_base = arr[:,1:] / arr[:,0].reshape(arr[:,1:].shape[0], 1)
    return (10**6 * ind_read_per_base / np.sum(ind_read_per_base, axis = 0).reshape(1, ind_read_per_base.shape[1]))



###########

## Extract data

data = pd.read_csv(sys.argv[1], sep='\t')

path_name = sys.argv[1]
path_name_rpkm = (sys.argv[1])[:-4]+"_rpkm.tsv" if sys.argv[2] == "single" else (sys.argv[1])[:-4]+"_fpkm.tsv"
# path_name_fpkm = (sys.argv[1])[:-4]+"_fpkm.tsv"
path_name_tpm = (sys.argv[1])[:-4]+"_tpm.tsv"




# Can use df.loc[a:b,c:d] to slice a dataframe (or maybe iloc works if loc doesn't)
# Can also use df.to_numpy() but remember that numpy arrays that don't hold ints are much slower to create or operate

## Take a slice of the data to only include length and counts


trimmed = np.array(data.iloc[:,5:].to_numpy(),dtype="float32")

## Call the RPKM and TPM functions

rpkm_np = to_rpkm_pre_sliced(trimmed)
tpm_np = to_tpm_pre_sliced(trimmed)

## Convert back to Pandas for tsv file conversion



data_only_rpkm = pd.DataFrame(data = rpkm_np, index = data.index, columns = list(data.columns)[6:])
data_only_tpm = pd.DataFrame(data = tpm_np, index = data.index, columns = list(data.columns)[6:])

rpkm_df = pd.concat([data.iloc[:,:6], data_only_rpkm], axis = 1)
tpm_df = pd.concat([data.iloc[:,:6], data_only_tpm], axis = 1)

rpkm_df.to_csv(path_name_rpkm, sep="\t", index=False)
tpm_df.to_csv(path_name_tpm, sep="\t", index=False)

print("Finished RPKM/FPKM and TPM Normalization")
