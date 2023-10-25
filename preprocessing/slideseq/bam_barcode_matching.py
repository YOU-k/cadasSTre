"""
This script cleans a BAM file, perform bead barcode matching, and build a whitelist based on barcode information.
"""
import os
from collections import Counter
import pysam
import gzip
import numpy as np
import editdistance

# function for cleaning a BAM file
def clean_bam(bam_fn,bam_out,bc_ref_dict,kmer_dict,bc_tag="XC"):
    matching_status = Counter()
    whitelist = {}  # {barcode:barcode_name}
    bc_cnt = Counter()
    mismatch_cnt = Counter()
    ixxx = 0
    mybc_cnt = 0
    # open the input BAM file
    with pysam.AlignmentFile(bam_fn, "rb") as bamfile:
        # create the output BAM file
        bam_o = pysam.AlignmentFile(bam_out, "wb", template=bamfile)
        # loop through each read in the input BAM file
        for read in bamfile:
            # extract the barcode from the read
            bc = read.get_tag(bc_tag, with_value_type=False)
            bc = bc[:14]  # truncate the barcode for slideseq
            # match the barcode to a reference barcode
            matched_bc, mismatch_num = barcode_matching(bc_ref_dict,kmer_dict,(bc,) )
            mismatch_cnt[mismatch_num] += 1
            # if no matching barcode -> update the matching status
            if len(matched_bc)==0:
                matching_status["barcode_matching_failed"] += 1
            else:
                # if matching barcode is found, update the whitelist, barcode count, and write the read to the output BAM file
                matching_status["barcode_matched"] += 1
                bc = matched_bc[0]
                if bc =="ACCGTGTTTTGCGG":
                    mybc_cnt += 1
                    print(ixxx, mybc_cnt)
                whitelist[bc] = bc
                bc_cnt[bc] += 1
                read.set_tag("BC", bc)
                bam_o.write(read)
            ixxx += 1
            # print progress and matching status
            if ixxx % 2000000 ==0:
                print(ixxx,"processed.   ",matching_status)
                print("\tmismatch   ",mismatch_cnt)
    # print final status
    print(matching_status)
    # return whitelist and barcode count
    return whitelist, bc_cnt

# function for budilding a kmer distribution
def build_kmer_dist(bc_dict,k_size):
    kmer_dict = {}
    no_unique_kmer = 0
    for bc in bc_dict:
        kmer_added = False
        for ix in range(len(bc)-k_size):
            sub_bc = bc[ix:(ix+k_size)]
            if sub_bc in kmer_dict:
                if len(kmer_dict[sub_bc])<10:
                    kmer_dict[sub_bc].append((bc,ix)) 
                    kmer_added = True
            else:
                kmer_dict[sub_bc] = [(bc,ix)] # kmer: [(barcode, kmer position),...]
                kmer_added = True
        if not kmer_added:
            no_unique_kmer += 1
    print(f"{no_unique_kmer} in {len(bc_dict)} barcode without unique kmer and not added in the kmer dict({len(kmer_dict)})")
    return kmer_dict

# function for kmer matching
def kmer_match(bc,km_dict,max_dist):
    k_size = len(next(iter(km_dict)))
    for ix in range(len(bc)-k_size):
        sub_bc = bc[ix:(ix+k_size)]
        if sub_bc in km_dict:
            fz = [(it[0], editdistance.eval(it[0], bc)) for it in km_dict[sub_bc]]
            if len(fz)>1:
                fz.sort(key=lambda x:x[1] )
            if fz[0][1]<=max_dist:
                return fz[0] # (barcode, offset if there are INDEL)
    return ""

# function for sub-barcode matching
def sub_barcode_matching(bc_ref_dict,kmer_dict,bc,max_dist):
    """
    """
    if bc in bc_ref_dict:
        return bc,0
    else:
        res = kmer_match(bc,kmer_dict,max_dist)
        return res

# function for barcode matching
def barcode_matching(bc_ref_dict,kmer_dict,bc,max_dist=1):
    """
    bc_ref_dict: {"bc1":{bc_x:bcx_name,...},{"bc2":{bc_y:bcy_name,...}..}
    kmer_dict: {"bc1":{kmer1:[(bc_x,offset),...],...},...}
    bc: (bc1,bc2,...)
    max_dist: max edit distance allowed
    return:
        matched_bc_name: (bcx_name,bcy_name,...)
        offset: offset of the last barcode, used to extract UMI if UMI is after last barcode
    """
    bc_seg_num = len(bc) # number of barcode segments
    assert(bc_seg_num==len(bc_ref_dict)) # ensure the number of barcode segments matches the reference dirctionary
    matched_bc_name = [] # list to store matched barcode names
    mismatch = 0 # initialize the mismatch counter
    for ix in range(bc_seg_num):
        bc_seg = bc[ix]
        res = sub_barcode_matching(bc_ref_dict[f"bc{ix+1}"],kmer_dict[f"bc{ix+1}"],bc_seg,max_dist)
        if len(res)>0:
            matched_bc_name.append(bc_ref_dict[f"bc{ix+1}"][res[0]])
            mismatch = res[1]
        else:
            return (),0  # one barcode cannot be matched, return empty tuple
    return matched_bc_name,mismatch


# define main function
def main(fq_dir,bam_fn,bam_out):
    barcode_csv = "/path/to/your/barcode.tsv"
    all_bc = [it.strip().split("\t")[0] for it in open(barcode_csv).readlines()[1:]]
    bc_ref_dict = {"bc1":{}} # create a reference dictionary for barcode segment 1
    for bc in all_bc:
        bc_ref_dict["bc1"][bc] = bc # populate the reference dictionary with barcode data
    print(f"Read {len(bc_ref_dict['bc1'])} reference barcode(bc1).")
    kmer_dict = {}
    kmer_dict["bc1"] = build_kmer_dist(bc_ref_dict["bc1"],7) # build the kmer distribution for barcode segment 1
    whitelist,bc_cnt = clean_bam(bam_fn,bam_out,bc_ref_dict,kmer_dict) # clean the BAM file and generate a whitelist
    whitelist_file = os.path.join(fq_dir,f"included_bc.csv") # define the path for the whitelist file
    print(f"write {len(whitelist)} barcode appeared in the fastq to {whitelist_file}") # get number of barcodes in the whitelist
    with open(whitelist_file,"w") as f:
        f.write("barcode_name,barcode,Num_of_reads\n") # header of output
        for bc in whitelist:
            f.write(f"{whitelist[bc]},{bc},{bc_cnt[bc]}\n")
    return True
    
# Check if the script is run as the main program
if __name__ == "__main__":
    fq_dir = "path/to/your/data"  # Define the directory containing your data
    bam_fn = os.path.join(fq_dir, "sample.bam")  # Define the path to the input BAM file
    bam_out = os.path.join(fq_dir, "sample_clean.bam")  # Define the path to the output BAM file
    main(fq_dir, bam_fn, bam_out)  # Call the main function with specified file paths