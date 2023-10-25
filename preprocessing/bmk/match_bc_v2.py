"""
### This script processes paired-end sequencing data to identify and match bead barcodes based on specified linkers and starting positions ###
"""
import os
from collections import Counter
import mappy as mp
import editdistance
import gzip
import numpy as np

### trim_bc: trim barcodes from sequencing reads ###
    # fn: Barcode reference file
    # fn1: Input paired-end sequencing file (R2)
    # fn1_out: Output file for trimmed sequencing reads
    # bc_ref_dict: Dictionary containing barcode references
    # kmer_dict: Dictionary containing k-mers for barcode matching
    # linker1, linker2: Linker sequences
    # linker1_start_range: Possible starting positions of linker1
def trim_bc(fn,fn1,fn1_out,bc_ref_dict,kmer_dict,linker1 =  "TACGT",linker2 = "CGACTC",linker1_start_range=[17,18,19,20,21]):
    fn1_out = open(fn1_out,"w") # open output file to write trimmed reads
    ixx = 0
    print(f"{fn}\n{fn1}\n{fn1_out}")
    matching_status = Counter() # counter to keep track of barcode matching status
    whitelist = {} # dictionary to store matched barcodes
    ixxx = 0
    # iterate through paired-end reads using Mappy's fastx_read
    for rec,rec2 in zip(mp.fastx_read(fn, read_comment=False),mp.fastx_read(fn1, read_comment=False)):
        for ix in linker1_start_range:
            # extract three parts of the barcode based on linker positions
            if rec[1][ix:(ix+len(linker1))] == linker1:
                #17+5+19+6+19
                bc1 = rec[1][(ix-17):ix]
                bc2 = rec[1][(ix+len(linker1)):(ix+len(linker1)+19)]
                bc3 = rec[1][(ix+len(linker1)+25):(ix+len(linker1)+44)]
                # call the barcode_matching function to match the barcode
                bc, offset = barcode_matching(bc_ref_dict,kmer_dict,(bc1, bc2, bc3))
                if len(bc)==0:
                    matching_status["barcode_matching_failed"] += 1
                else:
                    matching_status["barcode_matched"] += 1
                    UMI= rec[1][(ix+len(linker1)+44+offset):(ix+len(linker1)+56+offset)]
                    bc_name = f"{bc_ref_dict['bc1'][bc[0]]}-{bc_ref_dict['bc2'][bc[1]]}-{bc_ref_dict['bc3'][bc[2]]}"
                    bc = bc[0]+bc[1]+bc[2]
                    whitelist[bc] = bc_name
                    # write the trimmed read to the output file
                    fn1_out.write(f"@{bc}_{UMI}#{rec2[0]}\n{rec2[1]}\n+\n{rec2[2]}\n")
                break
        ixxx += 1
        # print progress every 10 million reads
        if ixxx % 10000000 ==0:
            print(ixxx,"processed.   ",matching_status)
    print(matching_status)
    return whitelist

### build_8mer_dist: build a dictionary of 8-mer distributions for each barcode ###
    # bc_dict: Dictionary containing barcodes
def build_8mer_dist(bc_dict):
    kmer_dict = {}
    no_unique_kmer = 0
    # iterate through barcodes in the dictionary
    for bc in bc_dict:
        kmer_added = False
        for ix in range(len(bc)-8):
            sub_bc = bc[ix:(ix+8)]
            # check if the 8-mer is in the kmer_dict
            if sub_bc in kmer_dict:
                if len(kmer_dict[sub_bc])<2:
                    kmer_dict[sub_bc].append((bc,ix)) 
                    kmer_added = True
            else:
                kmer_dict[sub_bc] = [(bc,ix)] # kmer: [(barcode, kmer position),...]
                kmer_added = True
        if not kmer_added:
            no_unique_kmer += 1
    print(f"{no_unique_kmer} in {len(bc_dict)} barcode without unique kmer and not added in the kmer dict({len(kmer_dict)})")
    return kmer_dict

### barcode_matching: perform barcode matching using reference barcodes and k-mers ###
    # bc_ref_dict: Dictionary containing barcode references
    # kmer_dict: Dictionary containing k-mers for barcode matching
    # bc: Tuple of three barcode sequences (bc1, bc2, bc3)
    # max_dist: Maximum allowed edit distance for matching
def barcode_matching(bc_ref_dict,kmer_dict,bc,max_dist=3):
    def kmer_match(bc,km_dict):
        for ix in range(len(bc)-8):
            sub_bc = bc[ix:(ix+8)]
            if sub_bc in km_dict:
                fz = [(it[0], it[1]-ix ) for it in km_dict[sub_bc]]
                if len(fz)>1:
                    fz.sort(key=lambda x:abs(x[1]) )
                if fz[0][1]<max_dist:
                    return fz[0] # (barcode, offset if there are INDEL)
        return ""
    bc1, bc2, bc3 = bc
    if True:
        match1 = bc1 in bc_ref_dict["bc1"]
        match2 = bc2 in bc_ref_dict["bc2"]
        match3 = bc3 in bc_ref_dict["bc3"]
        if match1 and match2 and match3:
            return (bc1, bc2, bc3),0
        else:
            offset=0
            mbc1,mbc2,mbc3 = "","",""
            if not match1:
                res = kmer_match(bc1,kmer_dict["bc1"])
                if len(res)>0:
                    mbc1, _ = res
            else:
                mbc1 = bc1
            if not match2:
                res = kmer_match(bc2,kmer_dict["bc2"])
                if len(res)>0:
                    mbc2, _ = res
            else:
                mbc2 = bc2
            if not match3:
                res = kmer_match(bc3,kmer_dict["bc3"])
                if len(res)>0:
                    mbc3, offset = res
            else:
                mbc3 = bc3
            if (mbc1 != "") and (mbc2 != "") and (mbc3 != ""):
                return (mbc1, mbc2, mbc3),offset
    return (),0

### Main function for processing paired-end sequencing data ###
    # fq_dir: Directory containing input files
    # fn1, fn2: Input paired-end sequencing files (R1 and R2)
def main(fq_dir,fn1,fn2):
    barcode_fa = "/mnt/disk2/yue_data/data/fastq/bmk/bcs.fa"
    bc_ref_dict = {"bc1":{},"bc2":{},"bc3":{}}
    for rec in mp.fastx_read(barcode_fa, read_comment=False):
        if rec[0][:3]=="bc1":
            bc_ref_dict["bc1"][rec[1]] = rec[0]
        elif rec[0][:3]=="bc2":
            bc_ref_dict["bc2"][rec[1]] = rec[0]
        elif rec[0][:3]=="bc3":
            bc_ref_dict["bc3"][rec[1]] = rec[0]
    print(f"Read {len(bc_ref_dict['bc1'])}-{len(bc_ref_dict['bc2'])}-{len(bc_ref_dict['bc3'])} reference barcode(bc1-bc2-bc3).")
    kmer_dict = {}
    kmer_dict["bc1"] = build_8mer_dist(bc_ref_dict["bc1"])
    kmer_dict["bc2"] = build_8mer_dist(bc_ref_dict["bc2"])
    kmer_dict["bc3"] = build_8mer_dist(bc_ref_dict["bc3"])
    fn1_out = os.path.join(fq_dir,f"combined_matched.fastq")
    whitelist = trim_bc(os.path.join(fq_dir,fn1),
                        os.path.join(fq_dir,fn2),
                        fn1_out,bc_ref_dict,kmer_dict)
    whitelist_file = os.path.join(fq_dir,f"whitelist.csv")
    print(f"write {len(whitelist)} barcode to whitelist")
    with open(whitelist_file,"w") as f:
        for bc in whitelist:
            f.write(f"{whitelist[bc]},{bc}\n")
    return True
    

if __name__=="__main__":
    fq_dir = "/mnt/disk2/yue_data/data/fastq/bmk/34-1"
    fn1 = "34-1_1.clean.fq.gz"
    fn2 = "34-1_2.clean.fq.gz"
    part_id_range = list(range(1,17))
    main(fq_dir,fn1,fn2)


