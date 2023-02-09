# Import
import pandas as pd
import numpy as np
import time
import os
import os.path
from os import path
from argparse import ArgumentParser

def bin_polish(segment_file, cpg_intersect_file, mappability_file, meth_unmeth_intersect_file, blacklist_file, output_path, region_output_path, im_state, nd_state, m_state, u_state, map_threshold, segment_min_sz, bin_sz):

    ############################################
    # READ IN CHROMHMM SEGMENTATION FILE
    df = pd.read_csv(open(segment_file), sep='\t', names=["chr", "segment_start", "segment_end", "state"])
    df['segment_length'] = df['segment_end'] - df['segment_start']
    data_length = df['segment_length'].sum()
    
    ###########################################
    # IMPORT HELPER COLUMNS
    # CpG density for original 200 bp segmentation
    cpg_data = pd.read_csv(open(cpg_intersect_file), sep='\t', names=["chr", "start", "end", "state", "cpg_density"])
    df["cpg_density"] = cpg_data['cpg_density'].values
    # Add mappability coverage column
    map_data = pd.read_csv(open(mappability_file), sep='\t', names=["chr",  "segment_start", "segment_end", "state", "map_coverage"])
    df = pd.merge(df, map_data, on=["chr", "segment_start", "segment_end", "state"])
    # Add methylated + unmethylated count columns
    meth_data = pd.read_csv(open(meth_unmeth_intersect_file), sep='\t', names=["chr",  "segment_start", "segment_end", "state", "meth_count", "unmeth_count"])
    df["meth_count"] = meth_data["meth_count"].values
    df["unmeth_count"] = meth_data["unmeth_count"].values
    # Add blacklist count column
    blacklist = pd.read_csv(open(blacklist_file), sep='\t', names=["chr", "start", "end", "state", "blacklist_overlap_count"])
    df["blacklist_count"] = blacklist["blacklist_overlap_count"].values
    
    ###########################################
    # CALCULATE INITIAL LIKELIHOOD
    
    liklihood_df = df.copy()
    initial_sse = calc_likelihood(liklihood_df)
    print("Initial SSE: ", initial_sse)
    
    ###########################################
    # REMOVE DATA
    
    og_chrom = df["chr"].unique()
    # States with low mappability and blacklist region overlap -> no_data states
    old_len = len(df)
    df.loc[(df['map_coverage'] < map_threshold) | (df['blacklist_count'] > 0) , 'state'] = nd_state 
    df = df.drop(columns=['map_coverage', 'blacklist_count'], axis=1)
    df = clean_df(df)
    print("Removed ", old_len - len(df), " regions due to low mappability + blacklist region overlap.")
    num_im_regions(df.copy(), im_state, og_chrom)
   
    ###########################################
    # POLISHING ALGORITHM
    
    medium_thresh = 2*bin_sz
    larger_thresh = 3*bin_sz

    # 1 - If 2 IM regions are separated by 1 200 bp bin, merge these regions
    df.loc[(df['state'] != im_state) & (df['state'] != nd_state) &(df["state"].shift(-1) == im_state) & (df["state"].shift(1) == im_state) & (df['segment_length'] == bin_sz) & (df["chr"].shift(1) == df['chr']) & (df["chr"].shift(-1) == df['chr']), 'state'] = im_state
    df = clean_df(df)
    num_im_regions(df.copy(), im_state, og_chrom)
    
    liklihood_df = df.copy()
    temp = calc_likelihood(liklihood_df)
    print("1 SSE: ", temp)  
   
    # 2 - If small IM region is within two large other regions => classify as earlier bin state
    df.loc[(df['state'] == im_state) & (df["state"].shift(1) != im_state) & (df['segment_length'].shift(-1) >= larger_thresh) & (df["state"].shift(-1) != im_state) & (df['segment_length'].shift(+1) >= larger_thresh) & (df['segment_length'] == bin_sz) & (df["chr"].shift(1) == df['chr']) & (df["chr"].shift(-1) == df['chr']) & (df["segment_length"].shift(+1) >= df["segment_length"].shift(-1)), 'state'] = df['state'].shift(+1)
    df.loc[(df['state'] == im_state) & (df["state"].shift(1) != im_state) & (df['segment_length'].shift(-1) >= larger_thresh) & (df["state"].shift(-1) != im_state) & (df['segment_length'].shift(+1) >= larger_thresh) & (df['segment_length'] == bin_sz) & (df["chr"].shift(1) == df['chr']) & (df["chr"].shift(-1) == df['chr']) & (df["segment_length"].shift(-1) > df["segment_length"].shift(+1)), 'state'] = df['state'].shift(-1)
    df = clean_df(df)
    num_im_regions(df.copy(), im_state, og_chrom)    
    liklihood_df = df.copy()
    temp = calc_likelihood(liklihood_df)
    print("2 SSE: ", temp)    

    # 3 - If 2 IM regions are separated by 1 400 bp bin that has less than 2 cpgs, merge these IM regions
    df.loc[(df['cpg_density'] <= 1) & (df['state'] != nd_state) & (df['state'] != im_state) & (df["state"].shift(+1) == im_state) & (df["state"].shift(-1) == im_state) & (df["chr"].shift(1) == df['chr']) & (df["chr"].shift(-1) == df['chr']), 'state'] = im_state
    df = clean_df(df)    
    liklihood_df = df.copy()
    temp = calc_likelihood(liklihood_df)
    print("3 SSE: ", temp)
    
    # 3b - If IM region is classified wrongly - re-classify
    df.loc[(df['state'] == im_state) & (df['meth_count'] > 0) & (df['unmeth_count'] == 0), 'state'] = m_state
    df.loc[(df['state'] == im_state) & (df['unmeth_count'] > 0) & (df['meth_count'] == 0), 'state'] = u_state    
    liklihood_df = df.copy()
    temp = calc_likelihood(liklihood_df)
    print("3b SSE: ", temp)
    
    # 4 - If 2 IM regions are separated by 1 200 bp bin, merge these regions
    df.loc[(df['state'] != nd_state) &(df['state'] != im_state) & (df["state"].shift(-1) == im_state) & (df["state"].shift(1) == im_state) & (df['segment_length'] == bin_sz) & (df["chr"].shift(1) == df['chr']) & (df["chr"].shift(-1) == df['chr']), 'state'] = im_state
    df = clean_df(df)
    num_im_regions(df.copy(), im_state, og_chrom)    
    liklihood_df = df.copy()
    temp = calc_likelihood(liklihood_df)
    print("4 SSE: ", temp)
    
    # 4b TRY - if <=600 state is within two IM states that have higher cpg density, merge
    df.loc[(df['state'] != im_state) & (df['state'] != nd_state) & (df["state"].shift(-1) == im_state) & (df["state"].shift(1) == im_state) & (df['segment_length'] <= larger_thresh) & (df["chr"].shift(1) == df['chr']) & (df["chr"].shift(-1) == df['chr'])& (df["cpg_density"].shift(+1) >= df['cpg_density'])& (df["cpg_density"].shift(-1) >= df['cpg_density']), 'state'] = im_state
    df = clean_df(df)
    num_im_regions(df.copy(), im_state, og_chrom)    
    liklihood_df = df.copy()
    temp = calc_likelihood(liklihood_df)
    print("4b SSE: ", temp)
    
    # 5 - If any state has low cpg density, make it no data
    df.loc[(df['cpg_density'] <= 2), 'state'] = nd_state
    df = clean_df(df)
    num_im_regions(df.copy(), im_state, og_chrom)    
    liklihood_df = df.copy()
    temp = calc_likelihood(liklihood_df)
    print("5 SSE: ", temp)
    
    # 6 - If small IM region is within two large other regions => classify as earlier bin state
    df.loc[(df['state'] == im_state) &(df['state'].shift(1) == df['state'].shift(-1))& (df["state"].shift(1) != im_state) & (df['segment_length'].shift(-1) >= larger_thresh) & (df["state"].shift(-1) != im_state) & (df['segment_length'].shift(+1) >= larger_thresh) & (df['segment_length'] <=medium_thresh) & (df["chr"].shift(1) == df['chr']) & (df["chr"].shift(-1) == df['chr'])& (df["segment_length"].shift(+1) >= df["segment_length"].shift(-1)), 'state'] = df['state'].shift(+1)
    df.loc[(df['state'] == im_state) &(df['state'].shift(1) == df['state'].shift(-1))& (df["state"].shift(1) != im_state) & (df['segment_length'].shift(-1) >= larger_thresh) & (df["state"].shift(-1) != im_state) & (df['segment_length'].shift(+1) >= larger_thresh) & (df['segment_length'] <=medium_thresh) & (df["chr"].shift(1) == df['chr']) & (df["chr"].shift(-1) == df['chr'])& (df["segment_length"].shift(-1) > df["segment_length"].shift(+1)), 'state'] = df['state'].shift(-1)
    df = clean_df(df)
    print("FINAL 200 metric")
    num_im_regions(df.copy(), im_state, og_chrom)    
    liklihood_df = df.copy()
    temp = calc_likelihood(liklihood_df)
    print("6 SSE: ", temp)
    
    #7 - Remove small 200bp bins
    df.loc[(df['segment_length'] <= segment_min_sz), 'state'] = nd_state
    df = clean_df(df)
    num_im_regions(df.copy(), im_state, og_chrom)
    
    ###########################################
    # CALCULATE FINAL LIKELIHOOD
    
    liklihood_df = df.copy()
    final_sse = calc_likelihood(liklihood_df)
    print("Final SSE: ", final_sse)    

    um_length = df[df["state"] == im_state]
    um_length_ = um_length["segment_length"].sum()
    print(data_length, um_length_)
    print("Percent of the genome removed: ", (100*(data_length - um_length_)/data_length), "%")
    print("Total length of UM segments: ", um_length_)     
    df = df.drop(columns=['segment_length', 'cpg_density',"unmeth_count","meth_count"], axis=1)    
    df.to_csv(output_path, sep="\t", index=False, header=False)   
    df = df.loc[df["state"] == im_state]    
    df["regions"] = df["chr"] + ":" +  df["segment_start"].astype(str) + "-" + df["segment_end"].astype(str)
    df["regions"].to_csv(region_output_path, index=False, header=False)    

def calc_likelihood(df):
    df["methylation_fraction"] = np.nan
    df["meth_sum"] = df["meth_count"] + df["unmeth_count"]
    df.loc[(df["meth_sum"] > 0), "methylation_fraction"] = df["meth_count"]/df["meth_sum"]
    df["state_methylation_fraction"] = np.nan
    df.loc[df["state"] == im_state, 'state_methylation_fraction'] = 0.5
    df.loc[df["state"] == u_state, 'state_methylation_fraction'] = 0.0
    df.loc[df["state"] == m_state, 'state_methylation_fraction'] = 1.0
    l_cal = df.copy()
    l_cal = l_cal[l_cal['state']!=nd_state]
    l_cal = l_cal.dropna().reset_index(drop=True)
    l_cal = l_cal.drop(["chr", "segment_start", "segment_end", "state","cpg_density","segment_length"], axis=1)
    l_cal["sse"] = (l_cal["methylation_fraction"] - l_cal["state_methylation_fraction"])**2
    sse = (1/len(l_cal))*l_cal["sse"].sum()
    return sse

def clean_df(df):    
    # "chr", "segment_start", "segment_end", "state", "segment_length"
    df["next_state"] = df["state"].shift()
    df["cumsum"] = ((df["state"] != df["next_state"])).cumsum()
    df = df.drop(columns=['next_state'])
    df = df.groupby(["cumsum", "chr"]).agg({'segment_start':'min', 'segment_end':'max','state':'first','segment_length':'sum', 'cpg_density':'sum', 'meth_count':'sum', 'unmeth_count':'sum'}).reset_index()
    df = df.drop(columns=['cumsum'])
    return df
    
def num_im_regions(df, im_state, og_chrom):
    # Number of IM states
    df['counts'] = 0
    df.loc[(df['state'] == im_state) & (df['state'].shift(-1) != im_state) & (df['chr'].shift(-1) == df['chr']), 'counts'] = 1
    count = df[df.counts == 1].shape[0]
    df = df.drop(columns=['counts'],axis=1)
    print("IM states found: ", count)

    # Number of small (200bp) IM states
    # make sure chromosomes are the same
    df['counts'] = 0
    df.loc[(df['state'].shift(-1) != im_state) & (df['state'].shift(+1) != im_state) & (df['state'] == im_state) & (df['segment_length'] ==segment_min_sz) & (df['chr'].shift(+1) == df['chr']) & (df['chr'].shift(-1) == df['chr']), 'counts'] = 1
    count = df[df.counts == 1].shape[0]
    df = df.drop(columns=['counts'],axis=1)
    print("Number of 200bp IM states: ", count)

def assign_states(state_labels_file):
    df = pd.read_csv(open(state_labels_file), sep='\t')
    im_state = df.loc[df['label'] == 'IM', 'state'].iloc[0]
    nd_state = df.loc[df['label'] == 'ND', 'state'].iloc[0]
    m_state = df.loc[df['label'] == 'M', 'state'].iloc[0]
    u_state = df.loc[df['label'] == 'U', 'state'].iloc[0]    
    return im_state, nd_state, m_state, u_state

parser = ArgumentParser()

parser.add_argument("-s", dest="seg")
parser.add_argument("-c", dest="cpg")
parser.add_argument("-b", dest="blacklist")
parser.add_argument("-m", dest="map")
parser.add_argument("-t", dest="um_track")
parser.add_argument("-l", dest="labels")
parser.add_argument("-k", dest="map_k")
parser.add_argument("-y", dest="min_bin")
parser.add_argument("-o", dest="output")
parser.add_argument("-r", dest="region_out")
parser.add_argument("-i", dest="bin_sz")

args = vars(parser.parse_args())
    
segment_file = args["seg"]
output_path = args["output"] 
cpg_intersect_file = args["cpg"]
mappability_file = args["map"]
blacklist_file = args["blacklist"]
state_labels_file = args["labels"]
meth_unmeth_intersect_file = args["um_track"] 
map_threshold = float(args["map_k"])
segment_min_sz = args["min_bin"]
region_output_path = args["region_out"]
bin_sz = args["bin_sz"]

im_state, nd_state, m_state, u_state = assign_states(state_labels_file)

bin_polish(segment_file, cpg_intersect_file, mappability_file, meth_unmeth_intersect_file, blacklist_file, output_path, region_output_path, im_state, nd_state, m_state, u_state, map_threshold, segment_min_sz, bin_sz)

