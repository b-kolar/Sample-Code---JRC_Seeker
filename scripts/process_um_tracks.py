# Imports
import pandas as pd
import numpy as np
import time
import os
import os.path
from os import path
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-m", dest="meth_track")
parser.add_argument("-u", dest="unmeth_track")
parser.add_argument("-o", dest="output")
args = vars(parser.parse_args())
meth_track = args["meth_track"]
unmeth_track = args["unmeth_track"]
output=args["output"]

m = pd.read_csv(open(meth_track), sep='\t', names=["seg_chr", "segment_start", "segment_end", "state", "cpg_chr", "cpg_start", "cpg_end", "meth_count"])
um = pd.read_csv(open(unmeth_track), sep='\t', names=["seg_chr", "segment_start", "segment_end", "state", "cpg_chr", "cpg_start", "cpg_end", "unmeth_count"])
m.drop(m.columns[[4,5,6]],axis=1,inplace=True)
um.drop(um.columns[[4,5,6]],axis=1,inplace=True)
m.loc[m["meth_count"] == ".", "meth_count"] = 0
um.loc[um["unmeth_count"] == ".", "unmeth_count"] = 0
m['meth_count']=m['meth_count'].astype(int)
um['unmeth_count']=um['unmeth_count'].astype(int)
m = m.groupby(["seg_chr", "segment_start", "segment_end", "state"], sort=False).meth_count.sum().reset_index()
um = um.groupby(["seg_chr", "segment_start", "segment_end", "state"], sort=False).unmeth_count.sum().reset_index()
m["unmeth_count"] = um["unmeth_count"]
m.to_csv(output, sep="\t", index=False, header=False)


