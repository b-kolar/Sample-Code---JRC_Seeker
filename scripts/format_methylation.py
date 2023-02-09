# Imports
import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="input",
                    help="input FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="output",
                    help="write output", metavar="OUTPUTFILE")
parser.add_argument("-m", dest="m_track")
parser.add_argument("-u", dest="um_track")
args = vars(parser.parse_args())

df = pd.read_csv(args["input"], sep='\t', names=[
                 "chr", "start", "end", "meth_val", "coverage"])
df['meth_count'] = round(df["coverage"]*df["meth_val"])
df['meth_count'] = df["meth_count"].astype(int)
df['unmeth_count'] = df["coverage"] - df['meth_count']
data = df[["chr", "start", "end", "meth_count", "unmeth_count"]].copy()

data.to_csv(args["output"], sep="\t", header=False, index=False)
methylated = data[["chr", "start", "end", "meth_count"]].copy()
unmethylated = data[["chr", "start", "end", "unmeth_count"]].copy()
methylated.to_csv(args["m_track"], sep="\t", header=False, index=False)
unmethylated.to_csv(args["um_track"], sep="\t", header=False, index=False)
