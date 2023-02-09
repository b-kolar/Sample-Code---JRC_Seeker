# Imports
import pandas as pd
import numpy as np
import time
import os.path
from os import path
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", dest="input")
parser.add_argument("-o", dest="output")
args = vars(parser.parse_args())
input_file=args["input"]
output_file=args["output"]

df = pd.read_csv(open(input_file), skiprows=[0], sep='\t', names=["state","methylated", "unmethylated"])
df = df.sort_values(["methylated"], ascending = (False)).reset_index(drop=True)
df["label"] = np.nan
m_index = -1
u_index = -1

if df.at[0,"unmethylated"] >= df.at[1,"unmethylated"]:
    df.at[0,"label"] = "IM"
    df.at[1,"label"] = "M"
    m_index = 1
else:
    df.at[1,"label"] = "IM"
    df.at[0,"label"] = "M"
    m_index = 0

if df.at[2,"unmethylated"] >= df.at[3,"unmethylated"]:
    df.at[2,"label"] = "U"
    df.at[3,"label"] = "ND"
    u_index = 2
else:
    df.at[3,"label"] = "U"
    df.at[2,"label"] = "ND"
    u_index = 3
    
df['state'] = 'E' + df['state'].astype(str)

if (df.at[u_index, "methylated"] > df.at[u_index, "unmethylated"]):
    print("ERROR: UNMETHYLATED STATE NOT CLASSIFIED CORRECTLY. RUN ON MORE DATA.")
 
if (df.at[m_index, "unmethylated"] > df.at[m_index, "methylated"]):
    print("ERROR: METHYLATED STATE NOT CLASSIFIED CORRECTLY. RUN ON MORE DATA.")

print("Saving state data")
print(df)
df.to_csv(output_file, sep="\t", index=False, header=True)

    


