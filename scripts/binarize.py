# Imports
import pandas as pd
import numpy as np
import time
import os
import os.path
from os import path
from argparse import ArgumentParser

def binarize_merged(sample, data_files, bin_size, lower_mv_bound, upper_mv_bound, data_assignment_setting, start_time, output_path, temp_files_path, chrom_length_file, k_threshold):  
    
    data_array = []
    data_array.append(data_files)

    for file in data_array:        
        print(file)
        df = pd.read_csv(open(file), sep='\t', names=["chr", "start", "end", "methylated", "unmethylated"])        
        length_dictionary, chromosomes = generate_length_data(df, chrom_length_file)        
        num_chrom = len(chromosomes)        
        print(f"{num_chrom} chromosomes found.")        
        inters_chrom = list(set(chromosomes) & set(length_dictionary.keys()))
        num_chrom = len(inters_chrom)        
        print(f"{num_chrom} chromosomes retained.")        
        missed_chromosomes = []
        chrom_count = 0
        
        for chromosome in chromosomes:            
            chrom_length = length_dictionary.get(chromosome)
            if chrom_length is None:
                missed_chromosomes.append(chromosome)
            else:
                chrom_count += 1
                df_chrom_extracted = df[df["chr"] == chromosome]
                
                # get number of bins
                num_bins = int(np.ceil(chrom_length / bin_size))

                print(f"Binning data for {chromosome}")
                # Bin the data
                binned_data = binning(df_chrom_extracted, chromosome, bin_size, num_bins, data_assignment_setting)

                # Extract columns of interest
                binned_data = binned_data.drop(binned_data.columns[[0,1]], axis=1)
                binned_data_cpy = binned_data.copy(deep=True)
                print(f"Binarizing data for {chromosome}")
                df_val = option_2(binned_data_cpy, sample, chromosome, lower_mv_bound, upper_mv_bound, k_threshold)

                # Start file
                data_header = pd.DataFrame([], columns=[sample, chromosome])
                data_header.to_csv(os.path.join(temp_files_path,f'{sample}_{chromosome}_binary.txt'), sep="\t", index=False)
                col = ["methylated", "unmethylated"]                    
                # convert all data to dataframe
                df_val = pd.DataFrame(df_val, columns=col)   

                # append to file
                df_val.to_csv(os.path.join(temp_files_path,f'{sample}_{chromosome}_binary.txt'), mode='a', sep="\t", index=False)                
                data_table = pd.read_csv(os.path.join(temp_files_path,f'{sample}_{chromosome}_binary.txt'), low_memory=False, sep='\t', header=0)
                filename = f'{sample}_{chromosome}_binary.txt'
                data_table.to_csv(os.path.join(output_path,filename), sep="\t", index=False)                    
                print(f"{chromosome} done ({chrom_count}/{num_chrom}).",(time.time() - start_time), " seconds elapsed.")
        
        
def generate_length_data(df, chrom_length_file):
    # Get names of chromosomes that have data in the file
    chroms = df.iloc[:, 0].unique()
    length_dictionary = {}
    
    index = 0
    
    # Find the length of chromosomes
    with open(chrom_length_file) as file:
        for line in file:
            (key, value) = line.split()
            if str(key) in chroms:
                length_dictionary[str(key)] = int(value)
                index += 1
    return length_dictionary, chroms
    
def binning(df, chrom_id, bin_size, num_bins, data_assignment_setting):
    
    col = df.columns.values.tolist()
    df = df.copy()
    df = df.to_numpy() # chrm, start, end, methylated, unmethylated
    binned_data = np.zeros((num_bins, 4),  dtype=int) # start, end, methylated, unmethylated 
    data_index = 0
    
    for i in range(num_bins):
        # start
        binned_data[i][0] = i*bin_size
        #end
        binned_data[i][1] = (i+1)*bin_size - 1
    
    while(data_index < len(df)):
        start_nucleotide = df[data_index][1]
        end_nucleotide = df[data_index][2]
        
        if data_assignment_setting == 'left' or data_assignment_setting == 'right':
            if data_assignment_setting == 'left':
                bin_index = int(np.floor(start_nucleotide / bin_size))
            elif data_assignment_setting == 'right':
                bin_index = int(np.floor(end_nucleotide / bin_size))
            else:
                print("ERROR IN BINNING")
            
            try:
                bin_start = binned_data[bin_index][0]
                bin_end = binned_data[bin_index][1]
            except:
                print("bin index: ", bin_index, np.shape(binned_data))
                
            if (data_assignment_setting == 'left' and (start_nucleotide >= bin_start and start_nucleotide <= bin_end)) or (data_assignment_setting == 'right' and (end_nucleotide >= bin_start and end_nucleotide <= bin_end)):
                #add methylated
                binned_data[bin_index][2] += df[data_index][3]
                # add unmethylated
                binned_data[bin_index][3] += df[data_index][4]
            else:
                print("ERROR", bin_index, bin_start, start_nucleotide)
        
        elif data_assignment_setting == 'both':
            # Add to two bins
            bin_index = [int(np.floor(start_nucleotide / bin_size)), int(np.floor(end_nucleotide / bin_size))]
            bin_start = binned_data[bin_index[0]][0]
            bin_end = binned_data[bin_index[1]][1]
            
            # Add to first bin:
            
            #add methylated
            binned_data[bin_index[0]][2] += df[data_index][3]
            # add unmethylated
            binned_data[bin_index[0]][3] += df[data_index][4]
            
            # Add to first bin:
            
            #add methylated
            binned_data[bin_index[1]][2] += df[data_index][3]
            # add unmethylated
            binned_data[bin_index[1]][3] += df[data_index][4]
        
        else:
            print("ERROR IN BINNING", bin_index, bin_start, start_nucleotide)
        data_index += 1
    binned_data = pd.DataFrame(binned_data, columns = ["start", "end", "methylated", "unmethylated"])
    return binned_data

count_removed = 0

def option_2(df, sample_name, chrom_name, lower_prob, upper_prob, k_threshold):        
    df["methylation_value"] = df["methylated"] / (df["methylated"] + df["unmethylated"])
    df_val = df.to_numpy()
    # PASS TO O2
    np.apply_along_axis(o2, axis=1, arr=df_val, lower_prob=lower_prob, upper_prob=upper_prob, k_threshold=k_threshold)    
    # remove methylation value, the second column
    df_val = np.delete(df_val, 2, 1)  # delete second row of A
    df_val = df_val.astype(int)    
    print(f"{count_removed} bins removed due to low count.")    
    # return only the array and columns    
    return df_val


def o2(x, lower_prob, upper_prob, k_threshold):
    # x[0] is methylated, x[1] is unmethylated, x[2] is methylation value    
    if ((x[0] + x[1]) >= k_threshold):
        # unmethylated condition
        if x[2] < lower_prob:
            x[0] = 0
            x[1] = 1
        # methylated condition
        elif x[2] > upper_prob:
            x[0] = 1
            x[1] = 0
        # intermediatelly methylated
        elif x[2] >= lower_prob and x[2] <= upper_prob:
            x[0] = 1
            x[1] = 1
        else:
            x[0] = 0
            x[1] = 0
    else:
        if (x[0] + x[1]) > 0:
            global count_removed
            count_removed +=1
        x[0] = 0
        x[1] = 0    
    return x

parser = ArgumentParser()
parser.add_argument("-i", dest="meth_data")
parser.add_argument("-b", dest="bin_size")
parser.add_argument("-l", dest="im_low")
parser.add_argument("-u", dest="im_high")
parser.add_argument("-s", dest="sample_name")
parser.add_argument("-d", dest="das")
parser.add_argument("-t", dest="temp_folder")
parser.add_argument("-k", dest="k_threshold")
parser.add_argument("-c", dest="chrom_lengths")
parser.add_argument("-o", dest="output_folder")

args = vars(parser.parse_args())

bin_size = int(args["bin_size"])
data_files = args["meth_data"]
sample = args["sample_name"]
lower_im_bound = float(args["im_low"])
upper_im_bound = float(args["im_high"])
data_assignment_setting = args["das"]
output_path = args["output_folder"]
temp_files_path = args["temp_folder"]
chrom_length_file = args["chrom_lengths"]
k_threshold = int(args["k_threshold"])

start_time = time.time()
   
print("Beginning binarization routine.")
binarize_merged(sample, data_files, bin_size, lower_im_bound, upper_im_bound, data_assignment_setting, start_time, output_path, temp_files_path, chrom_length_file, k_threshold)
print("Done. ",(time.time() - start_time), "seconds elapsed.")


