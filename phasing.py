import numpy as np
import pandas as pd


def phase(input_url, ref_ind=6, show_all=True):
    """
    Implements the first stage of phasing.
    Egg in the first trio is chosen as reference by default.
    Output the results into BED files, each per cell.
    :param input_url: The path of input genotype file.
    :param ref_ind: The index of column used as reference, starting from 0.
    :param show_all: If True, show all regions; If False, only show phases without breakpoints.
    :return: 1 if succeed, 0 otherwise.
    """
    # read in the input file and store all column names
    try:
        df = pd.read_csv(input_url, sep='\t')
    except:
        print("Can't read in the input file. File may not exist or tab-delimited.")
        return 0
    names = df.columns
    data = np.array(df).astype(str)
    if (data.shape[1] - 4) % 3 != 0:
        print("Incorrect number of columns.")
        return 0

    # define which regions to show
    if show_all:
        to_show = -2
    else:
        to_show = 0

    # loop over all SNPs
    for j in range(4, data.shape[1]):
        # open the file for output
        url = "output_" + names[j] + ".bed"
        f = open(url, 'w')

        # initialize variables
        chrom = data[0, 1]
        start = end = data[0, 2]
        old = phase = None
        for i in range(data.shape[0]):
            # calculate the phase value
            # Only heterozygous maternal SNPs are informative for phasing
            if data[i, 3] != 'AB':
                phase = -2
            # SNPs that are not called in reference or target cell
            elif data[i, ref_ind] == 'NC' or data[i, j] == 'NC':
                phase = -1
            # SNPs that are in phase with the reference
            elif data[i, ref_ind] == data[i, j]:
                phase = 1
            # SNPs that are heterozygous
            elif data[i, j] == 'AB':
                phase = 0.5
            # SNPs that are not in phase with the reference
            else:
                phase = 0

            # define the first region
            if old is None:
                # homozygous maternal SNPs are not counted
                if phase != -2:
                    old = phase
                    chrom = data[i, 1]
                continue

            # merge regions with the same phase
            # extend the region if the chromosome doesn't change and the phase remains
            if chrom == data[i, 1] and phase == old:
                end = data[i, 2]
            # update the region when phase or chromosome changes
            else:
                # output phases (and breakpoints)
                if old >= to_show:
                    f.write("chr" + chrom + "\t" + start + "\t" + end + "\t" + str(old) + "\n")
                start = end = data[i, 2]
                old = phase
                chrom = data[i, 1]
            # output the last region
            if i == data.shape[0] - 1:
                if old >= to_show:
                    f.write("chr" + chrom + "\t" + start + "\t" + end + "\t" + str(old) + "\n")

        f.close()
    return 1


if __name__ == '__main__':
    phase("input_data_trios.txt", show_all=False)

