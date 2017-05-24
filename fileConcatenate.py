## fileConcatenate.py
# Nipun Gunawardena
# Concatenate series of LEMS files into one CSV file
# After running this, should run the LEMS_Complete_Averaging.m script in matlab to 
# average and convert .csv files to .mat files (which are smaller)

# In future, might run this script directly from within matlab


import os


if __name__ == "__main__":
    # Get Files
    fileList = os.listdir(path='.')
    fileList.sort()

    # Filter file list so only LEMS data files remain
    fileList = [i for i in fileList if i[0] == 'L' and 'CSV' in i and 'Latest' not in i]

    # Setup variables
    fileGroup = []                  # List of lists. Each list is a set of files from one LEMS
    currentGroup = []               # Used to make file group
    currentLems = fileList[0][4]    # Start LEMS

    # Group files into list of lists. Each sublist has files from one LEMS
    for dataFile in fileList:
        if dataFile[4] == currentLems:
            currentGroup.append(dataFile)
        else:
            fileGroup.append(currentGroup)
            currentGroup = [dataFile]
            currentLems = dataFile[4]

    fileGroup.append(currentGroup)  # For last group

    # Iterate through groups and cat CSV files
    for group in fileGroup:
        # Prepare file to write
        currentLems = group[0][4]
        outFileName = "LEMS"+currentLems+"_Latest.CSV"
        print("Creating", outFileName)
        fOut = open(outFileName, "w")

        # Write first file with header
        for line in open(group[0]):
            fOut.write(line)

        # Write rest of files, skip header
        for dataFile in group[1:]:
            print("    Opening", dataFile, "and appending to", outFileName)
            inFile = open(dataFile)
            for index, line in enumerate(inFile):
                if (index == 0):
                    continue
                fOut.write(line)
            inFile.close()  # Close original file
        fOut.close()    # Close outfile (concatenated file)
        
