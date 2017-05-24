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

    # Setup variables
    fileGroup = []
    currentGroup = []
    currentLems = fileList[0][4]

    # Group files into list of lists. Each sublist has files from one LEMS
    for dataFile in fileList:
        if dataFile[4] == currentLems:
            currentGroup.append(dataFile)
        else:
            if dataFile[0] != '.':
                fileGroup.append(currentGroup)
                currentGroup = [dataFile]
                currentLems = dataFile[4]

    fileGroup.append(currentGroup)  # For last group

    # Iterate through groups and cat CSV files
    for group in fileGroup:
        # Prepare file to write
        if (group[0][-3:] != "CSV") or ('LEMS' not in group[0]):
            continue
        currentLems = group[0][4]
        outFileName = "LEMS"+currentLems+"_Latest.CSV"
        print("Creating", outFileName)
        fOut = open(outFileName, "w")

        # Write first file with header
        for line in open(group[0]):
            fOut.write(line)

        # Write rest of files, skip header
        for dataFile in group[1:]:
            if ('Latest' in dataFile):
                continue
            print("    Opening", dataFile, "and appending to", outFileName)
            inFile = open(dataFile)
            for index, line in enumerate(inFile):
                if (index == 0):
                    continue
                fOut.write(line)
            inFile.close()  # Close writee file
        fOut.close()    # Close outfile
        
