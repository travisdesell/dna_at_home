#!/usr/bin/python
"""combine the no repeat headers with the extended header files"""

def main():
    header = []
    lineDa = []
    count = 0
    with  open('wgEncodeOpenChromChipMcf7Pol2SerumstimRawDataRep1_deduplicated_slug_1000.fa', 'r') as headerData:
        for line in headerData:
            header.append(line)
    with  open('wgEncodeOpenChromChipMcf7Pol2SerumstimRawDataRep1_deduplicated_slug_1000_no_repeats.fa', 'r') as lineData:
        for line in lineData:
            lineDa.append(line)
            count = count+1

    with open ("output", "w") as output:
        for i in range (0, count):
            if count %3 == 0:
                output.write(header[i])
            else:
                output.write(lineDa[i])
                

main()
    



