#!usr/bin/python

"""creates random start points across the human genome and will make us a fastafile"""
import random
def main():
    """googogogo"""
    #container to hold gene size
    chromosome = []
    size = []
    data = []
    runningTotal = []
    randomStartPoints = []
    humanGenomeSize= 0
    j = 0
    needThousand = 0
    outputFile = open("randomValues.fa","w")
    in_data = open('hg19.fa.fai', 'r')
    for line in in_data:
        line_no_end = line.rstrip('\n')
        # print "please gimme the good stuff: %s" % line_no_end
        list_container= []
        list_container = line_no_end.split("\t")
        #print "work maybe?: %s" % list_container[1]
        data.append(line_no_end)
        chromosome.append(list_container[0])
        size.append(list_container[1])

    for chrmSize in size:
        humanGenomeSize = humanGenomeSize + int(chrmSize)
        runningTotal.append(humanGenomeSize)
        
    in_data.close()
    #print runningTotal[0]
    #print chromosome[0]
    #print runningTotal[24]
    #print "final size = %i" %  humanGenomeSize
    for i in range(2000): #added extra just to be able to throw out when hitting Ns in a file
        rand = random.randint(0,humanGenomeSize)
        #print "ints going into variable: %i" % rand
        randomStartPoints.append(rand)
    for value in randomStartPoints:
        j  = 0 #finds the chromosome name for the file as well as where we want to stop
        while 1:
            if value < runningTotal[j]:
                break;
            else:
                j = j+1
        
        #print value
        #print runningTotal[j]
        if j == 0:
            valueToFind = value
        else:
            valueToFind = value-runningTotal[j-1]
        chrCount = 0
        basePairs = ""
        firstLine = "> '"+str(chromosome[j]) + " " + str(valueToFind) + " " + str(valueToFind +500) + " " + str(valueToFind + 1000)+" +'\n"
        print firstLine
        end = 0 
        filename= '/home/groupdirs/epigenetics/sequences/hg19/'+  chromosome[j]+ ".fa"
        with open(filename) as file:
            while True:
                c = file.read(1)
                if not c:
                    print "end of file"
                    break
                if c == "\n":
                    continue
                chrCount = chrCount+1
                if chrCount == valueToFind:
                    break
            l = 0
            while l < 1000:
                c=file.read(1)
                if c =="\n":
                    continue
                if c == "N" or c == "n":
                    end = 1
                    break
                basePairs = basePairs + c
                l = l+1
            if end == 1:
                file.close()
                continue
            print basePairs
            outputFile.write(firstLine)
            outputFile.write(basePairs)
            outputFile.write("\n\n")
            needThousand = needThousand + 1
            if needThousand == 1000 :
                outputFile.close()
                file.close()
                quit()


        file.close()
        
                        #print chromosome[j]+"test" 

main ()
