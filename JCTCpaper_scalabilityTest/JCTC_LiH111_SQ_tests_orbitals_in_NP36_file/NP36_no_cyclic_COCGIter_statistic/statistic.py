import re
import sys
import numpy as np

def main(argv):
    if (len(argv) != 4) :
        print("input parameters does not fit requirement! processorNumber - NPnuChi0comm - testTime")
        quit()
    processorNumber = int(sys.argv[1])
    NPnuChi0comm = int(sys.argv[2])
    testTime = int(sys.argv[3])
    listCOCGIterRecordProcessors = []
    subspaceIterTime = [4, 3, 3, 2, 2, 2, 2, 2]
    N_r_perprocessor = int(576/(processorNumber/NPnuChi0comm))
    N_occ_perprocessor = int(8/NPnuChi0comm)
    for omegaIndex in range(8):
        nSternEquationsTheOmegaPerProcessor = subspaceIterTime[omegaIndex] * N_r_perprocessor * N_occ_perprocessor
        COCGIterRecordProcessors = np.zeros(shape = (processorNumber, nSternEquationsTheOmegaPerProcessor), dtype = np.int8) # sumSqrtNuMulVecsTime, sumSternheimerTime, sumIniGuessTime, sumCOCGTime
        listCOCGIterRecordProcessors.append(COCGIterRecordProcessors)
    
    patternL1 = re.compile(r"nuChi0Eigscomm (\d+), rank (\d+), omega (\d+), sumOperatorMulTime_nuChi0EigsComm (\d+\.\d+) s, composed of sumSqrtNuMulVecsTime (\d+\.\d+) ms, sumSternheimerTime (\d+\.\d+) ms, sumCollectdRhoTime (\d+\.\d+) ms")
    patternL2 = re.compile(r"start index (\d+), block size 1, block COCG iterated for (\d+) times; square of RHS_frobnorm (\d+\.\d+)E-(\d+), square of finalResidFrobnorm (\d+\.\d+)E-(\d+)")
    for i in range(processorNumber):
        fileName = "timingRecords_rank" + str(i) + "_" + str(int(i/NPnuChi0comm)) + "nuChi0Eigscomm.log"
        with open(fileName, 'r') as theTimingRecord:
            omegaIndex = 0
            equationIndex = 0
            COCGIterRecordProcessors_theOmega = listCOCGIterRecordProcessors[omegaIndex]
            for line in theTimingRecord:
                match2 = patternL2.search(line)
                if match2:
                    COCGIterRecordProcessors_theOmega[i, equationIndex] = int(match2.group(2))
                    equationIndex += 1
                    continue
                match1 = patternL1.search(line)
                if match1:
                    omegaIndex += 1
                    equationIndex = 0
                    if omegaIndex < 8:
                        COCGIterRecordProcessors_theOmega = listCOCGIterRecordProcessors[omegaIndex]
                
    for omegaIndex in range(8):
        COCGIterRecordProcessors_theOmega = listCOCGIterRecordProcessors[omegaIndex]    
        mean = np.mean(COCGIterRecordProcessors_theOmega)
        std = np.std(COCGIterRecordProcessors_theOmega)
        maxIter = np.max(COCGIterRecordProcessors_theOmega)
        minIter = np.min(COCGIterRecordProcessors_theOmega)
        print(f"omega {omegaIndex}: COCG iteration time mean {mean} dev {std}, max {maxIter}, min {minIter}\n")
    

if __name__ == "__main__":
    main(sys.argv)