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
    timeingRecordProcessors = np.zeros(shape = (processorNumber, 5), dtype = np.double) # sumSqrtNuMulVecsTime, sumSternheimerTime, sumIniGuessTime, sumCOCGTime
    
    for i in range(processorNumber):
        patternL1 = re.compile(r"nuChi0Eigscomm (\d+), rank (\d+), omega (\d+), sumOperatorMulTime_nuChi0EigsComm (\d+\.\d+) s, composed of sumSqrtNuMulVecsTime (\d+\.\d+) ms, sumSternheimerTime (\d+\.\d+) ms, sumCollectdRhoTime (\d+\.\d+) ms")
        patternL2 = re.compile(r"In sumSternheimerTime, sumIniGuessTime (\d+\.\d+) ms; sumCOCGTime (\d+\.\d+) ms, composed of sumLhsfunTime (\d+\.\d+) ms, sumSolveMuTime (\d+\.\d+) ms, sumMultipleTime (\d+\.\d+) ms")
        fileName = "timingRecords_rank" + str(i) + "_" + str(int(i/NPnuChi0comm)) + "nuChi0Eigscomm.log"
        with open(fileName, 'r') as theTimingRecord:
            sumOperatorMulTime = 0.0
            sumSqrtNuMulVecsTime = 0.0
            sumSternheimerTime = 0.0
            sumIniGuessTime = 0.0
            sumCOCGTime = 0.0
            for line in theTimingRecord:
                match1 = patternL1.search(line)
                if match1:
                    sumOperatorMulTime += float(match1.group(4))
                    sumSqrtNuMulVecsTime += float(match1.group(5))
                    sumSternheimerTime += float(match1.group(6))
                match2 = patternL2.search(line)
                if match2:
                    sumIniGuessTime += float(match2.group(1))
                    sumCOCGTime += float(match2.group(2))
        sumOperatorMulTime /= float(testTime)
        sumSqrtNuMulVecsTime /= float(testTime)
        sumSternheimerTime /= float(testTime)
        sumIniGuessTime /= float(testTime)
        sumCOCGTime /= float(testTime)
        timeingRecordProcessors[i, 0] = sumOperatorMulTime
        timeingRecordProcessors[i, 1] = sumSqrtNuMulVecsTime
        timeingRecordProcessors[i, 2] = sumSternheimerTime
        timeingRecordProcessors[i, 3] = sumIniGuessTime
        timeingRecordProcessors[i, 4] = sumCOCGTime

    mean = np.mean(timeingRecordProcessors, axis = 0)
    std = np.std(timeingRecordProcessors, axis = 0)
    longest = np.max(timeingRecordProcessors, axis = 0)
    print(f"sumSqrtNuMulVecsTime mean {mean[1]:.1f} dev {std[1]:.1f} ms, longest {longest[1]:.1f} ms\n"
        f"sumSternheimerTime mean {mean[2]:.1f} dev {std[2]:.1f} ms, longest {longest[2]:.1f} ms\n"
        f"sumIniGuessTime mean {mean[3]:.1f} dev {std[3]:.1f} ms, longest {longest[3]:.1f} ms\n"
        f"sumCOCGTime mean {mean[4]:.1f} dev {std[4]:.1f} ms, longest {longest[4]:.1f} ms\n")
    
    P_maxOperatorMulTime = np.argmax(timeingRecordProcessors[:, 0])
    print(f"longest OperatorMulTime is {P_maxOperatorMulTime:d}, value {timeingRecordProcessors[P_maxOperatorMulTime, 0]:.3f}s (not ms)\n"
        f"its SqrtNuMulVecsTime {timeingRecordProcessors[P_maxOperatorMulTime, 1]/1000.0:.3f}s, its IniGuessTime {timeingRecordProcessors[P_maxOperatorMulTime, 3]/1000.0:.3f}s, its COCGTime {timeingRecordProcessors[P_maxOperatorMulTime, 4]/1000.0:.3f}s\n")

if __name__ == "__main__":
    main(sys.argv)