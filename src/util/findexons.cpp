#include <set>
#include <MacTypes.h>
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "Sequence.h"

#ifdef OPENMP
#include <omp.h>
#endif
// to define structs
struct DpMatrixRow {
    DpMatrixRow() {}
    DpMatrixRow(size_t prevPotentialId, int pathScore)
            : prevPotentialId(prevPotentialId), pathScore(pathScore) {}
    size_t prevPotentialId; // prev potential ID
    long pathScore; //path score
};

struct ExonCandidates{
    ExonCandidates(){}
    ExonCandidates(long score, std::vector<Matcher::result_t> candidates) : score(score), candidates(candidates){}
    long score;
    std::vector<Matcher::result_t> candidates;
};
class ExonFinder{
    public:

        Matcher::result_t exonFlipper(Matcher::result_t inputAlignment){
            int temp; // for value switch
            temp = inputAlignment.qStartPos;
            inputAlignment.qStartPos = inputAlignment.qEndPos;
            inputAlignment.qEndPos = temp;

            temp = inputAlignment.dbStartPos;
            inputAlignment.dbStartPos = inputAlignment.dbEndPos;
            inputAlignment.dbEndPos = temp;

            temp = inputAlignment.queryOrfStartPos;
            inputAlignment.queryOrfStartPos = inputAlignment.queryOrfEndPos;
            inputAlignment.queryOrfEndPos = temp;

            //variables to correct cigar dataNM_001298432.1
            char cigarTempSymbol;
            std::string cigarTempNumber;
            std::string cigarFlipped;
            std::string backtrace = inputAlignment.backtrace;
            //to correct cigar data
            for(size_t c=0; c<backtrace.size();c++){
                bool foundM = backtrace[c] == 'M';
                bool foundD = backtrace[c] == 'D';
                bool foundI = backtrace[c] == 'I';
                if(foundM || foundD || foundI){
                    cigarTempSymbol = foundM?'M':cigarTempSymbol;
                    cigarTempSymbol = foundD?'D':cigarTempSymbol;
                    cigarTempSymbol = foundI?'I':cigarTempSymbol;
                    cigarFlipped = cigarTempNumber+cigarTempSymbol + cigarFlipped;
                    cigarTempNumber = std::string();
                } else {
                    cigarTempNumber = cigarTempNumber + backtrace[c];
                }
            }// end of for loop statement to correct cigar data
            inputAlignment.backtrace = cigarFlipped;
            return inputAlignment;
        }

        //class costructer
        ExonFinder(IndexReader * tDbr, IndexReader * qDbr, unsigned int queryKey) : tDbr(tDbr), qDbr(qDbr), queryKey(queryKey){}
        //to do dynamic programming
        void findOptimalExons(std::vector<Matcher::result_t> & optimalExonSolution,
                              std::vector<Matcher::result_t> & exonPath, //orfResults
                              unsigned int thread_idx,
                              long & orfScore) {
            if(exonPath.size()==0)
                return;
            std::sort(exonPath.begin(), exonPath.end(), Matcher::compareByDbposQpos);
            std::sort(exonPath.begin(), exonPath.end(), Matcher::compareByDbkeyAndStrand);

//            std::sort(exonPath.begin(), exonPath.end(), Matcher::compareByStrand);

            std::vector<ExonCandidates> candidates = createPotentialExonCombinations(exonPath);
            for(size_t candidateIdx = 0; candidateIdx < candidates.size(); candidateIdx++){
                ExonCandidates & candidate = candidates[candidateIdx];
                //DpMatrixRow * dpMatrixRow = new DpMatrixRow[maxSeqLen + 1]; // object used in DP
                dpMatrixRow.clear();
                //to construct <exonPath> which carries exon cadidate data, save the data of intron candidate into stretcheVec
                std::sort(candidate.candidates.begin(), candidate.candidates.end(), Matcher::compareByDbposQpos);
                tempExonVec.clear();
                trimmedExonResult.clear();
                findExonBoundaries(trimmedExonResult, candidate.candidates, tempExonVec, thread_idx);
                //to set up <dpMatrixRow>
                for (size_t id = 0; id < trimmedExonResult.size(); id++) {
                    long cost2 = trimmedExonResult[id].qEndPos>trimmedExonResult[id].queryOrfEndPos ? COST_MAX : 0;
                    long score1 = length(trimmedExonResult[id])*trimmedExonResult[id].seqId;
                    long score2 = (trimmedExonResult[id].qEndPos == trimmedExonResult[id].queryOrfEndPos) ? 100 : 0;
                    long score3 = (trimmedExonResult[id].queryOrfStartPos==trimmedExonResult[id].qStartPos) ? 100 : 0;
                    long score =  score1  - cost2 + score2 + score3;
                    dpMatrixRow.emplace_back(DpMatrixRow(id,score));
                }
                long bestPathScore = INT_MIN;
                size_t lastPotentialExonInBestPath = 0;
                size_t currId;

                //std::cout << trimmedExonResult.size() <<"\n";

                for (size_t currExon = 0; currExon < trimmedExonResult.size(); currExon++) {

//                    std::cout << trimmedExonResult[currExon].qStartPos << "\t" << trimmedExonResult[currExon].qEndPos << "\t";
//                    std::cout << trimmedExonResult[currExon].dbStartPos+1 << "\t" << trimmedExonResult[currExon].dbEndPos+1 << "\t";
//                    std::cout << trimmedExonResult[currExon].queryOrfStartPos << "_" << trimmedExonResult[currExon].queryOrfEndPos << "\t";
//                    std::cout << trimmedExonResult[currExon].dbOrfStartPos+1 << "_" << trimmedExonResult[currExon].dbOrfEndPos+1 << "\t";
//                    std::cout << trimmedExonResult[currExon].dbKey<< "\t" << trimmedExonResult[currExon].backtrace<< "\n";

                    for (size_t prevExon = 0; prevExon < currExon; prevExon++) {
                        bool strand = trimmedExonResult[currExon].dbEndPos>trimmedExonResult[currExon].dbStartPos;
                        int intronLength = strand?trimmedExonResult[currExon].dbStartPos - trimmedExonResult[prevExon].dbEndPos+1:trimmedExonResult[prevExon].dbEndPos-trimmedExonResult[currExon].dbStartPos+1 ;
                        bool isNotTooLongIntron = (intronLength < INTRON_MAX);
                        bool isNotTooShortIntron = intronLength > INTRON_MIN;
                        bool sameOrf = prevExon==0 || (trimmedExonResult[currExon].qStartPos - trimmedExonResult[prevExon].qEndPos)%3==1;
                        if (isNotTooLongIntron && isNotTooShortIntron && sameOrf){
                            long exonGap = std::abs(trimmedExonResult[currExon].qStartPos - trimmedExonResult[prevExon].qEndPos -1);
                            //the longer gap the bigger cost
                            long cost1 = std::pow(exonGap,2)*100;
                            cost1 = (cost1>COST_MAX)? COST_MAX : cost1;
                            long cost2 = trimmedExonResult[currExon].qEndPos>trimmedExonResult[currExon].queryOrfEndPos ? COST_MAX : 0;
                            long score1 = length(trimmedExonResult[currExon])*trimmedExonResult[currExon].seqId;
                            long score2 = (trimmedExonResult[currExon].qEndPos == trimmedExonResult[currExon].queryOrfEndPos) ? 100 : 0;
                            long bestScorePrev = dpMatrixRow[prevExon].pathScore;
                            long currScoreWithPrev = bestScorePrev + score1 - cost1- cost2 + score2;

//                            std::cout << bestScorePrev << "\t" << exonGap << "\t" <<cost1 << "\t"<< score1 << "\t"<<currScoreWithPrev << "\n";

                            // update row of currPotentialExon in case of improvement:
                            if (currScoreWithPrev > dpMatrixRow[currExon].pathScore  ) {
                                dpMatrixRow[currExon].prevPotentialId = prevExon;
                                dpMatrixRow[currExon].pathScore = currScoreWithPrev;
                            } //end of if statement to update
                        } //end of if conditional statement for avoid overlap and ...
                    } //end of DP 2nd for loop statement

//                    std::cout << dpMatrixRow[currExon].pathScore << "\n";

                    //to update the global max in case of improvement:

                    if (dpMatrixRow[currExon].pathScore > bestPathScore) {
                        lastPotentialExonInBestPath = currExon;
                        bestPathScore = dpMatrixRow[currExon].pathScore;
                    } //end of if conditional statement
                } //end of DP 1st for loop statement

//                std::cout << "----------" << "\n";
//                std::cout << bestPathScore << "\n\n";

                //end of Dynamic Progamming
                //to update <optimalExonSolution>
                currId = lastPotentialExonInBestPath;
                bool isBestScore = false;
//                std::cout << bestPathScore << "\n";
                for(size_t i=0; i<candidates.size(); i++){
                    if(bestPathScore > candidates[i].score) {
                        isBestScore = true;
                    } else {
                        isBestScore = false;
                        break;
                    }
                }

                if(isBestScore&&dpMatrixRow.size()>0){
                    orfScore = bestPathScore;
                    optimalExonSolution.clear();
                    while (dpMatrixRow[currId].prevPotentialId != currId) {
                        optimalExonSolution.emplace_back(trimmedExonResult[currId]);
                        currId = dpMatrixRow[currId].prevPotentialId;
                    }
                    optimalExonSolution.emplace_back(trimmedExonResult[currId]);
                    std::sort(optimalExonSolution.begin(),optimalExonSolution.end(),Matcher::compareByDbposQpos);
                    candidate.score = bestPathScore;
                    dpMatrixRow.clear();
                } //end of if conditional statement
            }//end of for loop statement
        }// end of function

    private:
    // class variable
        const int INTRON_MAX = 500000;
        const int INTRON_MIN = 20;
        const int COST_MAX = 5000;
        IndexReader * tDbr;
        IndexReader * qDbr;
        unsigned int queryKey;
        std::vector<Matcher::result_t> exonPath;
        std::vector<Matcher::result_t> trimmedExonResult;
        std::vector<Matcher::result_t> tempExonVec;
        std::vector<DpMatrixRow> dpMatrixRow;


        std::vector<ExonCandidates> createPotentialExonCombinations(std::vector<Matcher::result_t> exonPath){
            std::vector<ExonCandidates> exonCombination;
            std::vector<Matcher::result_t> tempVector;
            //
            int prevDBKey = exonPath[0].dbKey;
            bool prevStrand = exonPath[0].dbEndPos > exonPath[0].dbStartPos;
            for(size_t i = 0; i < exonPath.size(); i++){
                int currDBKey = exonPath[i].dbKey;
                bool currStrand = exonPath[i].dbEndPos > exonPath[i].dbStartPos;
                if(prevStrand==currStrand && prevDBKey == currDBKey){
                    tempVector.emplace_back(exonPath[i]);
                }else{
                    exonCombination.emplace_back(ExonCandidates(0, tempVector));
                    tempVector.clear();
                    prevDBKey = exonPath[i].dbKey;
                    prevStrand = exonPath[i].dbEndPos > exonPath[i].dbStartPos;
                    tempVector.emplace_back(exonPath[i]);
                }
            }
            if(tempVector.size()>0){
                exonCombination.emplace_back(0, tempVector);
            }
            tempVector.clear();
            return exonCombination;
        }



    // to find AG GT
        bool isAG(char * targetSeq, int index){
            char nt1 = std::toupper(targetSeq[index]);
            char nt2 = std::toupper(targetSeq[index+1]);
            if(nt1=='A'&&nt2=='G')
                return true;
            else
                return false;
        }
        bool isGT(char * targetSeq, int index){
            char nt1 = std::toupper(targetSeq[index]);
            char nt2 = std::toupper(targetSeq[index+1]);
            if((nt1=='G'&&nt2=='T') || (nt1=='G'&&nt2=='C') )
                return true;
            else
                return false;
        }
        bool isCT(char * targetSeq, int index){
            char nt1 = std::toupper(targetSeq[index]);
            char nt2 = std::toupper(targetSeq[index+1]);
            if(nt1=='C' && nt2=='T')
                return true;
            else
                return false;
        }
        bool isAC(char * targetSeq, int index){
            char nt1 = std::toupper(targetSeq[index]);
            char nt2 = std::toupper(targetSeq[index+1]);

            if((nt1=='A' && nt2=='C') || (nt1=='G' && nt2=='C'))
                return true;
            else
                return false;
        }

    // to build target and query sequence
        char * targetSequence(unsigned int dbKey, unsigned int thread_idx){
            size_t targetId = tDbr->sequenceReader->getId(dbKey);
            char * targetSeq = tDbr->sequenceReader->getData(targetId, thread_idx);
            return targetSeq;
        }
        char * querySequence(int thread_idx){
            unsigned int queryId = qDbr->sequenceReader->getId(queryKey);
            char * qSeq = qDbr->sequenceReader->getData(queryId, thread_idx);
            return qSeq;
        }

        // to update cigar
        int firstMatch(std::string cigar){
            int mPos = cigar.find('M');
            std::string tempNumber="";
            for(int c=0; c<mPos; c++){
                tempNumber = tempNumber + cigar[c];
            }
            return std::stoi(tempNumber);
        }
        int lastMatch(std::string cigar){
            int iPos =  cigar.find('I');
            int dPos =  cigar.find('D');
            int pos = iPos>dPos?iPos+1:dPos+1;
            std::string tempNumber="";
            for(size_t c=pos; c<cigar.size()-1; c++){
                tempNumber = tempNumber + cigar[c];
            }
            return std::stoi(tempNumber);
        }
        std::string firstCigar(std::string cigar){
            std::string newCigar = "";
            std::string cigarTempNumber = "";
            char  cigarTempSymbol;
            for(size_t c=0; c<cigar.size();c++){
                bool foundM = cigar[c] == 'M';
                bool foundD = cigar[c] == 'D';
                bool foundI = cigar[c] == 'I';
                if(foundM || foundD || foundI){
                    cigarTempSymbol = foundM?'M':cigarTempSymbol;
                    cigarTempSymbol = foundD?'D':cigarTempSymbol;
                    cigarTempSymbol = foundI?'I':cigarTempSymbol;
                    break;
                } else {
                    cigarTempNumber = cigarTempNumber + cigar[c];
                }
            }// end of for loop statement to correct cigar data
            newCigar = cigarTempNumber + cigarTempSymbol;
            return newCigar;
        }
        std::string lastCigar(std::string cigar){
            std::string newCigar = "";
            std::string cigarTempNumber = "";
            char  cigarTempSymbol;
            int cigarSize = cigar.size()-1;
            for(int c=cigarSize; c>=0; c--){
                bool foundM = cigar[c] == 'M';
                bool foundD = cigar[c] == 'D';
                bool foundI = cigar[c] == 'I';
                if(foundM || foundD || foundI) {
                    if(cigarTempSymbol=='M'||cigarTempSymbol=='D'||cigarTempSymbol=='I')
                        break;
                    cigarTempSymbol = foundM ? 'M' : cigarTempSymbol;
                    cigarTempSymbol = foundD ? 'D' : cigarTempSymbol;
                    cigarTempSymbol = foundI ? 'I' : cigarTempSymbol;
                } else {
                    cigarTempNumber = cigar[c] + cigarTempNumber;
                }
            }// end of for loop statement to correct cigar data
            newCigar = cigarTempNumber+cigarTempSymbol;
            return newCigar;
        }

        std::string  cigarUpdateAG(std::string cigar, int overlap){
            std::string newCigar = "";
            std::string cigarTempNumber = "";
            char  cigarTempSymbol;
            int tempNumber;
            for(size_t c=0; c<cigar.size();c++){
                bool foundM = cigar[c] == 'M';
                bool foundD = cigar[c] == 'D';
                bool foundI = cigar[c] == 'I';
                if(foundM || foundD || foundI){
                    cigarTempSymbol = foundM?'M':cigarTempSymbol;
                    cigarTempSymbol = foundD?'D':cigarTempSymbol;
                    cigarTempSymbol = foundI?'I':cigarTempSymbol;
                    tempNumber = std::stoi(cigarTempNumber);
                    if(tempNumber > overlap){
                        newCigar = newCigar + std::to_string(tempNumber-overlap)+cigarTempSymbol;
                        overlap = 0;
                    } else{
                        overlap -= tempNumber;
                    }
                    cigarTempNumber = "";
                } else {
                    cigarTempNumber = cigarTempNumber + cigar[c];
                }
            }// end of for loop statement to correct cigar data
            return newCigar;
        }
        std::string  cigarUpdateGT(std::string cigar, int overlap){
            std::string newCigar = "";
            std::string cigarTempNumber = "";
            char  cigarTempSymbol;
            int tempNumber;
            int cigarSize = cigar.size()-1;
            for(int c=cigarSize; c>=0; c--){
                bool foundM = cigar[c] == 'M';
                bool foundD = cigar[c] == 'D';
                bool foundI = cigar[c] == 'I';
                if(foundM || foundD || foundI) {
                    tempNumber = cigarTempNumber!=""?std::stoi(cigarTempNumber):0;
                    if(tempNumber>overlap && tempNumber>0){
                        newCigar = std::to_string(tempNumber - overlap)+cigarTempSymbol+ newCigar;
                        overlap = 0;
                    } else{
                        overlap -= tempNumber;
                    }
                    cigarTempSymbol = foundM ? 'M' : cigarTempSymbol;
                    cigarTempSymbol = foundD ? 'D' : cigarTempSymbol;
                    cigarTempSymbol = foundI ? 'I' : cigarTempSymbol;
                    cigarTempNumber = "";
                } else {
                    cigarTempNumber = cigar[c] + cigarTempNumber;
                }
            }// end of for loop statement to correct cigar data
            tempNumber = cigarTempNumber!=""?std::stoi(cigarTempNumber):0;
            newCigar = std::to_string(tempNumber - overlap)+cigarTempSymbol+ newCigar;
            return newCigar;
        }

        //to update identity
        int cigarLength(std::string cigar){
            std::string tempNumb = "";
            int length = 0;
            for (size_t i=0; i<cigar.size(); i++){
                if (cigar[i]=='M' ||cigar[i]=='I' ||cigar[i]=='D'){
                    length += std::stoi(tempNumb);
                    tempNumb = "";
                } else {
                    tempNumb = tempNumb + cigar[i];
                }
            }
            return length;

        }
        float identityUpdateAG(std::string tempCigar, float tempId, int overlapLength, int qPos, int targetPos, char * targetSeq, char  * qSeq){
            int originalCigarLength = cigarLength(tempCigar);
            int originalMatches = round(originalCigarLength*tempId);
            int newCigarLength = originalCigarLength - overlapLength;
            int newMatches = 0;
            int matchCnt = 0;
            float newid = tempId;

            if(overlapLength>0){
                for(int i=1; i <overlapLength+1; i++){
                    if(qSeq[qPos+i]==targetSeq[targetPos+i])
                        matchCnt++;
                }
                newMatches = originalMatches - matchCnt;
                newid = (float)newMatches/newCigarLength;
            } else{
                for(int i=overlapLength; i<0; i++){
                    if(qSeq[qPos+i]==targetSeq[targetPos+i])
                        matchCnt++;
                }
                newMatches = originalMatches + matchCnt;
                newid = (float)newMatches/newCigarLength;
            }
            newid = (newid>1)?1:newid;
            newid = (newid<0||isnan(newid))?0:newid;
            return newid;
        }
            float identityUpdateGT(std::string tempCigar, float tempId, int overlapLength, int qPos, int targetPos, char * targetSeq, char  * qSeq){
            int originalCigarLength = cigarLength(tempCigar);
            int originalMatches = round(originalCigarLength*tempId);
            int newCigarLength = originalCigarLength - overlapLength;
            int newMatches = 0;
            int matchCnt = 0;
            float newid = tempId;
            if(overlapLength>0){
                for(int i=1; i <overlapLength+1; i++){
                    if(qSeq[qPos-i]==targetSeq[targetPos-i])
                        matchCnt++;
                }
                newMatches = originalMatches - matchCnt;
                newid = (float)newMatches/newCigarLength;
            } else{
                for(int i=overlapLength; i<0; i++){
                    if(qSeq[qPos-i]==targetSeq[targetPos-i])
                        matchCnt++;
                }
                newMatches = originalMatches + matchCnt;
                newid = (float)newMatches/newCigarLength;
            }
            newid = (newid>1)?1:newid;
            newid = (newid<0||isnan(newid))?0:newid;
            return newid;
        }

        float identityUpdateCT(std::string tempCigar, float tempId, int overlapLength, int qPos, int targetPos, char * targetSeq, char  * qSeq){
            int originalCigarLength = cigarLength(tempCigar);
            int originalMatches = round(originalCigarLength*tempId);
            int newCigarLength = originalCigarLength - overlapLength;
            int newMatches = 0;
            int matchCnt = 0;
            float newid = tempId;
            if(overlapLength>0){
                for(int i=1; i <overlapLength+1; i++){
                    if(qSeq[qPos+i]==targetSeq[targetPos-i])
                        matchCnt++;
                }
                newMatches = originalMatches - matchCnt;
                newid = (float)newMatches/newCigarLength;
            } else{
                for(int i=overlapLength; i<0; i++){
                    if(qSeq[qPos+i]==targetSeq[targetPos-i])
                        matchCnt++;
                }
                newMatches = originalMatches + matchCnt;
                newid = (float)newMatches/newCigarLength;
            }
            newid = (newid>1)?1:newid;
            newid = (newid<0||isnan(newid))?0:newid;
            return newid;
        }

        float identityUpdateAC(std::string tempCigar, float tempId, int overlapLength, int qPos, int targetPos, char * targetSeq, char  * qSeq){
            int originalCigarLength = cigarLength(tempCigar);
            int originalMatches = round(originalCigarLength*tempId);
            int newCigarLength = originalCigarLength - overlapLength;
            int newMatches = 0;
            int matchCnt = 0;
            float newid = tempId;
            if(overlapLength>0){
                for(int i=1; i <overlapLength+1; i++){
                    if(qSeq[qPos-i]==targetSeq[targetPos+i])
                        matchCnt++;
                }
                newMatches = originalMatches - matchCnt;
                newid = (float)newMatches/newCigarLength;
            } else{
                for(int i=overlapLength; i<0; i++){
                    if(qSeq[qPos-i]==targetSeq[targetPos+i])
                        matchCnt++;
                }
                newMatches = originalMatches + matchCnt;
                newid = (float)newMatches/newCigarLength;
            }
            newid = (newid>1)?1:newid;
            newid = (newid<0||isnan(newid))?0:newid;
            return newid;
        }
        int length(Matcher::result_t exon){
            return  exon.qEndPos - exon.qStartPos +1;
        }
        int dbLength(Matcher::result_t exon){
            return  abs(exon.dbEndPos - exon.dbStartPos) +1;
        }


        // to cut off AT and GT
        void findExonBoundaries(
                std::vector<Matcher::result_t> & trimmedExonResult,
                std::vector<Matcher::result_t> & exonPath,
                std::vector<Matcher::result_t> & tempExonVec,
                unsigned int thread_idx){

            std::sort(exonPath.begin(), exonPath.end(), Matcher::compareByDbkeyAndStrand);
            int standardInScope = 100; //100
            int standardOutScope = 2;
            int inScope;
            int outScope;
            char * targetSeq = targetSequence(exonPath[0].dbKey, thread_idx);
            char * qSeq = querySequence(thread_idx);

//            std::cout << "orf\n";
            for(size_t exon=0; exon<exonPath.size(); exon++) {
                bool isForward = exonPath[exon].dbStartPos < exonPath[exon].dbEndPos;
//                char * targetSeq = targetSequence(exonPath[exon].dbKey, thread_idx);
//                char * qSeq = querySequence(thread_idx);
                outScope = standardOutScope;
                inScope = length(exonPath[exon])>standardInScope?standardInScope:length(exonPath[exon]);

//                std::cout << exonPath[exon].qStartPos << "\t" << exonPath[exon].qEndPos << "\t";
//                std::cout << exonPath[exon].dbStartPos+1 << "\t" << exonPath[exon].dbEndPos+1 << "\t";
//                std::cout << exonPath[exon].seqId << "\t" << exonPath[exon].queryOrfStartPos << "\t" << exonPath[exon].queryOrfEndPos << "\t";
//                std::cout << exonPath[exon].dbKey << "\t" << exonPath[exon].dbOrfStartPos << "\t" << exonPath[exon].dbOrfEndPos << "\n";

                if(exonPath[exon].qStartPos == exonPath[exon].queryOrfStartPos){
                    tempExonVec.emplace_back(exonPath[exon]);
                    outScope = 0;
                }
                //find AG
                if (isForward) {
                    int dbPointer = exonPath[exon].dbStartPos-outScope;
                    int qPointer = exonPath[exon].qStartPos-outScope;
                    int overlapLength = 0-outScope;
                    int dbEnder = dbPointer+inScope;
                    std::string tempCigar = exonPath[exon].backtrace;
                    float tempId = exonPath[exon].seqId;
                    int tempQueryPos = exonPath[exon].qStartPos;
                    int tempTargetPos = exonPath[exon].dbStartPos;
                    while(dbPointer < dbEnder){
                        if(!isAG(targetSeq, dbPointer-2)){
                            dbPointer++;
                            qPointer++;
                            overlapLength++;
                            continue;
                        }
                        exonPath[exon].dbStartPos = dbPointer;
                        exonPath[exon].qStartPos = qPointer;
                        exonPath[exon].backtrace = cigarUpdateAG(tempCigar, overlapLength);
                        exonPath[exon].seqId = identityUpdateAG(tempCigar, tempId,overlapLength, tempQueryPos
                        , tempTargetPos, targetSeq, qSeq);
                        if(exonPath[exon].qEndPos > exonPath[exon].qStartPos){
                            std::string firstCigarInfo = firstCigar(exonPath[exon].backtrace);
                            char firstCigarChar = firstCigarInfo[firstCigarInfo.size()-1];
                            int firstCigarSize = std::stoi(firstCigarInfo.substr(0,firstCigarInfo.size()-1));
                            if(firstCigarChar=='I'){
                                exonPath[exon].qStartPos += firstCigarSize;
                                exonPath[exon].backtrace = exonPath[exon].backtrace.substr(firstCigarInfo.size(),exonPath[exon].backtrace.size());
                            } else if(firstCigarChar=='D'){
                                exonPath[exon].dbStartPos += firstCigarSize;
                                exonPath[exon].backtrace = exonPath[exon].backtrace.substr(firstCigarInfo.size(),exonPath[exon].backtrace.size());
                            }
                            exonPath[exon].seqId = identityUpdateAG(tempCigar, tempId,overlapLength, tempQueryPos
                                    , tempTargetPos, targetSeq, qSeq);
                            tempExonVec.emplace_back(exonPath[exon]);
                        }
                        dbPointer++;
                        qPointer++;
                        overlapLength++;
                    }
                } else {
                    int dbPointer = exonPath[exon].dbStartPos+outScope;
                    int qPointer = exonPath[exon].qStartPos-outScope;
                    int overlapLength = 0-outScope;
                    int dbEnder = dbPointer-inScope;
                    std::string tempCigar = exonPath[exon].backtrace;
                    float tempId = exonPath[exon].seqId;
                    int tempQueryPos = exonPath[exon].qStartPos;
                    int tempTargetPos = exonPath[exon].dbEndPos;
                    while(dbPointer > dbEnder){
                        if(!isCT(targetSeq, dbPointer+1)){
                            dbPointer--;
                            qPointer++;
                            overlapLength++;
                            continue;
                        }
                        exonPath[exon].dbStartPos = dbPointer;
                        exonPath[exon].qStartPos = qPointer;
                        exonPath[exon].backtrace= cigarUpdateAG(tempCigar, overlapLength);
                        if(exonPath[exon].qEndPos > exonPath[exon].qStartPos){
                            std::string firstCigarInfo = firstCigar(exonPath[exon].backtrace);
                            char firstCigarChar = firstCigarInfo[firstCigarInfo.size()-1];
                            int firstCigarSize = std::stoi(firstCigarInfo.substr(0,firstCigarInfo.size()-1));
                            if(firstCigarChar=='I'){
                                exonPath[exon].qStartPos += firstCigarSize;
                                exonPath[exon].backtrace = exonPath[exon].backtrace.substr(firstCigarInfo.size(),exonPath[exon].backtrace.size());
                            } else if(firstCigarChar=='D'){
                                exonPath[exon].dbStartPos -= firstCigarSize;
                                exonPath[exon].backtrace = exonPath[exon].backtrace.substr(firstCigarInfo.size(),exonPath[exon].backtrace.size());
                            }
                            exonPath[exon].seqId = identityUpdateAG(tempCigar, tempId,overlapLength, tempQueryPos
                                    , tempTargetPos, targetSeq, qSeq);
                            tempExonVec.emplace_back(exonPath[exon]);
                        }
                        dbPointer--;
                        qPointer++;
                        overlapLength++;
                    }

                } //else statement end
            } //for loop end
//            std::cout << "<<\n";
            for(unsigned int trimmedExon=0; trimmedExon<tempExonVec.size(); trimmedExon++){

//                std::cout<<tempExonVec[trimmedExon].dbStartPos+1<<"_"<<tempExonVec[trimmedExon].dbEndPos+1 << "\t" << tempExonVec[trimmedExon].qStartPos<<"_"<<tempExonVec[trimmedExon].qEndPos << "\n";
//                std::cout<<tempExonVec[trimmedExon].queryOrfStartPos<<"_"<<tempExonVec[trimmedExon].queryOrfEndPos << "\t" << tempExonVec[trimmedExon].qStartPos<<"_"<<tempExonVec[trimmedExon].qEndPos << "++\n";

                bool isForward = tempExonVec[trimmedExon].dbStartPos < tempExonVec[trimmedExon].dbEndPos;
//                char * targetSeq = targetSequence(tempExonVec[trimmedExon].dbKey, thread_idx);
//                char * qSeq = querySequence(thread_idx);
                inScope =  standardInScope; //lastMatch(tempExonVec[trimmedExon].backtrace)>standardInScope?standardInScope:lastMatch(tempExonVec[trimmedExon].backtrace);
                inScope = length(tempExonVec[trimmedExon])>standardInScope?standardInScope:length(tempExonVec[trimmedExon]);
                outScope = standardOutScope;
//                if(isForward && isStpF(targetSeq, tempExonVec[trimmedExon].dbEndPos+1)){
//                    trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
//                    outScope = 0;
//                }else if(!isForward&&isStpR(targetSeq, tempExonVec[trimmedExon].dbEndPos-1)){
//                    trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
//                    outScope = 0;
//                }
                if(tempExonVec[trimmedExon].qEndPos == tempExonVec[trimmedExon].queryOrfEndPos){
                    trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
                    outScope = 0;
                }

                if(isForward){
                    int dbPointer = tempExonVec[trimmedExon].dbEndPos-inScope;
                    int qPointer = tempExonVec[trimmedExon].qEndPos-inScope;
                    int overlapLength = inScope;
                    int dbEnder = tempExonVec[trimmedExon].dbEndPos;
                    std::string tempCigar = tempExonVec[trimmedExon].backtrace;
                    float tempId = tempExonVec[trimmedExon].seqId;
                    int tempQueryPos = tempExonVec[trimmedExon].qEndPos;
                    int tempTargetPos = tempExonVec[trimmedExon].dbEndPos;
                    while(dbPointer < dbEnder+3+outScope){
                        if (dbPointer > tempExonVec[trimmedExon].dbOrfEndPos )
                            break;
                        if(!isGT(targetSeq, dbPointer+1)){
                            dbPointer++;
                            qPointer++;
                            overlapLength--;
                            continue;
                        }
                        tempExonVec[trimmedExon].dbEndPos = dbPointer;
                        tempExonVec[trimmedExon].qEndPos = qPointer;
                        tempExonVec[trimmedExon].backtrace = cigarUpdateGT(tempCigar, overlapLength);
                        if(tempExonVec[trimmedExon].qEndPos > tempExonVec[trimmedExon].qStartPos){
                            std::string lastCigarInfo = lastCigar(tempExonVec[trimmedExon].backtrace);
                            char lastCigarChar = lastCigarInfo[lastCigarInfo.size()-1];
                            int lastCigarSize = std::stoi(lastCigarInfo.substr(0,lastCigarInfo.size()-1));
                            if(lastCigarChar=='I'){
                                tempExonVec[trimmedExon].qEndPos -= lastCigarSize;
                                tempExonVec[trimmedExon].backtrace = tempExonVec[trimmedExon].backtrace.substr(0,tempExonVec[trimmedExon].backtrace.size()-lastCigarInfo.size());
                            } else if(lastCigarChar=='D'){
                                tempExonVec[trimmedExon].dbEndPos -= lastCigarSize;
                                tempExonVec[trimmedExon].backtrace = tempExonVec[trimmedExon].backtrace.substr(0,
                                        tempExonVec[trimmedExon].backtrace.size()-lastCigarInfo.size());
                            }
                            tempExonVec[trimmedExon].seqId = identityUpdateGT(tempCigar, tempId,overlapLength, tempQueryPos
                                    , tempTargetPos, targetSeq, qSeq);
                            trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
                        }
                        dbPointer++;
                        qPointer++;
                        overlapLength--;
                    }
                } else {
                    int dbPointer = tempExonVec[trimmedExon].dbEndPos+inScope;
                    int qPointer = tempExonVec[trimmedExon].qEndPos-inScope;
                    int overlapLength = inScope+outScope;
                    int dbEnder = tempExonVec[trimmedExon].dbEndPos;
                    std::string tempCigar = tempExonVec[trimmedExon].backtrace;
                    float tempId = tempExonVec[trimmedExon].seqId;
                    int tempQueryPos = tempExonVec[trimmedExon].qEndPos;
                    int tempTargetPos = tempExonVec[trimmedExon].dbEndPos;

                    while(dbPointer > dbEnder-3-outScope){
                        if (dbPointer < tempExonVec[trimmedExon].dbOrfEndPos)
                            break;
                        if(!isAC(targetSeq, dbPointer-2)){
                            dbPointer--;
                            qPointer++;
                            overlapLength--;
                            continue;
                        }
                        tempExonVec[trimmedExon].dbEndPos = dbPointer;
                        tempExonVec[trimmedExon].qEndPos = qPointer;
                        tempExonVec[trimmedExon].backtrace = cigarUpdateGT(tempCigar, overlapLength);
                        if(tempExonVec[trimmedExon].qEndPos > tempExonVec[trimmedExon].qStartPos){
                            std::string lastCigarInfo = lastCigar(tempExonVec[trimmedExon].backtrace);
                            char lastCigarChar = lastCigarInfo[lastCigarInfo.size()-1];
                            int lastCigarSize = std::stoi(lastCigarInfo.substr(0,lastCigarInfo.size()-1));
                            if(lastCigarChar=='I'){
                                tempExonVec[trimmedExon].qEndPos -= lastCigarSize;
                                tempExonVec[trimmedExon].backtrace = tempExonVec[trimmedExon].backtrace.substr(0,tempExonVec[trimmedExon].backtrace.size()-lastCigarInfo.size());
                            } else if(lastCigarChar=='D'){
                                tempExonVec[trimmedExon].dbEndPos += lastCigarSize;
                                tempExonVec[trimmedExon].backtrace = tempExonVec[trimmedExon].backtrace.substr(0,
                                        tempExonVec[trimmedExon].backtrace.size()-lastCigarInfo.size());
                            }
                            tempExonVec[trimmedExon].seqId = identityUpdateGT(tempCigar, tempId,overlapLength, tempQueryPos
                                    , tempTargetPos, targetSeq, qSeq);
                            trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
                        }
                        dbPointer--;
                        qPointer++;
                        overlapLength--;
                    }
                }//else statement end
            } // for loop end
        }//end of function


};//end of class

int findexons(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, (DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA) );
    IndexReader tDbr(par.db2, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, (DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA) );
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    Debug::Progress progress(alnDbr.getSize());
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> inputAlignments;
        std::vector<Matcher::result_t> orfResults;
        std::vector<Matcher::result_t> optimalExonSolution;
        std::vector<ExonCandidates> optimalSolutionWithScore;
        char buffer[2048];

#pragma omp  for schedule(dynamic, 10)
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            progress.updateProgress();
            const unsigned int queryKey = alnDbr.getDbKey(i);
            ExonFinder exonFinder(&tDbr, &qDbr, queryKey);
            char *data = alnDbr.getData(i, thread_idx);
            inputAlignments.clear();
            Matcher::readAlignmentResults(inputAlignments, data, true);

            std::sort(inputAlignments.begin(), inputAlignments.end(), Matcher::compareByDbposQpos);
            std::sort(inputAlignments.begin(), inputAlignments.end(), Matcher::compareByStrand);
            std::sort(inputAlignments.begin(), inputAlignments.end(), Matcher::compareOrfStartOrfEnd);


            int prevQueryOrfStartPos = inputAlignments[0].queryOrfStartPos;
            int prevQueryOrfEndPos = inputAlignments[0].queryOrfEndPos;
            resultWriter.writeStart(thread_idx);
            long orfScore = 0;
            long maxScore = 0;

//            if(queryKey!= 16)
//                continue;

//            if (inputAlignments.size()==0)
//                std::cout << queryKey << "\n";
            optimalSolutionWithScore.clear();
            for(size_t resIdx = 0; resIdx <inputAlignments.size(); resIdx++){

                // to correct flipped exon candidate data
                if(inputAlignments[resIdx].qStartPos>inputAlignments[resIdx].qEndPos){
                    inputAlignments[resIdx] = exonFinder.exonFlipper(inputAlignments[resIdx]);
                } // end of if conditional statement to correct flipped exons


                //subset orf

                bool currStrand = inputAlignments[resIdx].dbEndPos > inputAlignments[resIdx].dbStartPos;
                bool prevStrand = inputAlignments[resIdx-1].dbEndPos > inputAlignments[resIdx-1].dbStartPos;
                bool isSameStrand = currStrand == prevStrand;

                bool isSubsetQueryOrf = (inputAlignments[resIdx].queryOrfEndPos <= inputAlignments[resIdx-1].queryOrfEndPos) && (inputAlignments[resIdx].queryOrfStartPos >= inputAlignments[resIdx-1].queryOrfStartPos);
                bool isSubsetDBOrf = currStrand? (inputAlignments[resIdx].dbOrfEndPos <= inputAlignments[resIdx-1].dbOrfEndPos) && (inputAlignments[resIdx].dbOrfStartPos >= inputAlignments[resIdx-1].dbOrfStartPos):(inputAlignments[resIdx].dbOrfEndPos >= inputAlignments[resIdx-1].dbOrfEndPos) && (inputAlignments[resIdx].dbOrfStartPos <= inputAlignments[resIdx-1].dbOrfStartPos);
                bool isSameDBKey = inputAlignments[resIdx].dbKey == inputAlignments[resIdx-1].dbKey;
                if (isSubsetQueryOrf && isSubsetDBOrf && isSameDBKey && isSameStrand){
                    inputAlignments[resIdx].queryOrfEndPos = inputAlignments[resIdx-1].queryOrfEndPos;
                    inputAlignments[resIdx].queryOrfStartPos = inputAlignments[resIdx-1].queryOrfStartPos;
                }

//                std::cout << inputAlignments[resIdx].qStartPos << "\t" << inputAlignments[resIdx].qEndPos << "\t";
//                std::cout << inputAlignments[resIdx].dbStartPos+1 << "\t" << inputAlignments[resIdx].dbEndPos+1 << "\t";
//                std::cout << inputAlignments[resIdx].queryOrfStartPos << "_" << inputAlignments[resIdx].queryOrfEndPos << "\t";
//                std::cout << inputAlignments[resIdx].dbOrfStartPos+1 << "_" << inputAlignments[resIdx].dbOrfEndPos+1 << "\t";
//                std::cout << inputAlignments[resIdx].dbKey<< "\t" << inputAlignments[resIdx].backtrace<< "+++\n";


                bool querySameOrf =  prevQueryOrfStartPos == inputAlignments[resIdx].queryOrfStartPos && prevQueryOrfEndPos == inputAlignments[resIdx].queryOrfEndPos;
                if(querySameOrf){
                    orfResults.emplace_back(inputAlignments[resIdx]);
//                    std::cout << inputAlignments[resIdx].dbStartPos << "\t" << inputAlignments[resIdx].dbEndPos << "\t" << inputAlignments[resIdx].queryOrfStartPos << "\t" <<temp << "A\n";
                }else{
//                    std::cout << "new orf\n";
                    exonFinder.findOptimalExons(optimalExonSolution, orfResults, thread_idx, orfScore);
                    orfResults.clear();
                    if(orfScore>maxScore){
                        optimalSolutionWithScore.emplace_back(ExonCandidates(orfScore, optimalExonSolution));
                        maxScore = orfScore;
                    }
                    orfResults.emplace_back(inputAlignments[resIdx]);
                    prevQueryOrfStartPos = inputAlignments[resIdx].queryOrfStartPos;
                    prevQueryOrfEndPos = inputAlignments[resIdx].queryOrfEndPos;
                }
            }
            //last orf info -> optimal
            if(orfResults.size() > 0){
                exonFinder.findOptimalExons(optimalExonSolution, orfResults, thread_idx, orfScore);
                orfResults.clear();
                if(orfScore>maxScore){
                    optimalSolutionWithScore.emplace_back(ExonCandidates(orfScore, optimalExonSolution));
                    maxScore = orfScore;
                }
            }
            if(optimalSolutionWithScore.size()>0) {
                size_t optimestPathSize = optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates.size();
                for(size_t optIdx = 0; optIdx < optimestPathSize; optIdx++){

                    size_t len = Matcher::resultToBuffer(buffer, optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx], true, false, true);
                    resultWriter.writeAdd(buffer, len, thread_idx);//result buffer, len, thread_idx

//                    std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].dbStartPos+1 <<"_" << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].dbEndPos+1 << "  \t";
//                    std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].dbEndPos  - optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].dbStartPos   << "  \t";
//                    std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].qStartPos <<"_" << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].qEndPos << "  \t";
//                    std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].queryOrfStartPos <<"_" << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].queryOrfEndPos << "  \t";
//                    std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].dbOrfStartPos <<"_" << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].dbOrfEndPos << "  \t";
//                    std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].backtrace<<"_"<<optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].seqId<< "\t";
//                    std::cout << queryKey << "_" << optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx].dbKey << "\n";

                }
            }

            resultWriter.writeEnd(queryKey, thread_idx);
        } //the first for-loop finished
    }//end of the main scope
    // tsv output
    resultWriter.close(true);
    alnDbr.close();
    return EXIT_SUCCESS;
} //end of method <findexon>


