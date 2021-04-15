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

            inputAlignment.backtrace = Util::reverseCigar(inputAlignment.backtrace);
            return inputAlignment;
        }

        //class costructer
        ExonFinder(IndexReader * tDbr, IndexReader * qDbr, unsigned int queryKey) : tDbr(tDbr), qDbr(qDbr), queryKey(queryKey){}
        //to do dynamic programming
        void findOptimalExons(std::vector<Matcher::result_t> & optimalExonSolution,
                              std::vector<Matcher::result_t> & exonPath, //orfResults
                              unsigned int thread_idx,
                              long & orfScore) {
            std::sort(exonPath.begin(), exonPath.end(), Matcher::compareByDbkeyAndStrand);

            std::vector<ExonCandidates> candidates = createPotentialExonCombinations(exonPath);
            for(size_t candidateIdx = 0; candidateIdx < candidates.size(); candidateIdx++){
                ExonCandidates & candidate = candidates[candidateIdx];
                //DpMatrixRow * dpMatrixRow = new DpMatrixRow[maxSeqLen + 1]; // object used in DP
                dpMatrixRow.clear();
                //to construct <exonPath> which carries exon cadidate data, save the data of intron candidate into stretcheVec
                tempExonVec.clear();
                trimmedExonResult.clear();
                findExonBoundaries(trimmedExonResult, candidate.candidates, tempExonVec, thread_idx);
                //to set up <dpMatrixRow>
                for (size_t id = 0; id < trimmedExonResult.size(); id++) {
                    bool isFirstExon = trimmedExonResult[id].queryOrfStartPos==trimmedExonResult[id].qStartPos;
                    bool metStp = trimmedExonResult[id].qEndPos>trimmedExonResult[id].queryOrfEndPos;
                    int cost1 = metStp? COST_MAX : 0;
//                    int score1 = isFirstExon ? queryOrfLength(trimmedExonResult[id])*0.1 : 0;
                    int score2 = queryLength(trimmedExonResult[id]) * trimmedExonResult[id].seqId;
//                    long score = score1 + score2 - cost1;
                    long score = score2 - cost1;
                    dpMatrixRow.emplace_back(DpMatrixRow(id,score));
                }
                long bestPathScore = INT_MIN;
                size_t lastPotentialExonInBestPath = 0;
                size_t currId;

                for (size_t currExon = 0; currExon < trimmedExonResult.size(); currExon++) {
//                    std::cout << trimmedExonResult[currExon].qStartPos << "\t" << trimmedExonResult[currExon].qEndPos << "\t";
//                    std::cout << trimmedExonResult[currExon].dbStartPos+1 << "\t" << trimmedExonResult[currExon].dbEndPos+1 << "\t";
//                    std::cout << trimmedExonResult[currExon].queryOrfStartPos << "_" << trimmedExonResult[currExon].queryOrfEndPos << "\t";
//                    std::cout << trimmedExonResult[currExon].dbOrfStartPos+1 << "_" << trimmedExonResult[currExon].dbOrfEndPos+1 << "\t";
//                    std::cout << trimmedExonResult[currExon].backtrace<< "\t";
//                    std::cout << trimmedExonResult[currExon].dbKey<< "\t" << trimmedExonResult[currExon].seqId<< "\n";
                    for (size_t prevExon = 0; prevExon < currExon; prevExon++) {
                        bool strand = trimmedExonResult[currExon].dbEndPos>trimmedExonResult[currExon].dbStartPos;
                        int intronLength = strand?trimmedExonResult[currExon].dbStartPos - trimmedExonResult[prevExon].dbEndPos+1:trimmedExonResult[prevExon].dbEndPos-trimmedExonResult[currExon].dbStartPos+1 ;
                        bool isNotTooLongIntron = (intronLength < INTRON_MAX);
                        bool isNotTooShortIntron = intronLength > INTRON_MIN;
                        bool sameOrf = prevExon==0 || (trimmedExonResult[currExon].qStartPos - trimmedExonResult[prevExon].qEndPos)%3==1;
                        bool notMetStp = trimmedExonResult[currExon].qEndPos <= trimmedExonResult[currExon].queryOrfEndPos;
                        bool isGoingForward = trimmedExonResult[currExon].qStartPos > trimmedExonResult[prevExon].qEndPos ;
//                        bool isGoingForward = trimmedExonResult[prevExon].qEndPos != trimmedExonResult[prevExon].queryOrfEndPos ;
                        if (isNotTooLongIntron && isNotTooShortIntron  && notMetStp && isGoingForward){ //&& sameOrf
//                            int exonGap = std::abs(trimmedExonResult[currExon].qStartPos - trimmedExonResult[prevExon].qEndPos -1)*10;
//                            int cost1 = std::min(exonGap, COST_MAX); // DO WE NEED THIS?
                            int cost2 = trimmedExonResult[currExon].qEndPos>trimmedExonResult[currExon].queryOrfEndPos ? COST_MAX : 0;
                            int score1 = queryLength(trimmedExonResult[currExon])*trimmedExonResult[currExon].seqId;
//                            int score2 = (trimmedExonResult[currExon].qEndPos == trimmedExonResult[currExon].queryOrfEndPos) ? queryOrfLength(trimmedExonResult[currExon])*0.1 : 0;
                            long bestScorePrev = dpMatrixRow[prevExon].pathScore;
//                            long currScoreWithPrev = bestScorePrev + score1  - cost2 - cost1 + score2;
                            long currScoreWithPrev = bestScorePrev + score1  - cost2;
                            // update row of currPotentialExon in case of improvement:
                            if (currScoreWithPrev > dpMatrixRow[currExon].pathScore  ) {
                                dpMatrixRow[currExon].prevPotentialId = prevExon;
                                dpMatrixRow[currExon].pathScore = currScoreWithPrev;
                            } //end of if statement to update
                        } //end of if conditional statement for avoid overlap and ...
                    } //end of DP 2nd for loop statement
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
                        optimalExonSolution.emplace_back(stpCodonExtension(trimmedExonResult[currId]));
                        currId = dpMatrixRow[currId].prevPotentialId;
                    }
                    optimalExonSolution.emplace_back(stpCodonExtension(trimmedExonResult[currId]));
                    std::sort(optimalExonSolution.begin(),optimalExonSolution.end(),Matcher::compareHitsByPosAndStrand);
                    candidate.score = bestPathScore;
                    dpMatrixRow.clear();
                } //end of if conditional statement
            }//end of for loop statement
        }// end of function

    private:
    // class variable
    const int INTRON_MAX = 500000;
    const int INTRON_MIN = 30;
    const int COST_MAX = 5000;
    IndexReader * tDbr;
    IndexReader * qDbr;
    unsigned int queryKey;
    std::vector<Matcher::result_t> exonPath;
    std::vector<Matcher::result_t> trimmedExonResult;
    std::vector<Matcher::result_t> tempExonVec;
    std::vector<DpMatrixRow> dpMatrixRow;

    typedef std::pair<char, int> cigarTuple;

    Matcher::result_t stpCodonExtension(Matcher::result_t inputExon){
//        std::cout<<inputExon.dbEndPos<<"\t";
        size_t lenSTPCodon = inputExon.dbStartPos < inputExon.dbEndPos ? 3 : -3;
        bool haveSTPCodon = inputExon.qEndPos == inputExon.queryOrfEndPos;
        inputExon.dbEndPos = haveSTPCodon ? (inputExon.dbEndPos + lenSTPCodon) : inputExon.dbEndPos;
//        std::cout<<inputExon.qStartPos << "\t" << inputExon.qEndPos << lenSTPCodon << "\t";
//        std::cout<<inputExon.dbEndPos<<"\n";
        return inputExon;
    }

    std::vector<ExonCandidates> createPotentialExonCombinations(std::vector<Matcher::result_t> exonPath){
        std::vector<ExonCandidates> exonCombination;
        std::vector<Matcher::result_t> tempVector;
        int prevDBKey = exonPath[0].dbKey;
        bool prevStrand = exonPath[0].dbEndPos > exonPath[0].dbStartPos;
        tempVector.emplace_back(exonPath[0]);
        for(size_t i = 1; i < exonPath.size(); i++){
            int currDBKey = exonPath[i].dbKey;
            bool currStrand = exonPath[i].dbEndPos > exonPath[i].dbStartPos;
            int distBetweenExons = currStrand ? exonPath[i].dbStartPos-exonPath[i-1].dbEndPos : exonPath[i].dbEndPos - exonPath[i-1].dbStartPos;
            bool tooLongIntron = distBetweenExons>INTRON_MAX;
            if(prevStrand==currStrand && prevDBKey == currDBKey && !tooLongIntron){
                tempVector.emplace_back(exonPath[i]);
            }else{
                exonCombination.emplace_back(0, tempVector);
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
    // ???
    std::vector<cigarTuple> cigarToTuple(std::string backtrace){
        std::vector<cigarTuple> returnTuple;
        size_t tmpPos = 0;
        int cnt = 0;
        for (size_t pos = 0; pos < backtrace.size(); pos++) {
            if(!isdigit(backtrace[pos])){
                int num = std::stoi(backtrace.substr(tmpPos, cnt));
                char cha = backtrace[pos];
                returnTuple.emplace_back(std::pair<char, int>(cha,num));
                tmpPos = pos+1;
                cnt=0;
            } else{
                cnt++;
            }
        }
        return returnTuple;
    }

    // to find Donnor Sites and Acceptor sites
    bool isAcceptorSiteF(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index]);
        char nt2 = std::toupper(targetSeq[index+1]);
        return nt1=='A'&&nt2=='G';
    }
    bool isDonorSitF(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index]);
        char nt2 = std::toupper(targetSeq[index+1]);
        return  (nt1=='G'&&nt2=='T') || (nt1=='G'&&nt2=='C');
    }
    bool isAcceptorSiteR(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index]);
        char nt2 = std::toupper(targetSeq[index+1]);
        return  nt1=='C' && nt2=='T';
    }
    bool isDonorSiteR(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index]);
        char nt2 = std::toupper(targetSeq[index+1]);
        return (nt1=='A' && nt2=='C') || (nt1=='G' && nt2=='C');
    }
    bool isStpCodonF(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index]);
        char nt2 = std::toupper(targetSeq[index+1]);
        char nt3 = std::toupper(targetSeq[index+2]);
        return  (nt1=='T'&&nt2=='G'&&nt3=='A') || (nt1=='T'&&nt2=='A'&&nt3=='A') || (nt1=='T'&&nt2=='A'&&nt3=='G');
    }
    bool isStpCodonR(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index]);
        char nt2 = std::toupper(targetSeq[index+1]);
        char nt3 = std::toupper(targetSeq[index+2]);
        return  (nt1=='T'&&nt2=='T'&&nt3=='A') || (nt1=='T'&&nt2=='C'&&nt3=='A')|| (nt1=='C'&&nt2=='T'&&nt3=='A');
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
    std::pair<std::string, int> cigarQueryPosUpdateAcceptorSite(std::string cigar, int overlap){
        int returnNumber = 0;
        std::string returnString;
        std::vector<cigarTuple> tupleVector = cigarToTuple(cigar);
        for(size_t cnt=0; cnt<tupleVector.size(); cnt++){
            int tempOverlap = std::min(overlap,tupleVector[cnt].second);
            switch (tupleVector[cnt].first) {
                case 'I':
                    returnNumber += tempOverlap;
                    tupleVector[cnt].second -= tempOverlap;
                    break;
                case  'D':
                    tupleVector[cnt].second -= tempOverlap;
                    overlap -= tempOverlap;
                    break;
                default:
                    tupleVector[cnt].second -= tempOverlap;
                    returnNumber += tempOverlap;
                    overlap -= tempOverlap;
                    break;
            }
            if(overlap==0)
                break;
        }
        for(size_t cnt=0; cnt<tupleVector.size(); cnt++){
            if(tupleVector[cnt].second==0)
                continue;

            if(returnString == "") {
                switch (tupleVector[cnt].first) {
                    case 'I':
                        returnNumber += tupleVector[cnt].second;
                        break;
                    case 'D':
                        overlap -= tupleVector[cnt].second;
                        break;
                    default:
                        returnString = std::to_string(tupleVector[cnt].second) + tupleVector[cnt].first + returnString;
                        break;
                }
            } else {
                returnString = std::to_string(tupleVector[cnt].second) + tupleVector[cnt].first + returnString;
            }
        }
        return std::pair<std::string,int>(returnString, returnNumber);
    }

    std::pair<std::string , int>  cigarQueryPosUpdateDonorSite(std::string cigar, int overlap){
        int returnNumber = 0;
        std::string returnString;
        std::vector<cigarTuple> tupleVector = cigarToTuple(cigar);
        for(size_t cnt=0; cnt<tupleVector.size(); cnt++){
            int tempOverlap = std::min(overlap,tupleVector[cnt].second);
            switch (tupleVector[cnt].first) {
                case 'I':
                    returnNumber += tempOverlap;
                    tupleVector[cnt].second -= tempOverlap;
                    break;
                case  'D':
                    tupleVector[cnt].second -= tempOverlap;
                    overlap -= tempOverlap;
                    break;
                default:
                    tupleVector[cnt].second -= tempOverlap;
                    returnNumber += tempOverlap;
                    overlap -= tempOverlap;
                    break;
            }
            if(overlap==0)
                break;
        }
        for(size_t cnt=0; cnt<tupleVector.size(); cnt++){
            if(tupleVector[cnt].second==0)
                continue;

            if(returnString == "") {
                switch (tupleVector[cnt].first) {
                    case 'I':
                        returnNumber += tupleVector[cnt].second;
                        break;
                    case 'D':
                        overlap -= tupleVector[cnt].second;
                        break;
                    default:
                        returnString = std::to_string(tupleVector[cnt].second) + tupleVector[cnt].first + returnString;
                        break;
                }
            } else {
                returnString = std::to_string(tupleVector[cnt].second) + tupleVector[cnt].first + returnString;
            }
        }

        return std::pair<std::string,int>(returnString, returnNumber);
    }

    //to update identity
    int cigarLength(std::string cigar){
        std::vector<cigarTuple> cigarTupleVec =  cigarToTuple(cigar);
        int returnNum = 0;
        for (size_t i=0; i < cigarTupleVec.size(); i++){
            returnNum += cigarTupleVec[i].second;
        }
        return  returnNum;
    }
    float matchRatio(std::string cigar){
        std::vector<cigarTuple> cigarTupleVec =  cigarToTuple(cigar);
        int returnNum1 = 0;
        int returnNum2 = 0;
        for (size_t i=0; i < cigarTupleVec.size(); i++){
            returnNum2 += cigarTupleVec[i].second;
            if(cigarTupleVec[i].first == 'M')
                returnNum1 += cigarTupleVec[i].second;
        }
        return  (float)returnNum1/returnNum2;
    }
    int queryOrfLength(Matcher::result_t exon){
        return  exon.queryOrfEndPos - exon.queryOrfStartPos +1;
    }
    int queryLength(Matcher::result_t exon){
        return  exon.qEndPos - exon.qStartPos +1;
    }
    int dbLength(Matcher::result_t exon){
        return  abs(exon.dbEndPos - exon.dbStartPos) +1;
    }
    bool metContainF(char * qSeq, int qStartPos, int qEndPos){
        int currPos = qStartPos;
        while(currPos<qEndPos-2){
            char nt1 = std::toupper(qSeq[currPos]);
            char nt2 = std::toupper(qSeq[currPos+1]);
            char nt3 = std::toupper(qSeq[currPos+2]);
            if(nt1=='A' && nt2=='T' && nt3=='G'){
                return true;
            }
            currPos += 3;
        }
        return false;
    }
    bool metContainR(char * qSeq, int qStartPos, int qEndPos){
        int currPos = qStartPos;
        while(currPos<qEndPos-2){
            char nt1 = std::toupper(qSeq[currPos]);
            char nt2 = std::toupper(qSeq[currPos+1]);
            char nt3 = std::toupper(qSeq[currPos+2]);
            if(nt1=='T' && nt2=='A' && nt3=='C'){
                return true;
            }
            currPos += 3;
        }
        return false;
    }


    // to cut off AG and GT
    void findExonBoundaries(
            std::vector<Matcher::result_t> & trimmedExonResult,
            std::vector<Matcher::result_t> & exonPath,
            std::vector<Matcher::result_t> & tempExonVec,
            unsigned int thread_idx){

        int standardInScope = 150;
        int standardOutScope = 2;
        int inScope;
        int outScope;
        char * targetSeq = targetSequence(exonPath[0].dbKey, thread_idx);

        for(size_t exon=0; exon<exonPath.size(); exon++) {
            bool isForward = exonPath[exon].dbStartPos < exonPath[exon].dbEndPos;
            outScope = standardOutScope;
            inScope = std::min((int)(dbLength(exonPath[exon])*0.7), standardInScope);
            float matchIdentity = exonPath[exon].seqId / matchRatio(exonPath[exon].backtrace);
            std::cout<<exonPath[exon].queryOrfStartPos<<std::endl;
            if(exonPath[exon].qStartPos == exonPath[exon].queryOrfStartPos){
                tempExonVec.emplace_back(exonPath[exon]);
                outScope = 0;
//            } else if(exonPath[exon].queryOrfStartPos - exonPath[exon].qStartPos < 30){
//                int residueLength = tempExonVec[trimmedExon].queryOrfEndPos - tempExonVec[trimmedExon].qEndPos;
//                int dbEndPos = tempExonVec[trimmedExon].dbEndPos;
//                if(isForward){
//                    for (int dbPos=dbEndPos; dbPos<dbEndPos+50; dbPos++){
//                        if (isStpCodonF(targetSeq, dbPos)){
//                            tempExonVec[trimmedExon].dbEndPos = dbPos;
//                            trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
//                            outScope = 0;
//                        }
//                    }
//                }
//                else{
//                    for (int dbPos=dbEndPos; dbPos>dbEndPos-50; dbPos--){
//                        if (isStpCodonR(targetSeq, dbPos)){
//                            tempExonVec[trimmedExon].dbEndPos = dbPos;
//                            trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
//                            outScope = 0;
//                        }
//                    }
//                }
//            }

            //find AG
            if (isForward) {
                int currDbPos = exonPath[exon].dbStartPos - outScope;
                int overlapLength = -outScope;
                const int dbScopeEndPos = exonPath[exon].dbStartPos + inScope;
                while(currDbPos < dbScopeEndPos){
                    if(!isAcceptorSiteF(targetSeq, std::max(0, currDbPos - 2))){
                        currDbPos++;
                        overlapLength++;
                        continue;
                    }

                    exonPath[exon].dbStartPos = currDbPos;
                    std::pair<std::string, int> cigarQueryPos = cigarQueryPosUpdateAcceptorSite(exonPath[exon].backtrace,overlapLength);

                    exonPath[exon].backtrace = cigarQueryPos.first;
                    exonPath[exon].qStartPos +=cigarQueryPos.second;
                    exonPath[exon].seqId = matchRatio(exonPath[exon].backtrace) * matchIdentity;
                    //TODO: reduce computation by reuse already trimmed cigar string
                    overlapLength = 0;
                    tempExonVec.emplace_back(exonPath[exon]);
                    currDbPos++;
                }
            } else {
                int currDbPos = exonPath[exon].dbStartPos+outScope;
                int overlapLength = -outScope;
                int dbScopeEndPos = currDbPos - inScope;
                while(currDbPos > dbScopeEndPos){
                    if( !isAcceptorSiteR(targetSeq,currDbPos+1) ){
                        currDbPos--;
                        overlapLength++;
                        continue;
                    }
                    exonPath[exon].dbStartPos = currDbPos;
                    std::pair<std::string, int> cigarQueryPos = cigarQueryPosUpdateAcceptorSite(exonPath[exon].backtrace, overlapLength);
                    exonPath[exon].qStartPos += cigarQueryPos.second;
                    exonPath[exon].backtrace = cigarQueryPos.first;
                    exonPath[exon].seqId = matchRatio(exonPath[exon].backtrace) * matchIdentity;
                    overlapLength = 0;
                    tempExonVec.emplace_back(exonPath[exon]);
                    currDbPos--;
                }

            } //else statement end
        } //for loop end
        for(unsigned int trimmedExon=0; trimmedExon<tempExonVec.size(); trimmedExon++){
            bool isForward = tempExonVec[trimmedExon].dbStartPos < tempExonVec[trimmedExon].dbEndPos;
            inScope = std::min( (int)(dbLength(tempExonVec[trimmedExon])*0.7), standardInScope);
            outScope = standardOutScope;
            if(tempExonVec[trimmedExon].qEndPos == tempExonVec[trimmedExon].queryOrfEndPos){
                trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
                outScope = 0;
            }
            else if(tempExonVec[trimmedExon].queryOrfEndPos - tempExonVec[trimmedExon].qEndPos < 30){
                int residueLength = tempExonVec[trimmedExon].queryOrfEndPos - tempExonVec[trimmedExon].qEndPos;
                int dbEndPos = tempExonVec[trimmedExon].dbEndPos;
                if(isForward){
                    for (int dbPos=dbEndPos; dbPos<dbEndPos+50; dbPos++){
                        if (isStpCodonF(targetSeq, dbPos)){
                            tempExonVec[trimmedExon].dbEndPos = dbPos;
                            trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
                            outScope = 0;
                        }
                    }
                }
                else{
                    for (int dbPos=dbEndPos; dbPos>dbEndPos-50; dbPos--){
                        if (isStpCodonR(targetSeq, dbPos)){
                            tempExonVec[trimmedExon].dbEndPos = dbPos;
                            trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
                            outScope = 0;
                        }
                    }
                }
            }
            float matchIdentity = tempExonVec[trimmedExon].seqId/matchRatio(tempExonVec[trimmedExon].backtrace);
            if(isForward){
                int currDbPos = tempExonVec[trimmedExon].dbEndPos - inScope;
                int overlapLength = inScope;
                int dbScopeEndPos = tempExonVec[trimmedExon].dbEndPos+outScope;
                std::string tempCigar = tempExonVec[trimmedExon].backtrace;
                int tempQueryPos = tempExonVec[trimmedExon].qEndPos;
                while(currDbPos < dbScopeEndPos + 3 + outScope){
                    if (currDbPos > tempExonVec[trimmedExon].dbOrfEndPos + outScope )
                        break;
                    if(!isDonorSitF(targetSeq, std::max(currDbPos + 1, 0))){
                        currDbPos++;
                        overlapLength--;
                        continue;
                    }
                    tempExonVec[trimmedExon].dbEndPos = currDbPos;
                    std::pair<std::string , int> cigarQueryPos = cigarQueryPosUpdateDonorSite(tempCigar,
                                                                                              overlapLength);
                    tempExonVec[trimmedExon].qEndPos = tempQueryPos - cigarQueryPos.second;
                    tempExonVec[trimmedExon].backtrace = cigarQueryPos.first;
                    tempExonVec[trimmedExon].seqId = matchRatio(tempExonVec[trimmedExon].backtrace)*matchIdentity;
                    trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
                    currDbPos++;
                    overlapLength--;
                }
            } else {
                int currDbPos = tempExonVec[trimmedExon].dbEndPos + inScope;
                int overlapLength = inScope;
                int dbScopeEndPos = tempExonVec[trimmedExon].dbEndPos-outScope;
                std::string tempCigar = tempExonVec[trimmedExon].backtrace;
                int tempQueryPos = tempExonVec[trimmedExon].qEndPos;
                //int temp = tempExonVec[trimmedExon].dbStartPos;
                while(currDbPos > dbScopeEndPos - 3 - outScope){
                    if (currDbPos < tempExonVec[trimmedExon].dbOrfEndPos - outScope)
                        break;
                    if( !isDonorSiteR(targetSeq, currDbPos - 2) ){
                        currDbPos--;
                        overlapLength--;
                        continue;
                    }
                    tempExonVec[trimmedExon].dbEndPos = currDbPos;
                    std::pair<std::string , int> cigarQueryPos = cigarQueryPosUpdateDonorSite(tempCigar,
                                                                                              overlapLength);
                    tempExonVec[trimmedExon].qEndPos = tempQueryPos - cigarQueryPos.second;
                    tempExonVec[trimmedExon].backtrace = cigarQueryPos.first;
                    tempExonVec[trimmedExon].seqId = matchRatio(tempExonVec[trimmedExon].backtrace)*matchIdentity;
                    trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
                    currDbPos--;
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
            if(data[0]=='\0'){
                resultWriter.writeData("",0,queryKey,thread_idx);
                continue;
            }
            inputAlignments.clear();
            Matcher::readAlignmentResults(inputAlignments, data, true);
            std::sort(inputAlignments.begin(), inputAlignments.end(), Matcher::compareOrfStartOrfEnd);

            int prevQueryOrfStartPos = inputAlignments[0].queryOrfStartPos;
            int prevQueryOrfEndPos = inputAlignments[0].queryOrfEndPos;
            resultWriter.writeStart(thread_idx);
            long orfScore = 0;
            long maxScore = 0;


            optimalSolutionWithScore.clear();
            for(size_t resIdx = 0; resIdx <inputAlignments.size(); resIdx++){

                // to make sure that query position is never reversed, only db position can be reversed
                // In default we search only on the formward frame 1,2,3 so this function is not called
                // It is only important if search also on the backward frame!

//                std::cout << inputAlignments[resIdx].qStartPos << "\t" << inputAlignments[resIdx].qEndPos << "\t";
//                std::cout << inputAlignments[resIdx].dbStartPos+1 << "\t" << inputAlignments[resIdx].dbEndPos+1 << "\t";
//                std::cout << inputAlignments[resIdx].queryOrfStartPos << "_" << inputAlignments[resIdx].queryOrfEndPos << "\t";
//                std::cout << inputAlignments[resIdx].dbOrfStartPos+1 << "_" << inputAlignments[resIdx].dbOrfEndPos+1 << "\t";
//                std::cout << inputAlignments[resIdx].dbKey << "\t" << inputAlignments[resIdx].score << "\t" << inputAlignments[resIdx].eval <<"\t";
//                std::cout << inputAlignments[resIdx].seqId<< "\t" << inputAlignments[resIdx].backtrace<< "+++\n";

//                // SEQUENCE IDENTITY THNRESHOLD
//                if(inputAlignments[resIdx].seqId<0.65)
//                    continue;
//                if(inputAlignments[resIdx].dbKey!=0)
//                    continue;


                if(inputAlignments[resIdx].qStartPos>inputAlignments[resIdx].qEndPos){
                    inputAlignments[resIdx] = exonFinder.exonFlipper(inputAlignments[resIdx]);
                } // end of if conditional statement to correct flipped exon

                bool querySameOrf =  prevQueryOrfStartPos == inputAlignments[resIdx].queryOrfStartPos && prevQueryOrfEndPos == inputAlignments[resIdx].queryOrfEndPos;
                if(querySameOrf){
                    orfResults.emplace_back(inputAlignments[resIdx]);
                }else{
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
                // last vector in optimalSolutionWithScore contains highest scoring solution
                size_t optimestPathSize = optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates.size();
                for(size_t optIdx = 0; optIdx < optimestPathSize; optIdx++){

                    size_t len = Matcher::resultToBuffer(buffer, optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx], true, false, true);
                    resultWriter.writeAdd(buffer, len, thread_idx);//result buffer, len, thread_idx
//                    if (queryKey == 102){
//                        std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].dbStartPos + 1 << "_" <<
//                                  optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                           1].candidates[optIdx].dbEndPos + 1 << "  \t";
//                        std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].dbEndPos -
//                                     optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].dbStartPos << "  \t";
//                        std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].qStartPos << "_"
//                                  << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].qEndPos << "  \t";
//                        std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].queryOrfStartPos << "_"
//                                  << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].queryOrfEndPos << "  \t";
//                        std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].dbOrfStartPos << "_"
//                                  << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].dbOrfEndPos << "  \t";
//                        std::cout << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].backtrace << "_"
//                                  << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                              1].candidates[optIdx].seqId << "\t";
//                        std::cout << queryKey << "_" << optimalSolutionWithScore[optimalSolutionWithScore.size() -
//                                                                                 1].candidates[optIdx].dbKey << "\n";
//                    }

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


