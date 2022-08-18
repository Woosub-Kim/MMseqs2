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
        void findOptimalExons(
                    std::vector<Matcher::result_t> & optimalExonSolution,
                    std::vector<Matcher::result_t> & exonPath, //orfResults
                    unsigned int thread_idx,
                    long & orfScore,
                    float orfKeepingBonusRatio,
                    unsigned int trimmingSpliceSiteInScope,
                    unsigned int trimmingSpliceSiteOutScope,
                    unsigned int trimmingTerminusOutScope,
                    unsigned int trimmingTerminusInScope
                ) {
            std::sort(exonPath.begin(), exonPath.end(), Matcher::compareByDbkeyAndStrand);

            std::vector<ExonCandidates> candidates = createPotentialExonCombinations(exonPath);
            for(size_t candidateIdx = 0; candidateIdx < candidates.size(); candidateIdx++){
                ExonCandidates & candidate = candidates[candidateIdx];
                //DpMatrixRow * dpMatrixRow = new DpMatrixRow[maxSeqLen + 1]; // object used in DP
                dpMatrixRow.clear();
                //to construct <exonPath> which carries exon cadidate data, save the data of intron candidate into stretcheVec
                tempExonVec.clear();
                trimmedExonResult.clear();
                findExonBoundaries(trimmedExonResult, candidate.candidates, tempExonVec, thread_idx, trimmingSpliceSiteInScope, trimmingSpliceSiteOutScope, trimmingTerminusOutScope, trimmingTerminusInScope);
                //to set up <dpMatrixRow>
                for (size_t id = 0; id < trimmedExonResult.size(); id++) {
//                    bool isFirstExon = trimmedExonResult[id].queryOrfStartPos==trimmedExonResult[id].qStartPos;
                    bool metStp = trimmedExonResult[id].qEndPos>trimmedExonResult[id].queryOrfEndPos;
                    int cost = metStp? COST_MAX : 0;
                    int score1 = queryLength(trimmedExonResult[id]) * trimmedExonResult[id].seqId;
                    long score = score1 - cost;
                    dpMatrixRow.emplace_back(DpMatrixRow(id,score));
                }
                long bestPathScore = INT_MIN;
                size_t lastPotentialExonInBestPath = 0;
                size_t currId;
                for (size_t currExon = 0; currExon < trimmedExonResult.size(); currExon++) {
                    for (size_t prevExon = 0; prevExon < currExon; prevExon++) {
                        bool strand = trimmedExonResult[currExon].dbEndPos>trimmedExonResult[currExon].dbStartPos;
                        int intronLength = strand?trimmedExonResult[currExon].dbStartPos - trimmedExonResult[prevExon].dbEndPos+1:trimmedExonResult[prevExon].dbEndPos-trimmedExonResult[currExon].dbStartPos+1 ;
                        bool isNotTooLongIntron = (intronLength < INTRON_MAX);
                        bool isNotTooShortIntron = intronLength > INTRON_MIN;
                        bool isLastExon = trimmedExonResult[currExon].queryOrfEndPos==trimmedExonResult[currExon].qEndPos;
                        bool sameOrf = prevExon==0 || (trimmedExonResult[currExon].qStartPos - trimmedExonResult[prevExon].qEndPos)%3==1;
                        bool notMetStp = trimmedExonResult[currExon].qEndPos <= trimmedExonResult[currExon].queryOrfEndPos;
                        bool isGoingForward = trimmedExonResult[currExon].qStartPos > trimmedExonResult[prevExon].qEndPos ;
                        if (isNotTooLongIntron && isNotTooShortIntron  && notMetStp && isGoingForward){ //&& sameOrf
                            int cost = trimmedExonResult[currExon].qEndPos>trimmedExonResult[currExon].queryOrfEndPos ? COST_MAX : 0;
                            int score1 = queryLength(trimmedExonResult[currExon])*trimmedExonResult[currExon].seqId;
                            int score2 = sameOrf ? queryLength(trimmedExonResult[currExon])*orfKeepingBonusRatio:0;
                            long bestScorePrev = dpMatrixRow[prevExon].pathScore;
                            long currScoreWithPrev = bestScorePrev - cost + score1 + score2;
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
                //end of Dynamic Progamming
                //to update <optimalExonSolution>
                currId = lastPotentialExonInBestPath;
                bool isBestScore = false;
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
        size_t lenSTPCodon = inputExon.dbStartPos < inputExon.dbEndPos ? 3 : -3;
        bool haveSTPCodon = inputExon.qEndPos == inputExon.queryOrfEndPos;
        inputExon.dbEndPos = haveSTPCodon ? (inputExon.dbEndPos + lenSTPCodon) : inputExon.dbEndPos;
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
        //temp
        return  (nt1=='G'&&nt2=='T') ;//|| (nt1=='G'&&nt2=='C');
    }
    bool isAcceptorSiteR(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index]);
        char nt2 = std::toupper(targetSeq[index+1]);
        return  nt1=='C' && nt2=='T';
    }
    bool isDonorSiteR(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index]);
        char nt2 = std::toupper(targetSeq[index+1]);
        // temp
        return (nt1=='A' && nt2=='C');// || (nt1=='G' && nt2=='C');
    }
    bool isStpCodonF(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index+1]);
        char nt2 = std::toupper(targetSeq[index+2]);
        char nt3 = std::toupper(targetSeq[index+3]);
        return  (nt1=='T'&&nt2=='G'&&nt3=='A') || (nt1=='T'&&nt2=='A'&&nt3=='A') || (nt1=='T'&&nt2=='A'&&nt3=='G');
    }
    bool isStpCodonR(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index-3]);
        char nt2 = std::toupper(targetSeq[index-2]);
        char nt3 = std::toupper(targetSeq[index-1]);
        return  (nt1=='T'&&nt2=='T'&&nt3=='A') || (nt1=='T'&&nt2=='C'&&nt3=='A')|| (nt1=='C'&&nt2=='T'&&nt3=='A');
    }
    bool isMetCodonF(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index]);
        char nt2 = std::toupper(targetSeq[index+1]);
        char nt3 = std::toupper(targetSeq[index+2]);
        return (nt1=='A'&&nt2=='T'&&nt3=='G');
    }
    bool isMetCodonR(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index-2]);
        char nt2 = std::toupper(targetSeq[index-1]);
        char nt3 = std::toupper(targetSeq[index]);
        return (nt1=='C'&&nt2=='A'&&nt3=='T');
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

    std::string addCigar(std::string cigar, char symbol, int length){
        std::string newCigar = "";
        std::vector<cigarTuple> tupleVector = cigarToTuple(cigar);
        if (cigar != "" && tupleVector[tupleVector.size()-1].first == symbol) {
            tupleVector[tupleVector.size() - 1].second += length;
        } else if (length > 0) {
            tupleVector.emplace_back(std::pair<char, int>(symbol, length));
        }
        for(size_t cnt=0; cnt<tupleVector.size(); cnt++){
            newCigar = newCigar + std::to_string(tupleVector[cnt].second) + tupleVector[cnt].first;
        }
        return newCigar;
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
    bool firstExon(int qStartPos, int qOrfStartPos, int inScope, int outScope){
        return (qStartPos - outScope) < qOrfStartPos && (qOrfStartPos < qStartPos + inScope);
    }
    bool lastExon(int qEndPos, int qOrfEndPos, int inScope, int outScope){
        return (qEndPos-inScope < qOrfEndPos) && (qOrfEndPos < qEndPos + outScope);
    }
    // to cut off AG and GT
    void findExonBoundaries(
                std::vector<Matcher::result_t> & trimmedExonResult,
                std::vector<Matcher::result_t> & exonPath,
                std::vector<Matcher::result_t> & tempExonVec,
                unsigned int thread_idx,
                int trimmingSpliceSiteInScope,
                int trimmingSpliceSiteOutScope,
                int trimmingTerminusOutScope,
                int trimmingTerminusInScope
            ) {
        float maxTrimmingScopeRatio = 1;
        int inScope;
        int outScope;
        char * targetSeq = targetSequence(exonPath[0].dbKey, thread_idx);
        // find AGs
        for(size_t exon=0; exon<exonPath.size(); exon++) {
            bool isForward = exonPath[exon].dbStartPos < exonPath[exon].dbEndPos;
            outScope = trimmingSpliceSiteOutScope;
            inScope = std::min((int)(dbLength(exonPath[exon])*maxTrimmingScopeRatio), trimmingSpliceSiteInScope);
            float matchIdentity = exonPath[exon].seqId / matchRatio(exonPath[exon].backtrace);
            bool isFirst = firstExon(exonPath[exon].qStartPos, exonPath[exon].queryOrfStartPos, trimmingTerminusInScope, trimmingTerminusOutScope);
            bool tempFlag = false;
            if(exonPath[exon].qStartPos == exonPath[exon].queryOrfStartPos && ((isForward&&isMetCodonF(targetSeq,exonPath[exon].dbStartPos))||(!isForward&&isMetCodonR(targetSeq,exonPath[exon].dbStartPos)))){
                tempExonVec.emplace_back(exonPath[exon]);
                outScope = 0;
                //
                inScope = 0;
                //
                tempFlag = true;
            }
            if(isFirst && exonPath[exon].qStartPos != exonPath[exon].queryOrfStartPos){
                if (isForward){
                    int dbPos = exonPath[exon].dbStartPos + trimmingTerminusInScope;
                    // 0 base
//                    int dbPos = exonPath[exon].dbStartPos + trimmingTerminusInScope - exonPath[exon].qStartPos%3;
                    int originStart = exonPath[exon].dbStartPos;
                    while(dbPos >= exonPath[exon].dbStartPos - trimmingTerminusOutScope) {
                        if (isMetCodonF(targetSeq, dbPos)){
                            outScope = 0;
                            //
                            inScope = 0;
                            //
                            tempFlag = true;
                            if (originStart < dbPos){
                                std::pair<std::string, int> cigarQueryPos = cigarQueryPosUpdateAcceptorSite(exonPath[exon].backtrace, originStart - dbPos);
                                exonPath[exon].backtrace = cigarQueryPos.first;
                                exonPath[exon].seqId = matchRatio(exonPath[exon].backtrace) * matchIdentity;
                            }
                            exonPath[exon].dbStartPos = dbPos;
                            exonPath[exon].qStartPos = exonPath[exon].queryOrfStartPos;
                            tempExonVec.emplace_back(exonPath[exon]);
                        }
                        dbPos = dbPos - 3;
                    }
                } else {
                    int dbPos = exonPath[exon].dbStartPos - trimmingTerminusInScope;
                    // 0 base
//                    int dbPos = exonPath[exon].dbStartPos - trimmingTerminusInScope + exonPath[exon].qStartPos%3;
                    int originStart = exonPath[exon].dbStartPos;
                    while (dbPos <= exonPath[exon].dbStartPos  + trimmingTerminusOutScope) {
                        if (isMetCodonR(targetSeq, dbPos)){
                            outScope = 0;
                            //
                            inScope = 0;
                            //
                            tempFlag = true;
                            if (originStart > dbPos){
                                std::pair<std::string, int> cigarQueryPos = cigarQueryPosUpdateAcceptorSite(exonPath[exon].backtrace, dbPos - originStart);
                                exonPath[exon].backtrace = cigarQueryPos.first;
                                exonPath[exon].seqId = matchRatio(exonPath[exon].backtrace) * matchIdentity;
                            }
                            exonPath[exon].dbStartPos = dbPos;
                            exonPath[exon].qStartPos = exonPath[exon].queryOrfStartPos;
                            tempExonVec.emplace_back(exonPath[exon]);
                        }
                        dbPos = dbPos + 3;
                    }
                }
            }
            if (tempFlag)
                continue;
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
                    exonPath[exon].qStartPos += cigarQueryPos.second;
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
        // find GTs
        for(unsigned int trimmedExon=0; trimmedExon<tempExonVec.size(); trimmedExon++){
            bool isForward = tempExonVec[trimmedExon].dbStartPos < tempExonVec[trimmedExon].dbEndPos;
            inScope = std::min( (int)(dbLength(tempExonVec[trimmedExon])*maxTrimmingScopeRatio), trimmingSpliceSiteInScope);
            outScope = trimmingSpliceSiteOutScope;
            bool tempFlag = false;
            if(tempExonVec[trimmedExon].qEndPos == tempExonVec[trimmedExon].queryOrfEndPos&&((isForward&&isStpCodonF(targetSeq,tempExonVec[trimmedExon].dbEndPos))||(!isForward&&isStpCodonR(targetSeq,tempExonVec[trimmedExon].dbEndPos)))){
                trimmedExonResult.emplace_back(tempExonVec[trimmedExon]);
                outScope = 0;
                //
                inScope = 0;
                //
                tempFlag = true;
            }
            bool isLast = lastExon(tempExonVec[trimmedExon].qEndPos, tempExonVec[trimmedExon].queryOrfEndPos, trimmingTerminusInScope, trimmingTerminusOutScope);
            if(isLast && tempExonVec[trimmedExon].qEndPos != tempExonVec[trimmedExon].queryOrfEndPos){
                if (isForward){
                    int dbPos = tempExonVec[trimmedExon].dbEndPos;
                    // 0 base
//                    int dbPos = tempExonVec[trimmedExon].dbEndPos -1 -  tempExonVec[trimmedExon].qEndPos%3;
                    while(dbPos <= tempExonVec[trimmedExon].dbEndPos + trimmingTerminusOutScope) {
                        if (isStpCodonF(targetSeq, dbPos)){
                            outScope = 0;
                            //
                            inScope = 0;
                            //
                            tempFlag = true;
                            tempExonVec[trimmedExon].dbEndPos = dbPos;
                            tempExonVec[trimmedExon].qEndPos = tempExonVec[trimmedExon].queryOrfEndPos;
                            tempExonVec.emplace_back(tempExonVec[trimmedExon]);
                            break;
                        }
                        dbPos = dbPos + 3;
                    }
                } else {
                    int dbPos = tempExonVec[trimmedExon].dbEndPos;
                    // 0 base
//                    int dbPos = tempExonVec[trimmedExon].dbEndPos +1 + tempExonVec[trimmedExon].qEndPos%3;
                    while (dbPos >= tempExonVec[trimmedExon].dbEndPos  - trimmingTerminusOutScope) {
                        if (isStpCodonR(targetSeq, dbPos)){
                            outScope = 0;
                            //
                            inScope = 0;
                            //
                            tempFlag = true;
                            tempExonVec[trimmedExon].dbEndPos = dbPos;
                            tempExonVec[trimmedExon].qEndPos = tempExonVec[trimmedExon].queryOrfEndPos;
                            tempExonVec.emplace_back(tempExonVec[trimmedExon]);
                            break;
                        }
                        dbPos = dbPos - 3;
                    }
                }
            }
            if (tempFlag)
                continue;
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
                    std::pair<std::string , int> cigarQueryPos = cigarQueryPosUpdateDonorSite(tempCigar, overlapLength);
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
                while(currDbPos > dbScopeEndPos - 3 - outScope){
                    if (currDbPos < tempExonVec[trimmedExon].dbOrfEndPos - outScope)
                        break;
                    if( !isDonorSiteR(targetSeq, currDbPos - 2) ){
                        currDbPos--;
                        overlapLength--;
                        continue;
                    }
                    tempExonVec[trimmedExon].dbEndPos = currDbPos;
                    std::pair<std::string , int> cigarQueryPos = cigarQueryPosUpdateDonorSite(tempCigar, overlapLength);
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

    float orfKeepingBonusRatio = (float)par.orfBonusRatio/100;
    unsigned int trimmingSpliceSiteInScope = par.trimSpliceInScope;
    unsigned int trimmingSpliceSiteOutScope = par.trimSpliceOutScope;
    unsigned int trimmingTerminusInScope = par.trimTermInScope;
    unsigned int trimmingTerminusOutScope = par.trimTermOutScope;
    float falsePositiveFilteringRatio = (float)par.filteringRatio/100;

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
                if(inputAlignments[resIdx].qStartPos>inputAlignments[resIdx].qEndPos){
                    inputAlignments[resIdx] = exonFinder.exonFlipper(inputAlignments[resIdx]);
                } // end of if conditional statement to correct flipped exon

                bool querySameOrf =  prevQueryOrfStartPos == inputAlignments[resIdx].queryOrfStartPos && prevQueryOrfEndPos == inputAlignments[resIdx].queryOrfEndPos;
                if(querySameOrf){
                    orfResults.emplace_back(inputAlignments[resIdx]);
                }else{
                    exonFinder.findOptimalExons(optimalExonSolution, orfResults, thread_idx, orfScore, orfKeepingBonusRatio, trimmingSpliceSiteInScope, trimmingSpliceSiteOutScope, trimmingTerminusOutScope, trimmingTerminusInScope);
                    orfResults.clear();
                    if(orfScore>maxScore){
                        int length = 0;
                        for (size_t optExonIdx = 0; optExonIdx < optimalExonSolution.size(); optExonIdx++ ) {
                            length = length + abs(optimalExonSolution[optExonIdx].qEndPos - optimalExonSolution[optExonIdx].qStartPos) + 1;
                        }
                        float scoreLengthRatio = (float)orfScore/length;
                        if (scoreLengthRatio > falsePositiveFilteringRatio) {
                            optimalSolutionWithScore.emplace_back(ExonCandidates(orfScore, optimalExonSolution));
//                            maxScore = orfScore;
                        }
                        maxScore = orfScore;

                    }
                    orfResults.emplace_back(inputAlignments[resIdx]);
                    prevQueryOrfStartPos = inputAlignments[resIdx].queryOrfStartPos;
                    prevQueryOrfEndPos = inputAlignments[resIdx].queryOrfEndPos;
                }
            }
            //last orf info -> optimal
            if(orfResults.size() > 0){
                exonFinder.findOptimalExons(optimalExonSolution, orfResults, thread_idx, orfScore, orfKeepingBonusRatio, trimmingSpliceSiteInScope, trimmingSpliceSiteOutScope, trimmingTerminusOutScope, trimmingTerminusInScope);
                orfResults.clear();
                if(orfScore>maxScore){
                    int length = 0;
                    for (size_t optExonIdx = 0; optExonIdx < optimalExonSolution.size(); optExonIdx++ ) {
                        length = length + abs(optimalExonSolution[optExonIdx].qEndPos - optimalExonSolution[optExonIdx].qStartPos) + 1;
                    }
                    float scoreLengthRatio = (float)orfScore/length;
                    if (scoreLengthRatio > falsePositiveFilteringRatio) {
                        optimalSolutionWithScore.emplace_back(ExonCandidates(orfScore, optimalExonSolution));
//                        maxScore = orfScore;
                    }
                    maxScore = orfScore;

                }
            }

            // output
            if(optimalSolutionWithScore.size()>0) {
                // last vector in optimalSolutionWithScore contains highest scoring solution
                size_t optimestPathSize = optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates.size();
                for(size_t optIdx = 0; optIdx < optimestPathSize; optIdx++){
                    size_t len = Matcher::resultToBuffer(buffer, optimalSolutionWithScore[optimalSolutionWithScore.size()-1].candidates[optIdx], true, false, true);
                    resultWriter.writeAdd(buffer, len, thread_idx);//result buffer, len, thread_idx
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


