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

struct donorAcceptorSitesCandidate{
    donorAcceptorSitesCandidate(){}
    donorAcceptorSitesCandidate(int score, int dist, int dbDonorPos, int dbAcceptorPos, int qDonorPos, int qAcceptorPos)
        : score(score), dist(dist),  dbDdonorPos(dbDonorPos), dbAcceptorPos(dbAcceptorPos), qDonorPos(qDonorPos), qAcceptorPos(qAcceptorPos){}
    int score;
    int dist;
    int dbDdonorPos;
    int dbAcceptorPos;
    int qDonorPos;
    int qAcceptorPos;

};

struct splicingSiteCandidate{
    splicingSiteCandidate(){}
    splicingSiteCandidate(int score, int dbPos, int qPos) : score(score), dbPos(dbPos), qPos(qPos){}
    int score;
    int dbPos;
    int qPos;
};

bool compareSplicingSiteCandidate(splicingSiteCandidate &first, splicingSiteCandidate &second){
    return first.score <= second.score;
}
bool compareDonorAcceptorSitesCandidate(donorAcceptorSitesCandidate &first, donorAcceptorSitesCandidate &second){
    if (first.score != second.score){
        return first.score > second.score;
    }
    if (first.dist != second.dist){
        return first.dist < second.dist;
    }
    return true;
}

class ExonFinder{
public:

    //class costructer
//    ExonFinder(IndexReader * tDbr, IndexReader * qDbr, unsigned int queryKey) : tDbr(tDbr), qDbr(qDbr), queryKey(queryKey){}
    ExonFinder(IndexReader * tDbr) : tDbr(tDbr){}

    Matcher::result_t flipExons(Matcher::result_t inputAlignment){
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

    //Dynamic programming
    void findOptimalExons(
            std::vector<Matcher::result_t> & optimalExonSolution,
            std::vector<Matcher::result_t> & alignments,
            // ???
//            unsigned int thread_idx,
            long & orfScore
    ) {
        std::sort(alignments.begin(), alignments.end(), Matcher::compareByDbkeyAndStrand);
        std::vector<ExonCandidates> candidates = createPotentialExonCombinations(alignments);
        for(size_t candidateIdx = 0; candidateIdx < candidates.size(); candidateIdx++){
            ExonCandidates & candidate = candidates[candidateIdx];
            dpMatrixRow.clear();
            for (size_t id = 0; id < candidate.candidates.size(); id++) {
                long score = queryLength(candidate.candidates[id]) * candidate.candidates[id].seqId;
                dpMatrixRow.emplace_back(DpMatrixRow(id,score));
            }
            long bestPathScore = INT_MIN;
            size_t currId;
            for (size_t currExon = 0; currExon < candidate.candidates.size(); currExon++) {
                long score = queryLength(candidate.candidates[currExon])*candidate.candidates[currExon].seqId;
                bool strand = candidate.candidates[currExon].dbEndPos>candidate.candidates[currExon].dbStartPos;
                for (size_t prevExon = 0; prevExon < currExon; prevExon++) {
                    int tIntronLength = strand ? candidate.candidates[currExon].dbStartPos - candidate.candidates[prevExon].dbEndPos+1:candidate.candidates[prevExon].dbEndPos-candidate.candidates[currExon].dbStartPos+1 ;
                    int qIntronLength = candidate.candidates[currExon].qStartPos - candidate.candidates[prevExon].qEndPos+1;
                    bool isNotTooLongIntron = (tIntronLength < MAX_INTRON_LENGTH);
                    bool isNotTooShortIntron = tIntronLength > MIN_INTRON_LENGTH;
                    bool isNotOverlapped = qIntronLength > - MAX_ALIGNMENTS_OVERLAP_LENGTH;
                    bool prevStrand = candidate.candidates[prevExon].dbEndPos>candidate.candidates[prevExon].dbStartPos;
                    bool isTheSameStrand = strand==prevStrand;
                    bool isTheSameOrf = candidate.candidates[currExon].queryOrfStartPos==candidate.candidates[prevExon].queryOrfStartPos && candidate.candidates[currExon].queryOrfEndPos==candidate.candidates[prevExon].queryOrfEndPos;
                    if (isTheSameStrand && isTheSameOrf && isNotTooLongIntron && isNotTooShortIntron  && isNotOverlapped)  { //&& sameOrf
                            long bestScorePrev = dpMatrixRow[prevExon].pathScore;
                            long currScoreWithPrev = bestScorePrev + score;
                            if (currScoreWithPrev > dpMatrixRow[currExon].pathScore) {
                                dpMatrixRow[currExon].prevPotentialId = prevExon;
                                dpMatrixRow[currExon].pathScore = currScoreWithPrev;
                            } //end of if statement to update
                    } //end of if conditional statement for avoid overlap and ...
                } //end of DP 2nd for loop statement
                if (dpMatrixRow[currExon].pathScore > bestPathScore) {
                    currId = currExon;
                    bestPathScore = dpMatrixRow[currExon].pathScore;
                } //end of if conditional statement
            } //end of DP 1st for loop statement
            //end of Dynamic Progamming
            //to update <optimalExonSolution>
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
                    optimalExonSolution.emplace_back(candidate.candidates[currId]);
                    currId = dpMatrixRow[currId].prevPotentialId;
                }
                optimalExonSolution.emplace_back(candidate.candidates[currId]);
                std::sort(optimalExonSolution.begin(), optimalExonSolution.end(), Matcher::compareHitsByPosWithStrand); // Matcher::compareHitsByPosWithStrand
                candidate.score = bestPathScore;
                dpMatrixRow.clear();
            } //end of if conditional statement
        }//end of for loop statement
    }// end of function

    // trimming alignments
    std::vector<Matcher::result_t> trimExons(std::vector<Matcher::result_t> exonCandidate, unsigned int thread_idx){
        char * targetSeq = getTargetSequence(exonCandidate[0].dbKey, thread_idx);
        bool strand = exonCandidate[0].dbEndPos > exonCandidate[0].dbStartPos;
        if (std::strlen(targetSeq)==0)
            return exonCandidate;
        for (unsigned int currExon=0; currExon<exonCandidate.size(); currExon++){
//             First Exon
            if (currExon==0){
                // start codon
                int resLen = exonCandidate[currExon].qStartPos - exonCandidate[currExon].queryOrfStartPos;
                int qPos = exonCandidate[currExon].qStartPos - exonCandidate[currExon].qStartPos%3;
                int dbPos = strand? exonCandidate[currExon].dbStartPos - exonCandidate[currExon].qStartPos%3 : exonCandidate[currExon].dbStartPos + exonCandidate[currExon].qStartPos%3;
                int scope = resLen + BONUS_SCOPE;
                std::vector<splicingSiteCandidate> startCodonCands;
                int pos = -CODON_LENGTH;
                while (pos<scope){
                    int dbCurrPos = strand? dbPos-pos : dbPos+pos;
                    int qCurrPos = qPos - pos;
                    bool isMetCodon = strand? isMetCodonF(targetSeq, dbCurrPos) : isMetCodonR(targetSeq, dbCurrPos);
                    if (isMetCodon)
                        startCodonCands.emplace_back(splicingSiteCandidate(abs(qCurrPos), dbCurrPos, 0));
                    pos += 3;
                }
                if (startCodonCands.size()>0){
                    std::sort(startCodonCands.begin(), startCodonCands.end(), compareSplicingSiteCandidate); // score the less the better
                    exonCandidate[currExon].dbStartPos = startCodonCands[0].dbPos;
                    exonCandidate[currExon].qStartPos = exonCandidate[currExon].queryOrfStartPos;
                    startCodonCands.clear();
                    // Acceptor site
                } else {
                    int dbPos = exonCandidate[currExon].dbStartPos;
                    int qPos = exonCandidate[currExon].qStartPos;
                    std::vector<splicingSiteCandidate> acceptorSiteCands;
                    for (int pos = START_END_CODON_FIND_SCOPE; pos > -START_END_CODON_FIND_SCOPE; pos--){
                        int dbCurrPos = strand? dbPos+pos: dbPos-pos;
                        int qCurrPos = qPos+pos;
                        bool isAcceptorSite = strand? isAcceptorSiteF(targetSeq, dbCurrPos) : isAcceptorSiteR(targetSeq, dbCurrPos);
                        if (isAcceptorSite)
                            acceptorSiteCands.emplace_back(splicingSiteCandidate(abs(pos), dbCurrPos, qCurrPos));
                    }
                    if (acceptorSiteCands.size()>0){
                        std::sort(acceptorSiteCands.begin(), acceptorSiteCands.end(), compareSplicingSiteCandidate); // score the less the better
                        exonCandidate[currExon].dbStartPos = acceptorSiteCands[0].dbPos;
                        exonCandidate[currExon].qStartPos = acceptorSiteCands[0].qPos;
                    }
                }
                // middle Exons
            } else {
                // find donor and acceptor sites at once
                std::vector<donorAcceptorSitesCandidate> dornorAcceptorSiteCands = getDonorAcceptorSiteCands(exonCandidate[currExon - 1], exonCandidate[currExon], targetSeq);
                if (dornorAcceptorSiteCands.size()>0){
                    std::sort(dornorAcceptorSiteCands.begin(), dornorAcceptorSiteCands.end(), compareDonorAcceptorSitesCandidate); //score the more the better, dist the less the better
                    exonCandidate[currExon-1].qEndPos = dornorAcceptorSiteCands[0].qDonorPos;
                    exonCandidate[currExon-1].dbEndPos = dornorAcceptorSiteCands[0].dbDdonorPos;
                    exonCandidate[currExon].qStartPos = dornorAcceptorSiteCands[0].qAcceptorPos;
                    exonCandidate[currExon].dbStartPos = dornorAcceptorSiteCands[0].dbAcceptorPos;
                } else {
                    // donor site
                    int dbPos = exonCandidate[currExon-1].dbEndPos;
                    int qPos = exonCandidate[currExon-1].qEndPos;
                    std::vector<splicingSiteCandidate> donorSiteCands;
                    for (int pos = -START_END_CODON_FIND_SCOPE; pos < START_END_CODON_FIND_SCOPE; pos++) {
                        int dbTempPos = strand ? dbPos+pos : dbPos-pos;
                        int qTempPos = qPos + pos;
                        bool isDonorSite = strand ? isDonorSitF(targetSeq, dbTempPos) : isDonorSiteR(targetSeq, dbTempPos);
                        if (isDonorSite)
                            donorSiteCands.emplace_back(splicingSiteCandidate(abs(pos), dbTempPos, qTempPos));
                    }
                    if (donorSiteCands.size()>0){
                        std::sort(donorSiteCands.begin(), donorSiteCands.end(),compareSplicingSiteCandidate);// score the less the better
                        exonCandidate[currExon-1].dbEndPos = donorSiteCands[0].dbPos;
                        exonCandidate[currExon-1].qEndPos = donorSiteCands[0].qPos;
                    }
                    // acceptor site
                    dbPos = exonCandidate[currExon].dbStartPos;
                    qPos = exonCandidate[currExon].qStartPos;
                    std::vector<splicingSiteCandidate> acceptorSiteCands;
                    for (int pos=START_END_CODON_FIND_SCOPE; pos > -START_END_CODON_FIND_SCOPE; pos--) {
                        int dbTempPos = strand ? dbPos+pos : dbPos-pos;
                        int qTempPos = qPos + pos;
                        bool isAcceptorSite = strand ? isAcceptorSiteF(targetSeq, dbTempPos) : isAcceptorSiteR(targetSeq, dbTempPos);
                        if (isAcceptorSite)
                            acceptorSiteCands.emplace_back(abs(pos), dbTempPos, qTempPos);
                    }
                    if (acceptorSiteCands.size()>0){
                        std::sort(acceptorSiteCands.begin(), acceptorSiteCands.end(), compareSplicingSiteCandidate); // score the less the better
                        exonCandidate[currExon].dbStartPos = acceptorSiteCands[0].dbPos;
                        exonCandidate[currExon].qStartPos = acceptorSiteCands[0].qPos;
                    }
                }
            }
        } // for End
        // Last Exon
        // stop codon
        int resLen = exonCandidate[exonCandidate.size()-1].queryOrfEndPos - exonCandidate[exonCandidate.size()-1].qEndPos;
        int scope = resLen + BONUS_SCOPE;
        int dbPos = strand ? exonCandidate[exonCandidate.size()-1].dbEndPos-exonCandidate[exonCandidate.size()-1].qEndPos%3-1 : exonCandidate[exonCandidate.size()-1].dbEndPos+exonCandidate[exonCandidate.size()-1].qEndPos%3+1;
        bool doFindStpCodon = false;
        while (scope>0){
            bool isStpCodon = strand ? isStpCodonF(targetSeq, dbPos) : isStpCodonR(targetSeq, dbPos);
            if (isStpCodon) {
                exonCandidate[exonCandidate.size()-1].qEndPos = exonCandidate[exonCandidate.size()-1].queryOrfEndPos;
                exonCandidate[exonCandidate.size()-1].dbEndPos =  strand ? dbPos+CODON_LENGTH : dbPos-CODON_LENGTH; // dbPos;
                doFindStpCodon = true;
                break;
            }
            scope -= 3;
            dbPos = strand ? dbPos+3 : dbPos-3;
        }
        // DonorSite
        if (!doFindStpCodon){
            int dbPos = exonCandidate[exonCandidate.size()-1].dbEndPos;
            int qPos = exonCandidate[exonCandidate.size()-1].qEndPos;
            std::vector<splicingSiteCandidate> donorSiteCands;
            for (int pos = -START_END_CODON_FIND_SCOPE; pos < START_END_CODON_FIND_SCOPE; pos++){
                int dbTempPos = strand? dbPos+pos : dbPos-pos;
                int qTempPos = qPos+pos;
                bool isDonorSite = strand ? isDonorSitF(targetSeq, dbTempPos) : isDonorSiteR(targetSeq, dbTempPos);
                if (isDonorSite){
                    donorSiteCands.emplace_back(splicingSiteCandidate(abs(pos), dbTempPos, qTempPos));
                }
            }
            if (donorSiteCands.size()>0) {
                std::sort(donorSiteCands.begin(), donorSiteCands.end(), compareSplicingSiteCandidate);
                exonCandidate[exonCandidate.size()-1].dbEndPos = donorSiteCands[0].dbPos;
                exonCandidate[exonCandidate.size()-1].qEndPos = donorSiteCands[0].qPos;
            }
        }
        return exonCandidate;
    } // method End

private:
    // class variable
    const int MAX_INTRON_LENGTH = 200000;
    const int MIN_INTRON_LENGTH = 30;
    const int CODON_LENGTH = 3;
    const int BONUS_SCOPE = 60;
    const int MAX_RESIDUE_LENGTH = 30;
    const int START_END_CODON_FIND_SCOPE = 90;
    const int MAX_ALIGNMENTS_OVERLAP_LENGTH = 90;

    IndexReader * tDbr;
    // Do I need?
//    IndexReader * qDbr;
    // Do I need?
//    unsigned int queryKey;
    std::vector<DpMatrixRow> dpMatrixRow;
//    // Do I need?
//    typedef std::pair<char, int> cigarTuple;
//
//    // Do I need?
//    Matcher::result_t stpCodonExtension(Matcher::result_t inputExon){
//        size_t lenSTPCodon = inputExon.dbStartPos < inputExon.dbEndPos ? 3 : -3;
//        bool haveSTPCodon = inputExon.qEndPos == inputExon.queryOrfEndPos;
//        inputExon.dbEndPos = haveSTPCodon ? (inputExon.dbEndPos + lenSTPCodon) : inputExon.dbEndPos;
//        return inputExon;
//    }

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
            bool tooLongIntron = distBetweenExons > MAX_INTRON_LENGTH;
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

//    // Do I need?
//    std::vector<cigarTuple> cigarToTuple(std::string backtrace){
//        std::vector<cigarTuple> returnTuple;
//        size_t tmpPos = 0;
//        int cnt = 0;
//        for (size_t pos = 0; pos < backtrace.size(); pos++) {
//            if(!isdigit(backtrace[pos])){
//                int num = std::stoi(backtrace.substr(tmpPos, cnt));
//                char cha = backtrace[pos];
//                returnTuple.emplace_back(std::pair<char, int>(cha,num));
//                tmpPos = pos+1;
//                cnt=0;
//            } else{
//                cnt++;
//            }
//        }
//        return returnTuple;
//    }

    // to find Donnor Sites and Acceptor sites
    bool isAcceptorSiteF(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index-2]);
        char nt2 = std::toupper(targetSeq[index-1]);
        return nt1=='A'&&nt2=='G';
    }
    bool isDonorSitF(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index+1]);
        char nt2 = std::toupper(targetSeq[index+2]);
        return  (nt1=='G'&&nt2=='T') ; //|| (nt1=='G'&&nt2=='C');
    }
    bool isAcceptorSiteR(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index+1]);
        char nt2 = std::toupper(targetSeq[index+2]);
        return  nt1=='C' && nt2=='T';
    }
    bool isDonorSiteR(char * targetSeq, int index){
        char nt1 = std::toupper(targetSeq[index-2]);
        char nt2 = std::toupper(targetSeq[index-1]);
        return (nt1=='A' && nt2=='C'); // || (nt1=='G' && nt2=='C');
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
    char * getTargetSequence(unsigned int dbKey, unsigned int thread_idx){
        size_t targetId = tDbr->sequenceReader->getId(dbKey);
        char * targetSeq = tDbr->sequenceReader->getData(targetId, thread_idx);
        return targetSeq;
    }

//    // Do I need?
//    char * querySequence(int thread_idx){
//        unsigned int queryId = qDbr->sequenceReader->getId(queryKey);
//        char * qSeq = qDbr->sequenceReader->getData(queryId, thread_idx);
//        return qSeq;
//    }
//
//    // Do I need?
//    std::pair<std::string, int> cigarQueryPosUpdateAcceptorSite(std::string cigar, int overlap){
//        int returnNumber = 0;
//        std::string returnString;
//        std::vector<cigarTuple> tupleVector = cigarToTuple(cigar);
//        for(size_t cnt=0; cnt<tupleVector.size(); cnt++){
//            int tempOverlap = std::min(overlap,tupleVector[cnt].second);
//            switch (tupleVector[cnt].first) {
//                case 'I':
//                    returnNumber += tempOverlap;
//                    tupleVector[cnt].second -= tempOverlap;
//                    break;
//                case  'D':
//                    tupleVector[cnt].second -= tempOverlap;
//                    overlap -= tempOverlap;
//                    break;
//                default:
//                    tupleVector[cnt].second -= tempOverlap;
//                    returnNumber += tempOverlap;
//                    overlap -= tempOverlap;
//                    break;
//            }
//            if(overlap==0)
//                break;
//        }
//        for(size_t cnt=0; cnt<tupleVector.size(); cnt++){
//            if(tupleVector[cnt].second==0)
//                continue;
//
//            if(returnString == "") {
//                switch (tupleVector[cnt].first) {
//                    case 'I':
//                        returnNumber += tupleVector[cnt].second;
//                        break;
//                    case 'D':
//                        overlap -= tupleVector[cnt].second;
//                        break;
//                    default:
//                        returnString = std::to_string(tupleVector[cnt].second) + tupleVector[cnt].first + returnString;
//                        break;
//                }
//            } else {
//                returnString = std::to_string(tupleVector[cnt].second) + tupleVector[cnt].first + returnString;
//            }
//        }
//        return std::pair<std::string,int>(returnString, returnNumber);
//    }
//    // Do I need?
//    std::pair<std::string , int>  cigarQueryPosUpdateDonorSite(std::string cigar, int overlap){
//        int returnNumber = 0;
//        std::string returnString;
//        std::vector<cigarTuple> tupleVector = cigarToTuple(cigar);
//        for(size_t cnt=0; cnt<tupleVector.size(); cnt++){
//            int tempOverlap = std::min(overlap,tupleVector[cnt].second);
//            switch (tupleVector[cnt].first) {
//                case 'I':
//                    returnNumber += tempOverlap;
//                    tupleVector[cnt].second -= tempOverlap;
//                    break;
//                case  'D':
//                    tupleVector[cnt].second -= tempOverlap;
//                    overlap -= tempOverlap;
//                    break;
//                default:
//                    tupleVector[cnt].second -= tempOverlap;
//                    returnNumber += tempOverlap;
//                    overlap -= tempOverlap;
//                    break;
//            }
//            if(overlap==0)
//                break;
//        }
//        for(size_t cnt=0; cnt<tupleVector.size(); cnt++){
//            if(tupleVector[cnt].second==0)
//                continue;
//
//            if(returnString == "") {
//                switch (tupleVector[cnt].first) {
//                    case 'I':
//                        returnNumber += tupleVector[cnt].second;
//                        break;
//                    case 'D':
//                        overlap -= tupleVector[cnt].second;
//                        break;
//                    default:
//                        returnString = std::to_string(tupleVector[cnt].second) + tupleVector[cnt].first + returnString;
//                        break;
//                }
//            } else {
//                returnString = std::to_string(tupleVector[cnt].second) + tupleVector[cnt].first + returnString;
//            }
//        }
//
//        return std::pair<std::string,int>(returnString, returnNumber);
//    }
//    // Do I need?
//    std::string addCigar(std::string cigar, char symbol, int length){
//        std::string newCigar = "";
//        std::vector<cigarTuple> tupleVector = cigarToTuple(cigar);
//        if (cigar != "" && tupleVector[tupleVector.size()-1].first == symbol) {
//            tupleVector[tupleVector.size() - 1].second += length;
//        } else if (length > 0) {
//            tupleVector.emplace_back(std::pair<char, int>(symbol, length));
//        }
//        for(size_t cnt=0; cnt<tupleVector.size(); cnt++){
//            newCigar = newCigar + std::to_string(tupleVector[cnt].second) + tupleVector[cnt].first;
//        }
//        return newCigar;
//    }
//
//    // Do I need?
//    //to update identity
//    int cigarLength(std::string cigar){
//        std::vector<cigarTuple> cigarTupleVec =  cigarToTuple(cigar);
//        int returnNum = 0;
//        for (size_t i=0; i < cigarTupleVec.size(); i++){
//            returnNum += cigarTupleVec[i].second;
//        }
//        return  returnNum;
//    }
//    // Do I need?
//    float matchRatio(std::string cigar){
//        std::vector<cigarTuple> cigarTupleVec =  cigarToTuple(cigar);
//        int returnNum1 = 0;
//        int returnNum2 = 0;
//        for (size_t i=0; i < cigarTupleVec.size(); i++){
//            returnNum2 += cigarTupleVec[i].second;
//            if(cigarTupleVec[i].first == 'M')
//                returnNum1 += cigarTupleVec[i].second;
//        }
//        return  (float)returnNum1/returnNum2;
//    }

    // Do I need?
    int queryOrfLength(Matcher::result_t exon){
        return  exon.queryOrfEndPos - exon.queryOrfStartPos +1;
    }
    int queryLength(Matcher::result_t exon){
        return  exon.qEndPos - exon.qStartPos +1;
    }
    // Do I need?
    int dbLength(Matcher::result_t exon){
        return  abs(exon.dbEndPos - exon.dbStartPos) +1;
    }
//    // Do I need?
//    bool metContainF(char * qSeq, int qStartPos, int qEndPos){
//        int currPos = qStartPos;
//        while(currPos<qEndPos-2){
//            char nt1 = std::toupper(qSeq[currPos]);
//            char nt2 = std::toupper(qSeq[currPos+1]);
//            char nt3 = std::toupper(qSeq[currPos+2]);
//            if(nt1=='A' && nt2=='T' && nt3=='G'){
//                return true;
//            }
//            currPos += 3;
//        }
//        return false;
//    }
//    // Do I need?
//    bool metContainR(char * qSeq, int qStartPos, int qEndPos){
//        int currPos = qStartPos;
//        while(currPos<qEndPos-2){
//            char nt1 = std::toupper(qSeq[currPos]);
//            char nt2 = std::toupper(qSeq[currPos+1]);
//            char nt3 = std::toupper(qSeq[currPos+2]);
//            if(nt1=='T' && nt2=='A' && nt3=='C'){
//                return true;
//            }
//            currPos += 3;
//        }
//        return false;
//    }
//    // Do I need?
//    bool firstExon(int qStartPos, int qOrfStartPos, int inScope, int outScope){
//        return (qStartPos - outScope) < qOrfStartPos && (qOrfStartPos < qStartPos + inScope);
//    }
//    // Do I need?
//    bool lastExon(int qEndPos, int qOrfEndPos, int inScope, int outScope){
//        return (qEndPos-inScope < qOrfEndPos) && (qOrfEndPos < qEndPos + outScope);
//    }

    std::vector<donorAcceptorSitesCandidate> getDonorAcceptorSiteCands(Matcher::result_t prevExon, Matcher::result_t currExon, char * targetSeq){
        std::vector<donorAcceptorSitesCandidate> dornorAcceptorSiteCands;
        std::vector<splicingSiteCandidate> donorSiteCands;
        std::vector<splicingSiteCandidate> acceptorSiteCands;
        bool strand = currExon.dbEndPos > currExon.dbStartPos;
        int residueLength = currExon.qStartPos - prevExon.qEndPos - 1;
        if (abs(residueLength) > MAX_RESIDUE_LENGTH)
            return dornorAcceptorSiteCands;
        int dbPrevPos = prevExon.dbEndPos;
        int dbCurrPos = currExon.dbStartPos;
        int qPrevPos = prevExon.qEndPos;
        int qCurrPos = currExon.qStartPos;
        // donor
        int loopStartPos = residueLength < 0 ? -CODON_LENGTH + residueLength : -CODON_LENGTH;
        int loopEndPos = residueLength < 0 ? CODON_LENGTH + 1 : CODON_LENGTH + 1 + residueLength;
        for (int pos = loopStartPos; pos<loopEndPos; pos++) {
            int dbTempPos = strand? dbPrevPos+pos: dbPrevPos-pos;
            int qTempPos = qPrevPos + pos;
            bool isDonorSite = strand? isDonorSitF(targetSeq, dbTempPos) : isDonorSiteR(targetSeq, dbTempPos);
            if (isDonorSite)
                donorSiteCands.emplace_back(splicingSiteCandidate(0, dbTempPos, qTempPos));
        }
        // acceptor
        loopStartPos = residueLength > 0 ? -CODON_LENGTH - residueLength : -CODON_LENGTH;
        loopEndPos = residueLength > 0 ? CODON_LENGTH + 1 : CODON_LENGTH + 1 - residueLength;
        for(int pos = loopStartPos; pos < loopEndPos; pos++) {
            int dbTempPos = strand? dbCurrPos+pos: dbCurrPos-pos;
            int qTempPos = qCurrPos + pos;
            bool isAcceptorSite = strand ? isAcceptorSiteF(targetSeq, dbTempPos) : isAcceptorSiteR(targetSeq, dbTempPos);
            if (isAcceptorSite)
                acceptorSiteCands.emplace_back(splicingSiteCandidate(0, dbTempPos, qTempPos));
        }
        for (unsigned int donorCand = 0; donorCand<donorSiteCands.size(); donorCand++){
            for (unsigned int acceptorCand = 0; acceptorCand<acceptorSiteCands.size(); acceptorCand++) {
                int qDonorSiteCandPos = donorSiteCands[donorCand].qPos;
                int qAcceptorSiteCandPos = acceptorSiteCands[acceptorCand].qPos;
                int dbDonorSiteCandPos = donorSiteCands[donorCand].dbPos;
                int dbAcceptorSiteCandPos = acceptorSiteCands[acceptorCand].dbPos;
                int score;
                int dist = qAcceptorSiteCandPos - qDonorSiteCandPos;
                if (dist == 1){
                    score=3;
                } else if (dist%3 == 1){
                    score = 2;
                } else if (dist>0) {
                    score = 0;
                } else {
                    score = -1;
                }
                dornorAcceptorSiteCands.emplace_back(donorAcceptorSitesCandidate(score, dist,  dbDonorSiteCandPos, dbAcceptorSiteCandPos, qDonorSiteCandPos, qAcceptorSiteCandPos));
            }
        }
        return dornorAcceptorSiteCands;
    }
};//end of class


int findexons(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    // TODO
    float falsePositiveFilteringRatio = (float)par.filteringRatio/100;
    //TODO
    bool encludingStopCodon = (bool)par.encludingStopCodon==1;

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
        std::vector<Matcher::result_t> totalInputAlignments;
        std::vector<Matcher::result_t> orfAlignments;
        std::vector<Matcher::result_t> optimalExonSolution;
        std::vector<ExonCandidates> optimalSolutionsWithScores;
        char buffer[2048];

#pragma omp  for schedule(dynamic, 10)
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            progress.updateProgress();
            const unsigned int queryKey = alnDbr.getDbKey(i);
//            if (queryKey==345)
//                std::cout<<""<<std::endl;
//            ExonFinder exonFinder(&tDbr, &qDbr, queryKey);
            ExonFinder exonFinder(&tDbr);
            char *data = alnDbr.getData(i, thread_idx);
            if(data[0]=='\0'){
                resultWriter.writeData("",0,queryKey,thread_idx);
                continue;
            }
            totalInputAlignments.clear();
            Matcher::readAlignmentResults(totalInputAlignments, data, true);
            std::sort(totalInputAlignments.begin(), totalInputAlignments.end(), Matcher::compareOrfStartOrfEnd);
            int prevQueryOrfStartPos = totalInputAlignments[0].queryOrfStartPos;
            int prevQueryOrfEndPos = totalInputAlignments[0].queryOrfEndPos;
            resultWriter.writeStart(thread_idx);
            long orfScore = 0;
            long maxScore = 0;
            optimalSolutionsWithScores.clear();
            for(size_t resIdx = 0; resIdx < totalInputAlignments.size(); resIdx++){
                // to make sure that query position is never reversed, only db position can be reversed
                // In default we search only on the formward frame 1,2,3 so this function is not called
                // It is only important if search also on the backward frame!
                if(totalInputAlignments[resIdx].qStartPos > totalInputAlignments[resIdx].qEndPos){
                    totalInputAlignments[resIdx] = exonFinder.flipExons(totalInputAlignments[resIdx]);
                } // end of if conditional statement to correct flipped exon
                bool querySameOrf = prevQueryOrfStartPos == totalInputAlignments[resIdx].queryOrfStartPos && prevQueryOrfEndPos == totalInputAlignments[resIdx].queryOrfEndPos;
                if(querySameOrf){
                    orfAlignments.emplace_back(totalInputAlignments[resIdx]);
                }else{
                    // getOptimalAlignmentsCombinaitons
                    exonFinder.findOptimalExons(
                            optimalExonSolution,
                            orfAlignments,
//                            thread_idx,
                            orfScore);
                    orfAlignments.clear();
                    if(orfScore>maxScore){
                        optimalSolutionsWithScores.emplace_back(ExonCandidates(orfScore, optimalExonSolution));
                        maxScore = orfScore;
                    }
                    orfAlignments.emplace_back(totalInputAlignments[resIdx]);
                    prevQueryOrfStartPos = totalInputAlignments[resIdx].queryOrfStartPos;
                    prevQueryOrfEndPos = totalInputAlignments[resIdx].queryOrfEndPos;
                }
            }
            //last orf info -> optimal
            if(orfAlignments.size() > 0){
                // getOptimalAlignmentsCombinaitons
                exonFinder.findOptimalExons(
                        optimalExonSolution,
                        orfAlignments,
//                        thread_idx,
                        orfScore
                        );
                orfAlignments.clear();
                if(orfScore>maxScore){
                    optimalSolutionsWithScores.emplace_back(ExonCandidates(orfScore, optimalExonSolution));
                    maxScore = orfScore;
                }
            }
            // output
            if(optimalSolutionsWithScores.size() > 0) {
                // the last vector in optimalSolutionsWithScores contains the best solution
                optimalSolutionsWithScores[optimalSolutionsWithScores.size() - 1].candidates = exonFinder.trimExons(optimalSolutionsWithScores[optimalSolutionsWithScores.size() - 1].candidates, thread_idx);
                size_t optimestPathSize = optimalSolutionsWithScores[optimalSolutionsWithScores.size() - 1].candidates.size();
                for(size_t optIdx = 0; optIdx < optimestPathSize; optIdx++){
                    size_t len = Matcher::resultToBuffer(buffer, optimalSolutionsWithScores[optimalSolutionsWithScores.size() - 1].candidates[optIdx], true, false, true);
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




