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
    ExonFinder(IndexReader * tDbr, IndexReader * qDbr, unsigned int queryKey) : tDbr(tDbr), qDbr(qDbr), queryKey(queryKey){}

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
            unsigned int thread_idx,
            long & orfScore
    ) {
        std::sort(alignments.begin(), alignments.end(), Matcher::compareByDbkeyAndStrand);
        std::vector<ExonCandidates> scoresAndAlignmentVectors = createPotentialExonCombinations(alignments);
        for(size_t candidateIdx = 0; candidateIdx < scoresAndAlignmentVectors.size(); candidateIdx++){
            ExonCandidates & scoreAndAlignmentVector = scoresAndAlignmentVectors[candidateIdx];
            dpMatrixRow.clear();
            for (size_t id = 0; id < scoreAndAlignmentVector.candidates.size(); id++) {
                long score = queryLength(scoreAndAlignmentVector.candidates[id]) * scoreAndAlignmentVector.candidates[id].seqId;
                dpMatrixRow.emplace_back(DpMatrixRow(id,score));
            }
            long bestPathScore = INT_MIN;
            size_t currId;
            for (size_t currExon = 0; currExon < scoreAndAlignmentVector.candidates.size(); currExon++) {
                long score = queryLength(scoreAndAlignmentVector.candidates[currExon]) * scoreAndAlignmentVector.candidates[currExon].seqId;
                bool strand = scoreAndAlignmentVector.candidates[currExon].dbEndPos > scoreAndAlignmentVector.candidates[currExon].dbStartPos;
                for (size_t prevExon = 0; prevExon < currExon; prevExon++) {
                    int tIntronLength = strand ? scoreAndAlignmentVector.candidates[currExon].dbStartPos - scoreAndAlignmentVector.candidates[prevExon].dbEndPos + 1 : scoreAndAlignmentVector.candidates[prevExon].dbEndPos - scoreAndAlignmentVector.candidates[currExon].dbStartPos + 1 ;
                    int qIntronLength = scoreAndAlignmentVector.candidates[currExon].qStartPos - scoreAndAlignmentVector.candidates[prevExon].qEndPos + 1;
                    bool isNotTooLongIntron = (tIntronLength < MAX_INTRON_LENGTH);
                    bool isNotTooShortIntron = tIntronLength > MIN_INTRON_LENGTH;
                    bool isNotOverlapped = qIntronLength > - MAX_ALIGNMENTS_OVERLAP_LENGTH;
                    bool prevStrand = scoreAndAlignmentVector.candidates[prevExon].dbEndPos > scoreAndAlignmentVector.candidates[prevExon].dbStartPos;
                    bool isTheSameStrand = strand==prevStrand;
                    bool isTheSameOrf = scoreAndAlignmentVector.candidates[currExon].queryOrfStartPos == scoreAndAlignmentVector.candidates[prevExon].queryOrfStartPos && scoreAndAlignmentVector.candidates[currExon].queryOrfEndPos == scoreAndAlignmentVector.candidates[prevExon].queryOrfEndPos;
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
            for(size_t i=0; i < scoresAndAlignmentVectors.size(); i++){
                if(bestPathScore > scoresAndAlignmentVectors[i].score) {
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
                    optimalExonSolution.emplace_back(scoreAndAlignmentVector.candidates[currId]);
                    currId = dpMatrixRow[currId].prevPotentialId;
                }
                optimalExonSolution.emplace_back(scoreAndAlignmentVector.candidates[currId]);
                std::sort(optimalExonSolution.begin(), optimalExonSolution.end(), Matcher::compareHitsByPosWithStrand); // Matcher::compareHitsByPosWithStrand
                scoreAndAlignmentVector.score = bestPathScore;
                dpMatrixRow.clear();
                // 2nd aligner
                std::string targetSeq = std::string(getTargetSequence(scoreAndAlignmentVector.candidates[0].dbKey, thread_idx));
                std::string querySeq = std::string(getQuerySequence(queryKey, thread_idx));
                int s1 = getBlosum62Score(aminoAcidPair('M', 'M'));
                bool strand = optimalExonSolution[0].dbEndPos > optimalExonSolution[0].dbStartPos;
                int prevQueryEnd = optimalExonSolution[0].queryOrfStartPos -1;
                int finalQueryEnd = optimalExonSolution[0].queryOrfEndPos;
                int tScope1 = optimalExonSolution[0].qStartPos - (optimalExonSolution[0].queryOrfStartPos-1) + BONUS_SCOPE;
                int tScope2 = optimalExonSolution[0].queryOrfEndPos - optimalExonSolution[0].qEndPos + BONUS_SCOPE;
                int prevTargetEnd = optimalExonSolution[0].dbStartPos + (strand ? -tScope1 : +tScope1);
                int finalTargetEnd = optimalExonSolution[optimalExonSolution.size() - 1].dbEndPos + (strand ? +tScope2 : -tScope2);
                std::vector<std::string> queryMissingAlnsScopeSeqs;
                std::vector<std::string> targetMissingAlnsScopeSeqs;
                for (size_t i=0; i<optimalExonSolution.size();i++){
                    int currQueryStart = optimalExonSolution[i].qStartPos;
                    int currTargetStart = optimalExonSolution[i].dbStartPos;
                    if (currQueryStart - prevQueryEnd > 90) {
                        int qStart = prevQueryEnd + 1;
                        int qEnd = currQueryStart - 1;
                        int tStart = strand ? prevTargetEnd + 1 : prevTargetEnd - 1;
                        int tEnd = strand ? currTargetStart - 1 : currTargetStart + 1;
                        std::string querySubSeq = getSubSequence(querySeq, qStart, qEnd + 1);
                        std::string targetSubSeq = getSubSequence(targetSeq, tStart, strand ? tEnd + 1 : tEnd - 1);
                        std::vector<Matcher::result_t> foundMissingAlnCandidates = findMissingAlignments(
                                querySubSeq, getOrfs(targetSubSeq), optimalExonSolution[0].dbKey,
                                qStart, tStart,
                                optimalExonSolution[0].queryOrfStartPos, optimalExonSolution[0].queryOrfEndPos,
                                tStart, tEnd
                        );
                    }
                    prevQueryEnd = optimalExonSolution[i].qEndPos;
                    prevTargetEnd = optimalExonSolution[i].dbEndPos;
                }
                if (finalQueryEnd - prevQueryEnd > 90){
                    int qStart = prevQueryEnd+1;
                    int tStart =strand ? prevTargetEnd+1 : prevTargetEnd-1;
                    std::string querySubSeq = getSubSequence(querySeq, qStart, finalQueryEnd+1);
                    std::string targetSubSeq = getSubSequence(targetSeq, tStart, strand? finalTargetEnd+1:finalTargetEnd-1);
                    std::vector<std::string> qOrfVector = getOrfs(querySubSeq);
                    std::vector<std::string> dbOrfVector = getOrfs(targetSubSeq);

                    int x = querySubSeq.size()%3;
                    int y = targetSubSeq.size()%3;

                    int foo = 0;
                }
//                std::string str = "ATGCATGCATGCATGC";
//                std::string seq1 = getSubSequence(str,3,7);
//                std::string seq2 = getSubSequence(str,7,3);
//                std::string seq3 = getSubSequence(str,6,2);
//                std::string seq4 = getSubSequence(str,7,3);
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
                        // TEMP
                        startCodonCands.emplace_back(splicingSiteCandidate(abs(qCurrPos), dbCurrPos, 0));
//                        startCodonCands.emplace_back(splicingSiteCandidate(abs(pos), dbCurrPos, 0));
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
    IndexReader * qDbr;
    unsigned int queryKey;
    std::vector<DpMatrixRow> dpMatrixRow;
    typedef std::pair<char, int> cigarTuple;
    typedef std::pair<char, char> aminoAcidPair;

    std::string getSubSequence(std::string seq, size_t start, size_t end){
        if (end >= start) {
            size_t len = end - start; //+1
            return seq.substr(start, len);
        } else{
            size_t len = start - end; //+1
            std::string tempSubseq = seq.substr(end+1,len);
            std::string subSeq;
            for (size_t i=0; i<tempSubseq.size(); i++){
//                char nt = tempSubseq.substr(i,1);
                char nt = tempSubseq[i];
                if (toupper(nt)=='A'){
                    subSeq = "T"+subSeq;
                } else if(toupper(nt)=='T'){
                    subSeq = "A"+subSeq;
                } else if (toupper(nt)=='G'){
                    subSeq = "C"+subSeq;
                } else if (toupper(nt)=='C'){
                    subSeq = "G"+subSeq;
                } else {
                    subSeq = nt + subSeq;
                }
            }
            return subSeq;
        }
    }

    std::vector<std::string> getOrfs(std::string seq){
        std::vector<std::string> orfVector;
        size_t seqLen = seq.size();
        size_t len = seq.size();
        for (size_t i=0; i<3; i++){
            orfVector.emplace_back(getSubSequence(seq,i,seqLen-(seqLen-i)%3));
        }
        return orfVector;
    }

    std::string getCigarString(std::vector<cigarTuple> cList){
        std::string cigar;
        char prevCigar = cList[0].second;
        int prevNum = cList[0].first;
        for (size_t i=1; i<cList.size(); i++) {
            int currNum = cList[i].first;
            char currCigar = toupper(cList[i].second);
            if (currCigar != prevCigar){
                cigar = cigar + std::to_string(prevNum) + prevCigar;
                prevNum = currNum;
                prevCigar = currCigar;
            } else {
                prevNum += currNum;
            }
        }
        cigar = cigar + std::to_string(prevNum) + prevCigar;
        return cigar;
    }

    int getGapOpen(std::string cigarString) {
        int numGaps;
        for (size_t i=0; i<cigarString.size(); i++){
            char cigar = cigarString[i];
            if (cigar=='I' && cigar=='D')
                numGaps++;
        }
        return numGaps;
    }

    std::vector<Matcher::result_t> findMissingAlignments(
            std::string qOrfSeq,
            std::vector<std::string> tOrfVector,
            unsigned int dbkey,
            int qAlnStart,
            int dbStart,
            int queryOrfStartPos,
            int queryOrfEndPos,
            int dbOrfStartPos,
            int dbOrfEndPos
    ){
        std::vector<Matcher::result_t> foundMissingAlignments;
        int penalty = -4;
        for (size_t tOrfIdx = 0; tOrfIdx < 3; tOrfIdx++){
            std::string tOrfSeq = tOrfVector[tOrfIdx];
            std::map<std::pair<int,int>, int> dpMatrix;
            dpMatrix.insert({std::pair<int,int>(-3,-3), 0});
            int qPos = 0;
            int tPos = 0;

            std::vector<std::pair<int,int>> foundAlnsEndPositions;
            while (qPos < qOrfSeq.size()){
                dpMatrix.insert({std::pair<int,int>(qPos, -3), 0});
                qPos += 3;
            }
            while (tPos < tOrfSeq.size()){
                dpMatrix.insert({std::pair<int,int>(-3, tPos), 0});
                tPos += 3;
            }
            qPos = 0;
            int maxScore = -1;
            while (qPos < qOrfSeq.size()) {
                std::string qCodon = getSubSequence(qOrfSeq, qPos, qPos + CODON_LENGTH);
                char qAA = translateCodon(qCodon);
                tPos = 0;
                while (tPos < tOrfSeq.size()) {
                    std::string tCodon = getSubSequence(tOrfSeq, tPos, tPos + CODON_LENGTH);
                    char tAA = translateCodon(tCodon);
                    int score = std::max(
                            {
                                    dpMatrix.at(std::pair<int, int>(qPos - 3, tPos - 3)) +
                                    getBlosum62Score(aminoAcidPair(qAA, tAA)),
                                    dpMatrix.at(std::pair<int, int>(qPos - 3, tPos)) + penalty,
                                    dpMatrix.at(std::pair<int, int>(qPos, tPos - 3)) + penalty,
                                    0
                            }
                    );
                    dpMatrix.insert({std::pair<int,int>(qPos, tPos), score});
                    tPos += 3;
                    if (score > maxScore){
                        maxScore = score;
                        foundAlnsEndPositions.clear();
                        foundAlnsEndPositions.emplace_back(std::pair<int,int>(qPos, tPos));
                    } else if (score == maxScore){
                        foundAlnsEndPositions.emplace_back(std::pair<int,int>(qPos, tPos));
                    }
                }
                qPos += 3;
            } // DP End
            for (size_t endPosPair=0; endPosPair<foundAlnsEndPositions.size(); endPosPair++) {
                int qDpEndPos = foundAlnsEndPositions[endPosPair].first;
                int tDpEndPos = foundAlnsEndPositions[endPosPair].second;
                int qPrevPos = qDpEndPos;
                int tPrevPos = tDpEndPos;
                int matches;
                int mismarches;
                int gapOpen;
                unsigned int alnLen;
                qPos = qDpEndPos;
                tPos = tDpEndPos;
                std::vector<cigarTuple> cigarVector;
                float seqId;
                std::string cigarString;

                while(true){
                    std::string qCodon = getSubSequence(qOrfSeq, qPos, qPos + CODON_LENGTH);
                    char qAA = translateCodon(qCodon);
                    std::string tCodon = getSubSequence(tOrfSeq, tPos, tPos + CODON_LENGTH);
                    char tAA = translateCodon(tCodon);
                    int currScore = dpMatrix[std::pair<int,int>(qPos, tPos)];

                    if (currScore == 0) {
                        if (alnLen == 0)
                            break;
                        seqId = (float) matches / alnLen;
                        cigarString = getCigarString(cigarVector);
                        gapOpen = getGapOpen(cigarString);
                        int qStartPos = qPrevPos + qAlnStart;
                        int qEndPos = qDpEndPos + qAlnStart;
                        unsigned int qLen = qEndPos - qStartPos + 1;
                        int dbStartPos = tPrevPos + dbStart;
                        int dbEndPos = tDpEndPos + dbStart;
                        unsigned int dbLen = std::abs(dbStartPos - dbEndPos) + 1;
                        float qCov = (float)qLen/(queryOrfEndPos - queryOrfStartPos+1);
                        foundMissingAlignments.emplace_back(
                                Matcher::result_t(
                                        dbkey,
                                        0, //score
                                        qCov,//qCov
                                        0,//dbCov
                                        seqId,
                                        10000, //eValue = 1.0e+04
                                        alnLen,
                                        qStartPos,
                                        qEndPos,
                                        qLen,
                                        dbStartPos,
                                        dbEndPos,
                                        dbLen,
                                        queryOrfStartPos,
                                        queryOrfEndPos,
                                        dbOrfStartPos,
                                        dbOrfEndPos,
                                        cigarString
                                )
                        );
                        break;
                    } else if (currScore == dpMatrix[std::pair<int,int>(qPos-CODON_LENGTH, tPos-CODON_LENGTH)]+getBlosum62Score(aminoAcidPair(qAA, tAA))) {
                        alnLen += CODON_LENGTH;
                        matches = (qAA == tAA) ? matches + CODON_LENGTH : matches;
                        mismarches = (qAA == tAA) ? mismarches + CODON_LENGTH : mismarches;
                        cigarVector.emplace_back(cigarTuple(CODON_LENGTH, 'M'));
                        qPrevPos = qPos;
                        tPrevPos = tPos;
                        qPos -= 3;
                        tPos -= 3;
                    } else if (currScore == dpMatrix[std::pair<int,int>(qPos-3, tPos)]+penalty){
                        alnLen += CODON_LENGTH;
                        mismarches = mismarches+CODON_LENGTH;
                        cigarVector.emplace_back(cigarTuple(CODON_LENGTH, 'D'));
                        qPrevPos = qPos;
                        qPos-=3;
                    } else if (currScore == dpMatrix[std::pair<int,int>(qPos, tPos-3)]+penalty) {
                        alnLen += CODON_LENGTH;
                        mismarches = mismarches + CODON_LENGTH;
                        cigarVector.emplace_back(cigarTuple(CODON_LENGTH, 'I'));
                        tPrevPos = tPos;
                        tPos -= 3;
                    }

                } // backtracing
            } // for end iterating db results
        } // for end iterating db orfs
        return foundMissingAlignments;
    }

//
//    // Do I need?
//    Matcher::result_t stpCodonExtension(Matcher::result_t inputExon){
//        size_t lenSTPCodon = inputExon.dbStartPos < inputExon.dbEndPos ? 3 : -3;
//        bool haveSTPCodon = inputExon.qEndPos == inputExon.queryOrfEndPos;
//        inputExon.dbEndPos = haveSTPCodon ? (inputExon.dbEndPos + lenSTPCodon) : inputExon.dbEndPos;
//        return inputExon;
//    }
//    std::map<char, char> complementaryNT = {
//            {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}
//    };

    std::map<std::string, char> codonTable = {
            {"TCA",'S'}, {"TCC",'S'}, {"TCG",'S'}, {"TCT",'S'}, {"AGC",'S'}, {"AGT",'S'},
            {"TTC",'F'}, {"TTT",'F'},
            {"TTA",'L'}, {"TTG",'L'}, {"CTA",'L'}, {"CTC",'L'}, {"CTG",'L'}, {"CTT",'L'},
            {"TAC",'Y'}, {"TAT",'Y'},
            {"TAA",'*'}, {"TAG",'*'}, {"TGA",'*'},
            {"TGC",'C'}, {"TGT",'C'},
            {"TGG",'W'},
            {"CCA",'P'}, {"CCC",'P'}, {"CCG",'P'}, {"CCT",'P'},
            {"CAC",'H'}, {"CAT",'H'},
            {"CAA",'Q'}, {"CAG",'Q'},
            {"CGA",'R'}, {"CGC",'R'}, {"CGG",'R'}, {"CGT",'R'}, {"AGA",'R'}, {"AGG",'R'},
            {"ATA",'I'}, {"ATC",'I'}, {"ATT",'I'},
            {"ATG",'M'},
            {"ACA",'T'}, {"ACC",'T'}, {"ACG",'T'}, {"ACT",'T'},
            {"AAC",'N'}, {"AAT",'N'},
            {"AAA",'K'}, {"AAG",'K'},
            {"GTA",'V'}, {"GTC",'V'}, {"GTG",'V'}, {"GTT",'V'},
            {"GCA",'A'}, {"GCC",'A'}, {"GCG",'A'}, {"GCT",'A'},
            {"GAC",'D'}, {"GAT",'D'},
            {"GAA",'E'}, {"GAG",'E'},
            {"GGA",'G'}, {"GGC",'G'}, {"GGG",'G'}, {"GGT",'G'},
    };

    std::map<aminoAcidPair, int> blosum62 = {

            {aminoAcidPair('N','W'),-4}, {aminoAcidPair('D','L'),-4}, {aminoAcidPair('D','W'),-4}, {aminoAcidPair('C','E'),-4},
            {aminoAcidPair('E','C'),-4}, {aminoAcidPair('G','I'),-4}, {aminoAcidPair('G','L'),-4}, {aminoAcidPair('I','G'),-4},
            {aminoAcidPair('L','D'),-4}, {aminoAcidPair('L','G'),-4}, {aminoAcidPair('F','P'),-4}, {aminoAcidPair('P','F'),-4},
            {aminoAcidPair('P','W'),-4}, {aminoAcidPair('W','N'),-4}, {aminoAcidPair('W','D'),-4}, {aminoAcidPair('W','P'),-4},

            {aminoAcidPair('A','W'),-3}, {aminoAcidPair('R','C'),-3}, {aminoAcidPair('R','I'),-3}, {aminoAcidPair('R','F'),-3},
            {aminoAcidPair('R','W'),-3}, {aminoAcidPair('R','V'),-3}, {aminoAcidPair('N','C'),-3}, {aminoAcidPair('N','I'),-3},
            {aminoAcidPair('N','L'),-3}, {aminoAcidPair('N','F'),-3}, {aminoAcidPair('N','V'),-3}, {aminoAcidPair('D','C'),-3},
            {aminoAcidPair('D','I'),-3}, {aminoAcidPair('D','M'),-3}, {aminoAcidPair('D','F'),-3}, {aminoAcidPair('D','Y'),-3},
            {aminoAcidPair('D','V'),-3}, {aminoAcidPair('C','R'),-3}, {aminoAcidPair('C','N'),-3}, {aminoAcidPair('C','D'),-3},
            {aminoAcidPair('C','Q'),-3}, {aminoAcidPair('C','G'),-3}, {aminoAcidPair('C','H'),-3}, {aminoAcidPair('C','K'),-3},
            {aminoAcidPair('C','P'),-3}, {aminoAcidPair('Q','C'),-3}, {aminoAcidPair('Q','I'),-3}, {aminoAcidPair('Q','F'),-3},
            {aminoAcidPair('E','I'),-3}, {aminoAcidPair('E','L'),-3}, {aminoAcidPair('E','F'),-3}, {aminoAcidPair('E','W'),-3},
            {aminoAcidPair('G','C'),-3}, {aminoAcidPair('G','M'),-3}, {aminoAcidPair('G','F'),-3}, {aminoAcidPair('G','Y'),-3},
            {aminoAcidPair('G','V'),-3}, {aminoAcidPair('H','C'),-3}, {aminoAcidPair('H','I'),-3}, {aminoAcidPair('H','L'),-3},
            {aminoAcidPair('H','V'),-3}, {aminoAcidPair('I','R'),-3}, {aminoAcidPair('I','N'),-3}, {aminoAcidPair('I','D'),-3},
            {aminoAcidPair('I','Q'),-3}, {aminoAcidPair('I','E'),-3}, {aminoAcidPair('I','H'),-3}, {aminoAcidPair('I','K'),-3},
            {aminoAcidPair('I','P'),-3}, {aminoAcidPair('I','W'),-3}, {aminoAcidPair('L','N'),-3}, {aminoAcidPair('L','E'),-3},
            {aminoAcidPair('L','H'),-3}, {aminoAcidPair('L','P'),-3}, {aminoAcidPair('K','C'),-3}, {aminoAcidPair('K','I'),-3},
            {aminoAcidPair('K','F'),-3}, {aminoAcidPair('K','W'),-3}, {aminoAcidPair('M','D'),-3}, {aminoAcidPair('M','G'),-3},
            {aminoAcidPair('F','R'),-3}, {aminoAcidPair('F','N'),-3}, {aminoAcidPair('F','D'),-3}, {aminoAcidPair('F','Q'),-3},
            {aminoAcidPair('F','E'),-3}, {aminoAcidPair('F','G'),-3}, {aminoAcidPair('F','K'),-3}, {aminoAcidPair('P','C'),-3},
            {aminoAcidPair('P','I'),-3}, {aminoAcidPair('P','L'),-3}, {aminoAcidPair('P','Y'),-3}, {aminoAcidPair('S','W'),-3},
            {aminoAcidPair('W','A'),-3}, {aminoAcidPair('W','R'),-3}, {aminoAcidPair('W','E'),-3}, {aminoAcidPair('W','I'),-3},
            {aminoAcidPair('W','K'),-3}, {aminoAcidPair('W','S'),-3}, {aminoAcidPair('W','V'),-3}, {aminoAcidPair('Y','D'),-3},
            {aminoAcidPair('Y','G'),-3}, {aminoAcidPair('Y','P'),-3}, {aminoAcidPair('V','R'),-3}, {aminoAcidPair('V','N'),-3},
            {aminoAcidPair('V','D'),-3}, {aminoAcidPair('V','G'),-3}, {aminoAcidPair('V','H'),-3}, {aminoAcidPair('V','W'),-3},

            {aminoAcidPair('A','N'),-2}, {aminoAcidPair('A','D'),-2}, {aminoAcidPair('A','H'),-2}, {aminoAcidPair('A','F'),-2},
            {aminoAcidPair('A','Y'),-2}, {aminoAcidPair('R','D'),-2}, {aminoAcidPair('R','G'),-2}, {aminoAcidPair('R','L'),-2},
            {aminoAcidPair('R','P'),-2}, {aminoAcidPair('R','Y'),-2}, {aminoAcidPair('N','A'),-2}, {aminoAcidPair('N','M'),-2},
            {aminoAcidPair('N','P'),-2}, {aminoAcidPair('N','Y'),-2}, {aminoAcidPair('D','A'),-2}, {aminoAcidPair('D','R'),-2},
            {aminoAcidPair('C','F'),-2}, {aminoAcidPair('C','W'),-2}, {aminoAcidPair('C','Y'),-2}, {aminoAcidPair('Q','G'),-2},
            {aminoAcidPair('Q','L'),-2}, {aminoAcidPair('Q','W'),-2}, {aminoAcidPair('Q','V'),-2}, {aminoAcidPair('E','G'),-2},
            {aminoAcidPair('E','M'),-2}, {aminoAcidPair('E','Y'),-2}, {aminoAcidPair('E','V'),-2}, {aminoAcidPair('G','R'),-2},
            {aminoAcidPair('G','Q'),-2}, {aminoAcidPair('G','E'),-2}, {aminoAcidPair('G','H'),-2}, {aminoAcidPair('G','K'),-2},
            {aminoAcidPair('G','P'),-2}, {aminoAcidPair('G','T'),-2}, {aminoAcidPair('G','W'),-2}, {aminoAcidPair('H','A'),-2},
            {aminoAcidPair('H','G'),-2}, {aminoAcidPair('H','M'),-2}, {aminoAcidPair('H','P'),-2}, {aminoAcidPair('H','T'),-2},
            {aminoAcidPair('H','W'),-2}, {aminoAcidPair('I','S'),-2}, {aminoAcidPair('L','R'),-2}, {aminoAcidPair('L','Q'),-2},
            {aminoAcidPair('L','K'),-2}, {aminoAcidPair('L','S'),-2}, {aminoAcidPair('L','W'),-2}, {aminoAcidPair('K','G'),-2},
            {aminoAcidPair('K','L'),-2}, {aminoAcidPair('K','Y'),-2}, {aminoAcidPair('K','V'),-2}, {aminoAcidPair('M','N'),-2},
            {aminoAcidPair('M','E'),-2}, {aminoAcidPair('M','H'),-2}, {aminoAcidPair('M','P'),-2}, {aminoAcidPair('F','A'),-2},
            {aminoAcidPair('F','C'),-2}, {aminoAcidPair('F','S'),-2}, {aminoAcidPair('F','T'),-2}, {aminoAcidPair('P','R'),-2},
            {aminoAcidPair('P','N'),-2}, {aminoAcidPair('P','G'),-2}, {aminoAcidPair('P','H'),-2}, {aminoAcidPair('P','M'),-2},
            {aminoAcidPair('P','V'),-2}, {aminoAcidPair('S','I'),-2}, {aminoAcidPair('S','L'),-2}, {aminoAcidPair('S','F'),-2},
            {aminoAcidPair('S','Y'),-2}, {aminoAcidPair('S','V'),-2}, {aminoAcidPair('T','G'),-2}, {aminoAcidPair('T','H'),-2},
            {aminoAcidPair('T','F'),-2}, {aminoAcidPair('T','W'),-2}, {aminoAcidPair('T','Y'),-2}, {aminoAcidPair('W','C'),-2},
            {aminoAcidPair('W','Q'),-2}, {aminoAcidPair('W','G'),-2}, {aminoAcidPair('W','H'),-2}, {aminoAcidPair('W','L'),-2},
            {aminoAcidPair('W','T'),-2}, {aminoAcidPair('Y','A'),-2}, {aminoAcidPair('Y','R'),-2}, {aminoAcidPair('Y','N'),-2},
            {aminoAcidPair('Y','C'),-2}, {aminoAcidPair('Y','E'),-2}, {aminoAcidPair('Y','K'),-2}, {aminoAcidPair('Y','S'),-2},
            {aminoAcidPair('Y','T'),-2}, {aminoAcidPair('V','Q'),-2}, {aminoAcidPair('V','E'),-2}, {aminoAcidPair('V','K'),-2},
            {aminoAcidPair('V','P'),-2}, {aminoAcidPair('V','S'),-2},

            {aminoAcidPair('A','R'),-1}, {aminoAcidPair('A','Q'),-1}, {aminoAcidPair('A','E'),-1}, {aminoAcidPair('A','I'),-1},
            {aminoAcidPair('A','L'),-1}, {aminoAcidPair('A','K'),-1}, {aminoAcidPair('A','M'),-1}, {aminoAcidPair('A','P'),-1},
            {aminoAcidPair('R','A'),-1}, {aminoAcidPair('R','M'),-1}, {aminoAcidPair('R','S'),-1}, {aminoAcidPair('R','T'),-1},
            {aminoAcidPair('D','G'),-1}, {aminoAcidPair('D','H'),-1}, {aminoAcidPair('D','K'),-1}, {aminoAcidPair('D','P'),-1},
            {aminoAcidPair('D','T'),-1}, {aminoAcidPair('C','I'),-1}, {aminoAcidPair('C','L'),-1}, {aminoAcidPair('C','M'),-1},
            {aminoAcidPair('C','S'),-1}, {aminoAcidPair('C','T'),-1}, {aminoAcidPair('C','V'),-1}, {aminoAcidPair('Q','A'),-1},
            {aminoAcidPair('Q','P'),-1}, {aminoAcidPair('Q','T'),-1}, {aminoAcidPair('Q','Y'),-1}, {aminoAcidPair('E','A'),-1},
            {aminoAcidPair('E','P'),-1}, {aminoAcidPair('E','T'),-1}, {aminoAcidPair('G','D'),-1}, {aminoAcidPair('H','D'),-1},
            {aminoAcidPair('H','K'),-1}, {aminoAcidPair('H','F'),-1}, {aminoAcidPair('H','S'),-1}, {aminoAcidPair('I','A'),-1},
            {aminoAcidPair('I','C'),-1}, {aminoAcidPair('I','T'),-1}, {aminoAcidPair('I','Y'),-1}, {aminoAcidPair('L','A'),-1},
            {aminoAcidPair('L','C'),-1}, {aminoAcidPair('L','T'),-1}, {aminoAcidPair('L','Y'),-1}, {aminoAcidPair('K','A'),-1},
            {aminoAcidPair('K','D'),-1}, {aminoAcidPair('K','H'),-1}, {aminoAcidPair('K','M'),-1}, {aminoAcidPair('K','P'),-1},
            {aminoAcidPair('K','T'),-1}, {aminoAcidPair('M','A'),-1}, {aminoAcidPair('M','R'),-1}, {aminoAcidPair('M','C'),-1},
            {aminoAcidPair('M','K'),-1}, {aminoAcidPair('M','S'),-1}, {aminoAcidPair('M','T'),-1}, {aminoAcidPair('M','W'),-1},
            {aminoAcidPair('M','Y'),-1}, {aminoAcidPair('F','H'),-1}, {aminoAcidPair('F','V'),-1}, {aminoAcidPair('P','A'),-1},
            {aminoAcidPair('P','D'),-1}, {aminoAcidPair('P','Q'),-1}, {aminoAcidPair('P','E'),-1}, {aminoAcidPair('P','K'),-1},
            {aminoAcidPair('P','S'),-1}, {aminoAcidPair('P','T'),-1}, {aminoAcidPair('S','R'),-1}, {aminoAcidPair('S','C'),-1},
            {aminoAcidPair('S','H'),-1}, {aminoAcidPair('S','M'),-1}, {aminoAcidPair('S','P'),-1}, {aminoAcidPair('T','R'),-1},
            {aminoAcidPair('T','D'),-1}, {aminoAcidPair('T','C'),-1}, {aminoAcidPair('T','Q'),-1}, {aminoAcidPair('T','E'),-1},
            {aminoAcidPair('T','I'),-1}, {aminoAcidPair('T','L'),-1}, {aminoAcidPair('T','K'),-1}, {aminoAcidPair('T','M'),-1},
            {aminoAcidPair('T','P'),-1}, {aminoAcidPair('W','M'),-1}, {aminoAcidPair('Y','Q'),-1}, {aminoAcidPair('Y','I'),-1},
            {aminoAcidPair('Y','L'),-1}, {aminoAcidPair('Y','M'),-1}, {aminoAcidPair('Y','V'),-1}, {aminoAcidPair('V','C'),-1},
            {aminoAcidPair('V','F'),-1}, {aminoAcidPair('V','Y'),-1},

            {aminoAcidPair('A','C'),0}, {aminoAcidPair('A','G'),0}, {aminoAcidPair('A','T'),0}, {aminoAcidPair('A','V'),0},
            {aminoAcidPair('R','N'),0}, {aminoAcidPair('R','E'),0}, {aminoAcidPair('R','H'),0}, {aminoAcidPair('N','R'),0},
            {aminoAcidPair('N','Q'),0}, {aminoAcidPair('N','E'),0}, {aminoAcidPair('N','G'),0}, {aminoAcidPair('N','K'),0},
            {aminoAcidPair('N','T'),0}, {aminoAcidPair('D','Q'),0}, {aminoAcidPair('D','S'),0}, {aminoAcidPair('C','A'),0},
            {aminoAcidPair('Q','N'),0}, {aminoAcidPair('Q','D'),0}, {aminoAcidPair('Q','H'),0}, {aminoAcidPair('Q','M'),0},
            {aminoAcidPair('Q','S'),0}, {aminoAcidPair('E','R'),0}, {aminoAcidPair('E','N'),0}, {aminoAcidPair('E','H'),0},
            {aminoAcidPair('E','S'),0}, {aminoAcidPair('G','A'),0}, {aminoAcidPair('G','N'),0}, {aminoAcidPair('G','S'),0},
            {aminoAcidPair('H','R'),0}, {aminoAcidPair('H','Q'),0}, {aminoAcidPair('H','E'),0}, {aminoAcidPair('I','F'),0},
            {aminoAcidPair('L','F'),0}, {aminoAcidPair('K','N'),0}, {aminoAcidPair('K','S'),0}, {aminoAcidPair('M','Q'),0},
            {aminoAcidPair('M','F'),0}, {aminoAcidPair('F','I'),0}, {aminoAcidPair('F','L'),0}, {aminoAcidPair('F','M'),0},
            {aminoAcidPair('S','D'),0}, {aminoAcidPair('S','Q'),0}, {aminoAcidPair('S','E'),0}, {aminoAcidPair('S','G'),0},
            {aminoAcidPair('S','K'),0}, {aminoAcidPair('T','A'),0}, {aminoAcidPair('T','N'),0}, {aminoAcidPair('T','V'),0},
            {aminoAcidPair('V','A'),0}, {aminoAcidPair('V','T'),0},

            {aminoAcidPair('A','S'),1}, {aminoAcidPair('R','Q'),1}, {aminoAcidPair('N','D'),1}, {aminoAcidPair('N','H'),1},
            {aminoAcidPair('N','S'),1}, {aminoAcidPair('D','N'),1}, {aminoAcidPair('Q','R'),1}, {aminoAcidPair('Q','K'),1},
            {aminoAcidPair('E','K'),1}, {aminoAcidPair('H','N'),1}, {aminoAcidPair('I','M'),1}, {aminoAcidPair('L','V'),1},
            {aminoAcidPair('K','Q'),1}, {aminoAcidPair('K','E'),1}, {aminoAcidPair('M','I'),1}, {aminoAcidPair('M','V'),1},
            {aminoAcidPair('F','W'),1}, {aminoAcidPair('S','A'),1}, {aminoAcidPair('S','N'),1}, {aminoAcidPair('S','T'),1},
            {aminoAcidPair('T','S'),1}, {aminoAcidPair('W','F'),1}, {aminoAcidPair('V','L'),1}, {aminoAcidPair('V','M'),1},

            {aminoAcidPair('R','K'),2}, {aminoAcidPair('D','E'),2}, {aminoAcidPair('Q','E'),2}, {aminoAcidPair('E','D'),2},
            {aminoAcidPair('E','Q'),2}, {aminoAcidPair('H','Y'),2}, {aminoAcidPair('I','L'),2}, {aminoAcidPair('L','I'),2},
            {aminoAcidPair('L','M'),2}, {aminoAcidPair('K','R'),2}, {aminoAcidPair('M','L'),2}, {aminoAcidPair('W','Y'),2},
            {aminoAcidPair('Y','H'),2}, {aminoAcidPair('Y','W'),2},

            {aminoAcidPair('I','V'),3}, {aminoAcidPair('F','Y'),3}, {aminoAcidPair('Y','F'),3}, {aminoAcidPair('V','I'),3},

            {aminoAcidPair('A','A'),4}, {aminoAcidPair('I','I'),4}, {aminoAcidPair('L','L'),4}, {aminoAcidPair('S','S'),4},
            {aminoAcidPair('V','V'),4},

            {aminoAcidPair('R','R'),5}, {aminoAcidPair('Q','Q'),5}, {aminoAcidPair('E','E'),5}, {aminoAcidPair('K','K'),5},
            {aminoAcidPair('M','M'),5}, {aminoAcidPair('T','T'),5},

            {aminoAcidPair('N','N'),6}, {aminoAcidPair('D','D'),6}, {aminoAcidPair('G','G'),6}, {aminoAcidPair('F','F'),6},

            {aminoAcidPair('P','P'),7}, {aminoAcidPair('Y','Y'),7},

            {aminoAcidPair('H','H'),8},

            {aminoAcidPair('C','C'),9},

            {aminoAcidPair('W','W'),11},
    };

    char translateCodon(std::string codon){
        if (toupper(codon[0])=='N' || toupper(codon[1])=='N' || toupper(codon[2])=='N'){
            return ' ';
        } else {
            return codonTable.at(codon);
        }
    }

    int getBlosum62Score(aminoAcidPair aaPair){
        if (aaPair.first==' ' || aaPair.second==' '){
            return -4; // minimum blosum62 score
        } else {
            return blosum62.at(aaPair);
        }
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

    int getNumMatches(std::string backtrace){
        int startPos = 0;
        int digits = 0;
        int numMatches = 0;
        for (size_t pos = 0; pos < backtrace.size(); pos++) {
            if(!isdigit(backtrace[pos])){
                int num = std::stoi(backtrace.substr(startPos, digits));
                char cha = backtrace[pos];
                if (cha == 'M')
                    numMatches += num;
                startPos = pos + 1;
                digits=0;
            } else{
                digits++;
            }
        }
        return numMatches;
    }

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

    char * getQuerySequence(unsigned int queryKey, int thread_idx){
        unsigned int queryId = qDbr->sequenceReader->getId(queryKey);
        char * qSeq = qDbr->sequenceReader->getData(queryId, thread_idx);
        return qSeq;
    }
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
//    bo
//    ol metContainF(char * qSeq, int qStartPos, int qEndPos){
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
// TODO
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
    // DO I NEED?
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
//            if (queryKey==498)
//                std::cout<<""<<std::endl;
            ExonFinder exonFinder(&tDbr, &qDbr, queryKey);
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
                            thread_idx,
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
                        thread_idx,
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



