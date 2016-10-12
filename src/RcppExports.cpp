// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// computeWinCount
List computeWinCount(IntegerVector startPos, IntegerVector endPos, IntegerVector startNeg, IntegerVector endNeg, int end, int win, int step, double logitThreshold, int minR);
RcppExport SEXP rnaCleanR_computeWinCount(SEXP startPosSEXP, SEXP endPosSEXP, SEXP startNegSEXP, SEXP endNegSEXP, SEXP endSEXP, SEXP winSEXP, SEXP stepSEXP, SEXP logitThresholdSEXP, SEXP minRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type startPos(startPosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type endPos(endPosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type startNeg(startNegSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type endNeg(endNegSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type step(stepSEXP);
    Rcpp::traits::input_parameter< double >::type logitThreshold(logitThresholdSEXP);
    Rcpp::traits::input_parameter< int >::type minR(minRSEXP);
    rcpp_result_gen = Rcpp::wrap(computeWinCount(startPos, endPos, startNeg, endNeg, end, win, step, logitThreshold, minR));
    return rcpp_result_gen;
END_RCPP
}
// computeWinCov
List computeWinCov(IntegerVector covPosLen, IntegerVector covPosVal, IntegerVector covNegLen, IntegerVector covNegVal, int end, int win, int step, double logitThreshold, int minR);
RcppExport SEXP rnaCleanR_computeWinCov(SEXP covPosLenSEXP, SEXP covPosValSEXP, SEXP covNegLenSEXP, SEXP covNegValSEXP, SEXP endSEXP, SEXP winSEXP, SEXP stepSEXP, SEXP logitThresholdSEXP, SEXP minRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type covPosLen(covPosLenSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type covPosVal(covPosValSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type covNegLen(covNegLenSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type covNegVal(covNegValSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type step(stepSEXP);
    Rcpp::traits::input_parameter< double >::type logitThreshold(logitThresholdSEXP);
    Rcpp::traits::input_parameter< int >::type minR(minRSEXP);
    rcpp_result_gen = Rcpp::wrap(computeWinCov(covPosLen, covPosVal, covNegLen, covNegVal, end, win, step, logitThreshold, minR));
    return rcpp_result_gen;
END_RCPP
}
// getIndex
List getIndex(CharacterVector strand);
RcppExport SEXP rnaCleanR_getIndex(SEXP strandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type strand(strandSEXP);
    rcpp_result_gen = Rcpp::wrap(getIndex(strand));
    return rcpp_result_gen;
END_RCPP
}
// keepReadCov
IntegerVector keepReadCov(IntegerVector startAl, IntegerVector endAl, CharacterVector strand, IntegerVector keepWinPos, IntegerVector keepWinNeg, int end, int win, int step);
RcppExport SEXP rnaCleanR_keepReadCov(SEXP startAlSEXP, SEXP endAlSEXP, SEXP strandSEXP, SEXP keepWinPosSEXP, SEXP keepWinNegSEXP, SEXP endSEXP, SEXP winSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type startAl(startAlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type endAl(endAlSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type strand(strandSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type keepWinPos(keepWinPosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type keepWinNeg(keepWinNegSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(keepReadCov(startAl, endAl, strand, keepWinPos, keepWinNeg, end, win, step));
    return rcpp_result_gen;
END_RCPP
}
// keepRead
List keepRead(int nbSample, IntegerVector startPosAl, IntegerVector endPosAl, IntegerVector groupNamePos, IntegerVector samplePos, IntegerVector startNegAl, IntegerVector endNegAl, IntegerVector groupNameNeg, IntegerVector sampleNeg, IntegerVector keepWinPos, IntegerVector keepWinNeg, int end, int win, int step);
RcppExport SEXP rnaCleanR_keepRead(SEXP nbSampleSEXP, SEXP startPosAlSEXP, SEXP endPosAlSEXP, SEXP groupNamePosSEXP, SEXP samplePosSEXP, SEXP startNegAlSEXP, SEXP endNegAlSEXP, SEXP groupNameNegSEXP, SEXP sampleNegSEXP, SEXP keepWinPosSEXP, SEXP keepWinNegSEXP, SEXP endSEXP, SEXP winSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nbSample(nbSampleSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type startPosAl(startPosAlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type endPosAl(endPosAlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groupNamePos(groupNamePosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type samplePos(samplePosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type startNegAl(startNegAlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type endNegAl(endNegAlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groupNameNeg(groupNameNegSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sampleNeg(sampleNegSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type keepWinPos(keepWinPosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type keepWinNeg(keepWinNegSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(keepRead(nbSample, startPosAl, endPosAl, groupNamePos, samplePos, startNegAl, endNegAl, groupNameNeg, sampleNeg, keepWinPos, keepWinNeg, end, win, step));
    return rcpp_result_gen;
END_RCPP
}
// keepReadOne
List keepReadOne(IntegerVector startPosAl, IntegerVector endPosAl, IntegerVector groupNamePos, IntegerVector startNegAl, IntegerVector endNegAl, IntegerVector groupNameNeg, IntegerVector keepWinPos, IntegerVector keepWinNeg, int end, int win, int step);
RcppExport SEXP rnaCleanR_keepReadOne(SEXP startPosAlSEXP, SEXP endPosAlSEXP, SEXP groupNamePosSEXP, SEXP startNegAlSEXP, SEXP endNegAlSEXP, SEXP groupNameNegSEXP, SEXP keepWinPosSEXP, SEXP keepWinNegSEXP, SEXP endSEXP, SEXP winSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type startPosAl(startPosAlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type endPosAl(endPosAlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groupNamePos(groupNamePosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type startNegAl(startNegAlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type endNegAl(endNegAlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groupNameNeg(groupNameNegSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type keepWinPos(keepWinPosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type keepWinNeg(keepWinNegSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(keepReadOne(startPosAl, endPosAl, groupNamePos, startNegAl, endNegAl, groupNameNeg, keepWinPos, keepWinNeg, end, win, step));
    return rcpp_result_gen;
END_RCPP
}
