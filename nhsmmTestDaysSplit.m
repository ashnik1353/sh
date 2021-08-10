function outTesting = nhsmmTestDaysSplit(curExp, learnedParams)

%Hidden semi-Markov model, aka variable duration Markov model

[inferedLabels, fwProbs, bestDur, inferedSegments, obsModel] = nhsmmInferenceDaysSplit(curExp.modelInfo, learnedParams, curExp.testFeatMat);

%[dummy,bestDur] = hsmmFwBkDaysSplit(modelParams, learnedParams, testFeatMat);
% [Y,I] = max(bestDur);
% inferedLabels = I;

outTesting.inferedLabels = inferedLabels;
outTesting.inferedSegments = inferedSegments;
outTesting.fwProbs = fwProbs;
outTesting.bestDur = bestDur;
outTesting.obsModel = obsModel;
