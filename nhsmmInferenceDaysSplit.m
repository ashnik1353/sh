function [inferedLabelsDays, bestProbDays, bestDurDays, inferedSegmentsDays, obsModelDays] = nhsmmInferenceDaysSplit(modelInfo, learnedParams, testFeatMatDays ,varargin)
initial = learnedParams.initial;
obsModelHSMM = learnedParams.obsModelHSMM;
transModel = learnedParams.transModel;



numAct = modelInfo.numAct;
numVals = modelInfo.numVals;
%%%%%%%%%%%ezf%%%%%%%%%%%
useCluster = modelInfo.useCluster;

obsModelDays = cell(length(testFeatMatDays),1);
% obsModelDays1 = cell(length(testFeatMatDays)*8,1);
% h = length(testFeatMatDays);
%%%%%%%%%%%%%%%%%%%%%%%%%
inferedLabelsDays = cell(length(testFeatMatDays),1);
bestProbDays = cell(length(testFeatMatDays),1);
bestDurDays = cell(length(testFeatMatDays),1);
inferedSegmentsDays = cell(length(testFeatMatDays),1);

% inferedLabelsDays1 = cell(length(testFeatMatDays)*8,1);
% bestProbDays1 = cell(length(testFeatMatDays)*8,1);
% bestDurDays1 = cell(length(testFeatMatDays)*8,1);
% inferedSegmentsDays1 = cell(length(testFeatMatDays)*8,1);


%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin >= 4)
    testFeatMatDays1 = varargin{1};
else
    testFeatMatDays1 = testFeatMatDays;
end
% for i=1:length(testFeatMatDays1)
%     s= size(testFeatMatDays1{i},2);
%     for j=1:8
%         
%         testFeatMatDays{j+8*(i-1)} = testFeatMatDays1{i}(:,ceil(s*(j-1)/8)+1:ceil(s*j/8));
%     end
% end
% testFeatMatDays1 = testFeatMatDays;

%%%%%%%%%%%%%%%%%%%%%%%%%%%5


for n=1:length(testFeatMatDays),    
    testFeatMat = testFeatMatDays{n};
    numTimeSteps = size(testFeatMat,2); 

    maxDur = learnedParams.maxDur;
	%%%%%%%%%%%%%%%%%%%ezf%%%%%%%%%%%%%%%%%%
	
	testFeatMat1 = testFeatMatDays1{n};
    if useCluster 
        
            obsModel = npreComputeProbs(modelInfo, testFeatMat1, learnedParams.medoids, maxDur,learnedParams.sigmas, learnedParams.S);
        
    else
        obsModel = npreComputeProbs1(modelInfo, testFeatMat1, learnedParams.S, maxDur);
        %obsModel = npreComputeProbs2(modelInfo, testFeatMat, learnedParams.S, maxDur, learnedParams.sigmas);
    end
    %obsModel = normalise(obsModel,1);
	obsModel = log(obsModel);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Viterbi algorithm as described in Rabiner 1989 (page 264)
    inferedSegments = zeros(1,numTimeSteps);
    inferedLabels = zeros(1,numTimeSteps);
    bestProb = zeros(numAct,numTimeSteps);
    bestState = zeros(numAct,numTimeSteps);
    bestDur = zeros(numAct,numTimeSteps);

    probObs = ones(numAct,numTimeSteps);
    probLogDurs= zeros(numAct,maxDur);

    % Precompute duration probabilities
    for d=1:maxDur,
        probLogDurs(:,d) = hsmmDurationModel(d,modelInfo, learnedParams);
    end
    probLogDurs = log(normalise(probLogDurs,2));

%     h = waitbar(0,'Please wait...');
    for k=1:numTimeSteps,
%         waitbar(k/numTimeSteps,h);
        %%% 2) Recursion
        
        % Calculate observation probability
        tempProbObs = ones(1,numAct);
        for i=1:numVals,
            idxVal = find(testFeatMat(:,k)==modelInfo.obsList(i));
            probCurObs = prod(obsModelHSMM(idxVal,i,:),1);
            tempProbObs=tempProbObs .* reshape(probCurObs,1,numAct);
        end
        probObs(:,k) = log(tempProbObs)';

        tempCumObsProb = 0;
        
        
        probAlphaNoTrans = zeros(numAct,min(k,maxDur));
        bestState4Dur = zeros(numAct,min(k,maxDur));

        D = min(k,maxDur);
        for d=1:D,
            if (k==d)        %%% 1) Initialization
                tempCumObsProb = tempCumObsProb + probObs(:, 1);
                
                    
                        probAlphaNoTrans(:,d) = log(initial) + obsModel(:,k,d) + tempCumObsProb + probLogDurs(:,d);                
                    
                
                
            else  %% Transition             
                % Determine best transition
                repBestProbs = repmat(bestProb(:,k-d),1,numAct);
                [bestProb4Dur(:,1),bestState4Dur(:,d)]= max(repBestProbs + log(transModel),[],1);

                % Add observation and duration
                tempCumObsProb = tempCumObsProb + probObs(:, k-d+1);
                
                   
                        probAlphaNoTrans(:,d) = bestProb4Dur + obsModel(:,k,d) + tempCumObsProb + probLogDurs(:,d);
                    
                
                
            end
        end
%         assert(~any(isnan(probAlphaNoTrans)));

        % Determine best duration for each state
        [bestProb(:,k), bestDur(:,k)] = max(probAlphaNoTrans, [], 2);
        % Store best state to get from. NOTE! if zero it means it's the initial state
        for i=1:numAct,
            bestState(i,k)=bestState4Dur(i,bestDur(i,k));
        end
    end
%     close(h);

    %%% 3) Termination
    [P,lastState] = max(bestProb(:,numTimeSteps));
    inferedLabels(numTimeSteps) = lastState;
    curSegment = -1;
    inferedSegments(numTimeSteps)= curSegment;
    duration = bestDur(inferedLabels(numTimeSteps),numTimeSteps);
    tempDur = bestDur(inferedLabels(numTimeSteps),numTimeSteps);
    tempDur = tempDur -1;
    %%% 4) Path backtracking

    for k=numTimeSteps-1:-1:1
        if (tempDur > 0)
            inferedLabels(k) = lastState;
            inferedSegments(k)=curSegment;
            tempDur = tempDur -1;
        else
            lastState = bestState(inferedLabels(k+duration),k+duration);
            inferedLabels(k) = lastState;
            curSegment = curSegment - 1;
            inferedSegments(k)= curSegment;

            duration = bestDur(inferedLabels(k),k);
            tempDur = bestDur(inferedLabels(k),k)-1;
        end
    end

    inferedSegments = inferedSegments+abs(inferedSegments(1))+1;

    inferedLabelsDays{n} = inferedLabels;
    bestProbDays{n} = bestProb;
    bestDurDays{n} = bestDur;
    inferedSegmentsDays{n} = inferedSegments;
    obsModelDays{n} = obsModel;
end

    


    