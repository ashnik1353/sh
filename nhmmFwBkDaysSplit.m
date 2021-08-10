function [logLikelihood, fwbkProbsDays,fwbkNormedProbsDays,alpha,betha] = nhmmFwBkDaysSplit(modelInfo, learnedParams, T)

dummy='dummy';
prior = learnedParams.prior;
obsModel = learnedParams.obsModel;
transModel = learnedParams.transModel;

% numStates = modelInfo.numStates;
numStates = 4;
% numStates = learnedParams.numStates;
numSense = modelInfo.numSense;
numVals = modelInfo.numVals;

trainFeatMat = T.points;
trainIndex = T.ind;
numTrain = size(trainIndex,2)-1;

fwbkProbsDays = cell(numTrain,1);
fwbkNormedProbsDays = cell(numTrain,1);
alpha = cell(numTrain,1);
betha = cell(numTrain,1);
px = cell(numTrain,1);
logLikelihood = log(1);

for n=1:numTrain,
    testFeatMat = trainFeatMat(:,trainIndex(n):trainIndex(n+1)-1);
    numTimeSteps = size(testFeatMat,2); 

    % Forward Backward Procedure as  as described in Rabiner 1989 (page 262)

    fwbkProbs = zeros(numStates, numTimeSteps);
    forwVar = zeros(numStates, numTimeSteps);
    backVar = zeros(numStates, numTimeSteps);

    %%% Forward
    % alpha_t(i) = p(o_1, o_2, .. o_t, q_t= s_i|model_parameters)
    % o = observation, q = state

    % Calculate observation probability
    
    
    P= repmat(testFeatMat(:,1),1,numStates);
    tempProbObs = prod(P.*obsModel+(1-P).*(1-obsModel),1);
    probObs = log(tempProbObs);

    %%% 1) Initialization
    forwVar(:,1) = log(prior) + probObs';

    for k=2:numTimeSteps,

        %%% 2) Induction

        % Calculate bestProb * transModel
        repBestProbs = repmat(forwVar(:,k-1),1,numStates);
        maxLogTransProb = logsum(repBestProbs + log(transModel),1);

        % Calculate observation probability
        
        P= repmat(testFeatMat(:,k),1,numStates);
        tempProbObs = prod(P.*obsModel+(1-P).*(1-obsModel),1);
        probObs = log(tempProbObs);

        forwVar(:,k) = maxLogTransProb' + probObs';
    end

    %%% 3) Termination
    %probOfObservations = sum(forwVar(:,numTimeSteps)); % Normalized, doesn't
                                                        % make sense

    %%% Backward
    % Beta_t(i) = P(o_{t+1}, o_{t+2}, .., o_T | q_t=s_i)
    %  o = observation, q = state

    %%% 1) Initialization
    backVar(:,numTimeSteps)=log(1);

    for k=(numTimeSteps-1):-1:1,
        %%% 2) Induction

        % Calculate bestProb * transModel
        repBestProbs = repmat(backVar(:,k+1),1,numStates);
        maxLogTransProb = logsum(repBestProbs + log(transModel),1);

        % Calculate observation probability
        P= repmat(testFeatMat(:,k),1,numStates);
        tempProbObs = prod(P.*obsModel+(1-P).*(1-obsModel),1);
        probObs = log(tempProbObs);

        backVar(:,k) = maxLogTransProb' + probObs';
    end

    %%% Forward-Backward combined
    % gamma_t(i) = P(q_t=S_i|O, model_parameters)
    % O = all observations (t=1:T) q = state

    for k=1:numTimeSteps,
        fwbkProbs(:,k) = forwVar(:,k)+backVar(:,k);

        % Normalize 
        fwbkNormedProbs(:,k) = normalise(exp(fwbkProbs(:,k)),1);
    end
    px{n} = log(sum(exp(forwVar(:,numTimeSteps))));
    alpha{n} = forwVar;
    betha{n} = backVar;
    logLikelihood = logLikelihood + px{n};
    fwbkProbsDays{n} = fwbkProbs;
    fwbkNormedProbsDays{n} = fwbkNormedProbs;
end