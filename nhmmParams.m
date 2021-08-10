function outTrain = nhmmParams(curExp, modelInfo)

smallValue = modelInfo.smallValue;
nState = modelInfo.nState;


trainFeatMat = curExp.points;
trainIndex = curExp.ind;

% Used to avoid zero probabilties

% [~,idx] = sort(std(trainFeatMat'),'descend');
[T,idx] = sort(sum(trainFeatMat,2),'descend');
T = normalise(T,1);
j=0;
s=0;
if nState ~=1
    while s<nState && j<size(trainFeatMat,1)
        j = j+1;
        s = s + T(j);
    end
end
if nState == 1;
    states=1:size(trainFeatMat,1);
else
     states = idx(1:j);
end
numStates = length(states)+1;

% Learning of HMM parameters through maximum likelihood. 
% By counting parameters will be estimated
% T = nfstate(trainFeatMat,states);

numTrain = size(trainIndex,2)-1;

% Parameters
prior = zeros(numStates,1);
obsModel = zeros(numStates,1);
transModel = zeros(numStates, numStates);

for n=1:numTrain,
    % initial state distribution: sum[Label(n,1)]/totalNumberDays
    trainPoint = trainFeatMat(:,trainIndex(n):trainIndex(n+1)-1);
%     I = T(n);
    I = nfstate(trainPoint(:,1),states);
    prior(I) = prior(I)+1;
    obsModel(I) = obsModel(I) + 1;    
    t=I;
    for i=2:size(trainPoint,2)
%         s = T(n+i-1);
        s = nfstate(trainPoint(:,i), states);
        obsModel(s)= obsModel(s) +1;
        transModel(t,s) = transModel(t,s) + 1;        
        t = s;
    end
end

% La place smoothing
prior = (prior+smallValue)/(numTrain+smallValue);
% totSense = sum(trainFeatMat,2);
obsModel=(obsModel+smallValue)/(size(trainFeatMat,2)+smallValue);
transModel = normalise(transModel+smallValue,2);

outTrain.states = states;
outTrain.prior = prior;
outTrain.obsModel = obsModel;
outTrain.transModel = transModel;




   
    
    

