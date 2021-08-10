function outTrain = nhsmmParams(curExp, modelInfo)

nBin = modelInfo.nBin;
nState = modelInfo.nState;


trainFeatMat = curExp.points;
trainIndex = curExp.ind;

% Used to avoid zero probabilties

% [~,idx] = sort(std(trainFeatMat'),'descend');
[~,idx] = sort(sum(trainFeatMat,2),'descend');
smallValue = 0.01;
if nState == 1;
    states=1:size(trainFeatMat,1);
else
     states = idx(1:round(nState * size(trainFeatMat,1)));
end
numStates = length(states)+1;

% Learning of HMM parameters through maximum likelihood. 
% By counting parameters will be estimated

numTrain = size(trainIndex,2)-1;

% Parameters
prior = zeros(numStates,1);
obsModel = zeros(numStates,1);
transModel = zeros(numStates, numStates);




durs = cell(numStates,1);
for s=1:numStates
        durs{s} = [];
end
for n=1:numTrain,
    % initial state distribution: sum[Label(n,1)]/totalNumberDays
    trainPoint = trainFeatMat(:,trainIndex(n):trainIndex(n+1)-1);
    I = nfstate(trainPoint(:,1),states);
    prior(I) = prior(I)+1;
    obsModel(I) = obsModel(I) + 1;    
    t=I;
    j=1;
    if size(trainPoint,2)==1
        durs{s} = [durs{s},1];
    end
    for i=2:size(trainPoint,2)
        s = nfstate(trainPoint(:,i), states);
        
        obsModel(s)= obsModel(s) +1;
        if (s~=t)
            transModel(t,s) = transModel(t,s) + 1;
            t = s;
            durs{s} = [durs{s},j];
            j=1;
        else
            j=j+1;
        end
        
    end
    if j~=1
       durs{s} = [durs{s},j];
    end
    
    
end
for s=1:numStates
    if isempty(durs{s})
        durs{s}=0;
    end
end    
%     if ~isempty(durs{s})
        
        
        
        durModel=nhsmmDurModel(durs, numStates, nBin);
        
%     end


% La place smoothing
prior = (prior+smallValue)/(numTrain+smallValue);
% totSense = sum(trainFeatMat,2);
obsModel=(obsModel+smallValue)/(size(trainFeatMat,2)+smallValue);
transModel = normalise(transModel+smallValue,2);

outTrain.states = states;
outTrain.prior = prior;
outTrain.obsModel = obsModel;
outTrain.transModel = transModel;
outTrain.durModel = durModel;




   
    
    

