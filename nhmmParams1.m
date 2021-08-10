function outTrain = nhmmParams1(curExp, modelInfo)

smallValue = modelInfo.smallValue;
nState = modelInfo.nState;

numAct = modelInfo.numAct;
numSense = modelInfo.numSense;
numVals = modelInfo.numVals;


T.trainFeatMat = curExp.points;
T.trainIndex = curExp.ind;
% row = size(T.trainIndex,2)-1;
% if  row>=500
%         b = ceil(0.05*row);
%         e = row-b;
%         T.trainFeatMat = [T.trainFeatMat(:,1:b), T.trainFeatMat(:,e+1:row)];
%         T.trainIndex = 1:2*b+1;
%         curExp.points = T.trainFeatMat;
%         curExp.ind = T.trainIndex;
% end

% Used to avoid zero probabilties

% [~,idx] = sort(std(curExp.points'),'descend');
[I,idx] = sort(sum(curExp.points,2),'descend');
I = normalise(I,1);
j=0;
s=0;
if nState ~=1
    while s<nState && j<size(curExp.points,1)
        j = j+1;
        s = s + I(j);
    end
end
% if nState == 1;
%     states=1:size(trainFeatMat,1);
% else
%      states = idx(1:j);
% end
% numStates = length(states)+1;

% Learning of HMM parameters through maximum likelihood. 
% By counting parameters will be estimated
% T = nfstate(trainFeatMat,states);

numStates = 4;
% numStates = j;

numTrain = size(T.trainIndex,2)-1;

% initial Parameters
L.numStates = numStates;
L.prior = normalise(rand(numStates,1),1);%/numState;
L.obsModel = rand(numSense,numStates);
L.transModel = normalise(rand(numStates, numStates),2);%/numState;

logLikelihood = 0;
tmp = 2;
j=0;
while abs(tmp-logLikelihood)>0.1 && j<15
    
    j=j+1;
     tmp = logLikelihood;
    [logLikelihood,fwbkProbsDays,fwbkNormedProbsDays,alpha,betha] = nhmmFwBkDaysSplit(modelInfo, L, curExp);
   
    numTrain = length(fwbkProbsDays);
    a = zeros(numStates,1);
    b = 0;
    for n=1:numTrain,
        a = a + fwbkNormedProbsDays{n}(:,1);
        b = b + sum(fwbkNormedProbsDays{n}(:,1),1);
    end
    L.prior = (a+0.01)/(b+0.01*numStates);
    L.prior = normalise(L.prior,1);
    ksi = zeros(numStates,numStates);
    b=zeros(1,numStates);
    for n=1:numTrain,
        trainPoint = T.trainFeatMat(:,T.trainIndex(n):T.trainIndex(n+1)-1);
        s = zeros(numStates,numStates);
        m = zeros(1,numStates);
        for i=2:size(trainPoint,2)
            point = trainPoint(:,i);
%             aa = ntr(L.transModel,point,L.obsModel,alpha{n}(:,i-1),betha{n}(:,i));
            s = s + ntr(L.transModel,point,L.obsModel,alpha{n}(:,i-1),betha{n}(:,i));
            m= m + sum (ksi,2)';            
        end 
        ksi = ksi + s;
        b = b + m;
    end
    L.transModel = (ksi+0.01) ./(repmat(b,numStates,1)+0.01*numStates);
    L.transModel = normalise(L.transModel,2);
    
    s=zeros(numSense,numStates);
    m=zeros(1,numStates);
    for n=1:numTrain,
        trainPoint = T.trainFeatMat(:,T.trainIndex(n):T.trainIndex(n+1)-1);
        s1=zeros(numSense,numStates);
        m1=zeros(1,numStates);
        for i=1:size(trainPoint,2)            
            s1 = s1 + repmat(fwbkNormedProbsDays{n}(:,i)',numSense,1).*repmat(trainPoint(:,i),1,numStates);
            m1 = m1 + fwbkNormedProbsDays{n}(:,i)';           
        end
        s = s + s1;
        m = m + m1;
    end
    L.obsModel = (s+0.01)./(repmat(m,numSense,1)+0.01*numStates);
    
end
% [j,v]


% La place smoothing

outTrain.numStates = numStates;
outTrain.prior = L.prior;
outTrain.obsModel = L.obsModel;
outTrain.transModel = L.transModel;




   
    
    

