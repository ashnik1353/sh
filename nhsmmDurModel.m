function mvParms = nhsmmDurModel(durCell, numStates, nBin)

stateList = 1:numStates;
%use maxdur from modelInfo



    
    numParms = nBin;

mvParms.values= zeros(length(stateList), numParms);

for i=1:length(durCell),
      
        mvParms.maxDur4State = max(max(durCell{i}), nBin);
        mvParms.binSize(i) = max(1,mvParms.maxDur4State/nBin);
        mvParms.values(i,:)=histc(durCell{i},1:mvParms.binSize(i):mvParms.maxDur4State);
    
end

mvParms.values=mvParms.values + 0.01;

mvParms.values = mvParms.values./repmat(sum(mvParms.values,2),1,numParms);