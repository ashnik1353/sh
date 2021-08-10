function [obsModel] = npreComputeProbs(modelInfo, testFeatMat, S, maxDur ,sigmas , LP )
	numTimeSteps = size(testFeatMat,2);
	numAct = modelInfo.numAct;
	thau = modelInfo.thau;
    size1 = modelInfo.size1;
    distType = modelInfo.distType;
    relative = modelInfo.relative;
    subsWeight = modelInfo.subsWeight;
    Itakura = modelInfo.Itakura;
    
    
    numSense = modelInfo.numSense;
    
    
	%sigma = modelInfo.sigma;
	obsModel = zeros(numAct, numTimeSteps, maxDur);
    parfor i=1:numAct
        parObsModel = zeros(numTimeSteps, maxDur);
        issize1 = ismember(i, size1);
        Sindex = S{i}.ind;
        Spoints = S{i}.points;
         %%%%%%%%%%%%%%ezf%%%%%%%%%%%%%%%%%%%%%%%
        %pcas = LP{i}.pca;
        %testFeatMatpca = pcas(:,1:4)'* testFeatMat;
        
        
            testFeatMatpca =  testFeatMat;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if length(Sindex)==1        
            parObsModel(:,:)=0.01;
            continue;
        end 
		numLernedPoints = size(Sindex,2)-1;
        for k=1:numTimeSteps-maxDur+1
            distance = zeros(numLernedPoints, maxDur);
			
            testPoint = testFeatMatpca(:,k:k+maxDur-1);
            
            for j=1:numLernedPoints,					
				modelPoint = Spoints(:,Sindex(j):Sindex(j+1)-1);
                %Hammings = bsxfun(@plus,sum(modelPoint,1)',sum(testFeatMatpca,1))-2*(modelPoint'*testFeatMatpca);
                switch distType 
                    case 0
                        distance(j,:) = nLevenstein(testPoint, modelPoint, thau, relative , subsWeight);
                    case 1
                        distance(j,:) = nDTW(testPoint, modelPoint, thau, Itakura);
                    case 2
                        distance(j,:) = nLevenstein2(testPoint, modelPoint, thau, relative,subsWeight);
                end
                
            end
            
			probs = ncompProbs(distance, sigmas{i}, issize1, modelInfo);
          
            
			for l = 1:maxDur
				parObsModel(k+l-1, l) = probs(l);
			end
        end
        for m=k+1:numTimeSteps
            distance = zeros(numLernedPoints, numTimeSteps-m+1);
			
            testPoint = testFeatMatpca(:,m:numTimeSteps);
            
            for j=1:numLernedPoints					
				modelPoint = Spoints(:,Sindex(j):Sindex(j+1)-1);
                switch distType 
                    case 0
                        distance(j,:) = nLevenstein(testPoint, modelPoint, thau, relative, subsWeight);
                    case 1
                        distance(j,:) = nDTW(testPoint, modelPoint, thau,Itakura);
                    case 2
                        distance(j,:) = nLevenstein2(testPoint, modelPoint, thau, relative,subsWeight);
                end
                
            end
            
			probs = ncompProbs(distance, sigmas{i}, issize1, modelInfo);
            
            
			for l = 1:numTimeSteps-m+1
				parObsModel(m+l-1, l) = probs(l);
			end
        end
        obsModel(i,:,:)=parObsModel;
    end
    
end
