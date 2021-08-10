function [obsModel] = npreComputeProbs1(modelInfo, testFeatMat, S, maxDur)
	numTimeSteps = size(testFeatMat,2);
	numAct = modelInfo.numAct;
	thau = modelInfo.thau;
	sigma = modelInfo.sigma;
	obsModel = zeros(numAct, numTimeSteps, maxDur);
    for i=1:numAct
        Sindex = S{i}.ind;
        Spoints = S{i}.points;
        %%%%%%%%%%%%%%ezf%%%%%%%%%%%%%%%%%%%%%%%
        Sprobs = S{i}.probs;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if length(Sindex)==1        
            obsModel(i,:,:)=0.00001;
            continue;
        end 
		numLernedPoints = size(Sindex,2)-1;
        for k=1:numTimeSteps-maxDur+1
			distance = zeros(numLernedPoints, maxDur);
            testPoint = testFeatMat(:,k:k+maxDur-1);
			for j=1:numLernedPoints					
				modelPoint = Spoints(:,Sindex(j):Sindex(j+1)-1);
				distance(j,:) = nLevenstein(testPoint, modelPoint, thau);
                %distance(j,:) = 1:maxDur;
			end
			probs = ncompProbs1(distance, sigma, Sprobs);
			for l = 1:maxDur
				obsModel(i,k+l-1, l) = probs(l);
			end
        end
        for m=k+1:numTimeSteps
			distance = zeros(numLernedPoints, numTimeSteps-m+1);
            testPoint = testFeatMat(:,m:numTimeSteps);
			for j=1:numLernedPoints					
				modelPoint = Spoints(:,Sindex(j):Sindex(j+1)-1);
				distance(j,:) = nLevenstein(testPoint, modelPoint, thau);	
                %distance(j,:) = 1:numTimeSteps-m;
			end
			probs = ncompProbs1(distance, sigma, Sprobs);
			for l = 1:numTimeSteps-m+1
				obsModel(i,m+l-1, l) = probs(l);
			end
        end
    end
end
