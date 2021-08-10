function [medoids,Distances] = nkmedoids(S, numCluster, modelInfo)

thau = modelInfo.thau;
size1 = modelInfo.size1;
distType = modelInfo.distType;
relative = modelInfo.relative;
subsWeight = modelInfo.subsWeight;
Itakura = modelInfo.Itakura;


medoids = cell(length(S), 1);
Distances = cell(length(S), 1);

for act=1:length(S)
    Stemp = S{act};
    row = size(Stemp.ind,2)-1;    
    k = numCluster;
    
    %Stemp.points =  Stemp.pca(:,1:4)'*Stemp.points;
%     if modelInfo.useLastCahange
%         Stemp.points =  Stemp.points(Stemp.idx(1:featSelect),:);
%     end
    if row<= numCluster 
        medoids{act}.ind = Stemp.ind;
        medoids{act}.points = Stemp.points;
        medoids{act}.probs = normalise(1:row,2);
        Dnew = zeros(row, row);
        for i=1:row
            A = Stemp.points(:,Stemp.ind(i):Stemp.ind(i+1)-1);        
            for j=i+1:row            
                B = Stemp.points(:,Stemp.ind(j):Stemp.ind(j+1)-1); 
                if distType == 0
                    L = nLevenstein(A,B,thau, relative, subsWeight);
                else
                    L = nDTW(A,B,thau, Itakura);
                end
                Dnew(i,j) = L(length(L));            
            end
        end
        Dnew = Dnew + Dnew';
        
        Distances{act}.dists = Dnew;
        Distances{act}.Labels = eye(row);
        continue;
    end
    if  row>=500
        b = ceil(0.1*row);
        e = row-b;
        Stemp.points = [Stemp.points(:,1:b), Stemp.points(:,e:row)];
        row = row-e+1+b;
    end
    
    col = row;
    
    if ismember(act, size1)
        Stemp.points = unique(Stemp.points','rows')';
        row = size(Stemp.points,2);
        
        Stemp.ind = 1:row+1;
        v = sum(Stemp.points,1);
        if distType ~= 1
            D = subsWeight*(bsxfun(@plus,v,v')-2*(Stemp.points'*Stemp.points))>thau;
        else
            D = bsxfun(@plus,v,v')-2*(Stemp.points'*Stemp.points);
        end
                
    else
        D = zeros(row, col);
        for i=1:row
            A = Stemp.points(:,Stemp.ind(i):Stemp.ind(i+1)-1);        
            for j=i+1:col            
                B = Stemp.points(:,Stemp.ind(j):Stemp.ind(j+1)-1);
                if distType == 0
                    L = nLevenstein(A,B,thau, relative, subsWeight);
                else
                    L = nDTW(A,B,thau, Itakura);
                end
                D(i,j) = L(length(L));            
            end
        end
        D = D + D';
    end 
     
    %v = dot(X,X,1);
    %D = bsxfun(@plus,v,v')-2*(X'*X);
    n = row;
    k = min(n,k);
    [~, label] = min(D(randsample(n,k),:),[],1);
    last = 0;
    while any(label ~= last)
        [~, index] = min(D*sparse(1:n,label,1,n,k,n),[],1);
        last = label;
        [~, label] = min(D(index,:),[],1);
    end
    %energy = sum(val);    
    index2 = 1;
    medPoints = [];
    medIndex = [];
    medHist = histc(label,1:k);
    for i=1:length(index)
        medIndex = [medIndex, index2];
        j = index(i);
        beginning = Stemp.ind(j);
        ending =  Stemp.ind(j+1)-1;
        medPoints = [medPoints, Stemp.points(:,beginning:ending)];
        index2 = index2 + ending-beginning+1; 
    end
    medIndex = [medIndex, index2];
    
    Dnew = zeros(row, size(medIndex,2)-1);
    if ismember(act, size1)       
        w = sum(medPoints,1);
        if distType ~=1
            Dnew =  subsWeight*(bsxfun(@plus,v',w)-2*(Stemp.points'*medPoints))>thau;
        else
            Dnew = bsxfun(@plus,v',w)-2*(Stemp.points'*medPoints);
        end        
    else
        for i=1:row
            A = Stemp.points(:,Stemp.ind(i):Stemp.ind(i+1)-1);        
            for j=1:size(medIndex,2)-1           
                B = medPoints(:,medIndex(j):medIndex(j+1)-1); 
                if distType == 0
                    L = nLevenstein(A,B,thau, relative, subsWeight);
                else
                    L = nDTW(A,B,thau, Itakura);
                end
                Dnew(i,j) = L(length(L));            
            end
        end
       % Dnew = Dnew + Dnew';
    end 
    Distances{act}.dists = Dnew;    
    medoids{act}.ind = medIndex;
    medoids{act}.points = medPoints;
    medoids{act}.probs = normalise(medHist,2);
    Distances{act}.Labels = sparse(1:n,label,1,n,k,n)*eye(k);
end