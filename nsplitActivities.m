function [S] = nsplitActivities(modelInfo, Labels, Segments, trainFeatMat)
actList = modelInfo.actList;
numAct = modelInfo.numAct;

S = cell(numAct,1);

parfor i=1:length(actList),
    start = 1;
    index = 1;       
    points=[];    
    for n=1:length(Labels),        
        idx = find(Labels{n}==actList(i));
        if (~isempty(idx))
            seg = Segments{n}(idx);
            points = [points, trainFeatMat{n}(:,idx)];        
            idx2 = find(seg(2:end)-seg(1:end-1)~=0);               
            if (~isempty(idx2))
                idx2 = idx2 + start;
            end
            index = [index, idx2, size(points,2) + 1];
            start = size(points,2) + 1;           
        end
    end       
	S{i}.ind = index;
%     S{i}.points = points;
    %[~,idx] = sort(std(points'),'descend');
%     [~,idx] = sort(sum(points,2),'descend');
    idx = 1:size(points,1);
    S{i}.idx = idx;
    %S{i}.pca = pca(points');
    
    S{i}.points = points;   
    
end

end