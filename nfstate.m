function [s] = nfstate(A,states)

% s = 1;
% [M,~] = find(A==1);
% [~,~,T] = intersect(M,states);
% if ~isempty(T)
%     s = T(1) + 1;
% end
Y = [];
for i=1:length(A)
    if A(i) == 1
        Y = [Y , i];
    end
end
if isempty(Y)
    s=1;
else
% [Y,I] = find(A>0);
    s=1;
    for k=1:length(Y)
        for j=1:length(states)
            if Y(k)==states(j)
                s=j+1;
                break;
            end
        end
    end
end


end
