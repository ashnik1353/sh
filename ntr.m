function ksi = ntr(A, point, obsModel,alpha,betha)
k = size(alpha,1);

alpha = repmat(alpha,1,k);

P= repmat(point,1,k);
tempProbObs = prod(P.*obsModel+(1-P).*(1-obsModel),1);

tmp = log(tempProbObs)+betha';
betha = repmat(tmp,k,1);
ksi = exp(alpha+log(A)+betha);
        
