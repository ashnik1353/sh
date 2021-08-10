function [sigma] = ncalcSigma(D , modelInfo)

    useBias = modelInfo.useBias;    
    
            phi = D.dists;
            if useBias
                phi = [ones(size(phi,1),1) phi];
            end
    
            K=size(phi,2);
            N = size(phi,1);
            lambdas = ones(1,K)/2;
            temp1 = zeros(1,K);
            pis = ones(1,K)/K;
    
            j=1;
            while (max(abs(lambdas-temp1),[],2)>0.1)&&(j<15)
                temp1 = lambdas;
                for i=1:K
                    pois(:,i) = poisspdf(phi(:,i),lambdas(i));
                end
                A = repmat(pis,N,1).*pois;
                gamas = normalise(A,2);
                NK = sum(gamas,1);
                lambdas = sum((gamas.*phi),1)./NK;
                pis = NK/N;
                j=j+1;
            end
            sigma.lambdas = lambdas;
            sigma.pis = pis;    
    
end