function [probs] = ncompProbs(distance, sigmas, size1  , modelInfo)


diffSize1 = modelInfo.diffSize1;
useBias = modelInfo.useBias;


        if size1 && diffSize1
                probs = ones(1,size(distance,2));
        else
            
            sigma = sigmas.lambdas;
            coef = sigmas.pis;
            if useBias
                distance = [ones(1,size(distance,2));  distance];
            end
            n = size(distance,1);
            p2 = zeros(n,size(distance,2));
            for i=1:n
                p2(i,:)= poisspdf(distance(i,:),sigma(i));
            end
            p2 = repmat(coef',1,size(distance,2)).*p2;
            probs = sum(p2,1);
        
  
        end


end