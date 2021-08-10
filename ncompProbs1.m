function [probs] = ncompProbs1(distance, sigma , coef)
    
 	%p2 = distance.*distance;
    
    p2 = poisspdf(distance,sigma);

    p2 = repmat(coef',1,size(distance,2)).*p2;
    probs = sum(p2,1);

end