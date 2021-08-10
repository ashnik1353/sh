function [probs] = ncompLogProbs(distance, sigma, n)
	p2 = distance.*distance;
	p2 = (-1/(2*sigma*sigma))*p2;
	p2 = exp(p2);
	coef = (1/n)*(1/(((2*pi)^(0.5))*sigma));
	probs = sum(p2,1);
	probs = coef*probs;
    probs = log(probs);
end