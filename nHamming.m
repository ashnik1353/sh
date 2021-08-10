function [hammingDistance] = nHamming(s1,s2)
	hammingDistance = sum(s1+s2-2*s1.*s2);
end