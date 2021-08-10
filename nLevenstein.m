function [Levenstein] = nLevenstein(A,B,thau,relative,subsWeight)
	r= size(B,2)+1;
    c= size(A,2)+1;
    m = zeros(r, c);
	m(1,:) = 0:c-1;
	m(:,1) = 0:r-1;
    Hammings = bsxfun(@plus,sum(B,1)',sum(A,1))-2*(B'*A);
    
    Hammings = subsWeight*(Hammings > thau);   

    %Hammings = B'*(1-A)>0;
    for i=2:r
		for j=2:c
		    m(i,j) = min([m(i-1,j)+1, m(i,j-1)+1, m(i-1,j-1)+Hammings(i-1,j-1)]);
		end
    end
    if relative
        n = ones(1,c-1);
        for i = 1:c-1
            n(i) = max(r-1,i);
        end
        Levenstein = m(r,2:c)./n;
    else
        Levenstein = m(r,2:c);	
    end
end