function [DTW] = nDTW(A,B,thau, Itakura)
	r= size(B,2);
    c= size(A,2);
    dissim = bsxfun(@plus,sum(B,1)',sum(A,1))-2*(B'*A);
%     dissim = 1*(dissim>0);
    m = zeros(r, c);
    m(1,1) = dissim(1,1);
	m(1,:) = 0:c-1;
	m(:,1) = 0:r-1;
    for i=2:c
        m(1,i) = m(1,i-1)+dissim(1,i);
    end
    for j=2:r
        m(j,1) = m(j-1,1)+dissim(j,1);
    end
    

%     for i=2:r
% 		for j=2:c
% 		   %match = (nHamming(B(:,i-1),A(:,j-1))>thau);
% 		   %m(i,j) = min([m(i-1,j)+1, m(i,j-1)+1, m(i-1,j-1)+match]);
%            m(i,j) = min([m(i-1,j), m(i,j-1), m(i-1,j-1)])+dissim(i,j);
% 		end
%     end
w=round(Itakura*min(r,c));
    for i=2:r
        m(i,1:max(i-w,2)) = 1000;
        m(i,min(i+w,c):c)=1000;
		for j=max(i-w,2):min(i+w,c)		   
           m(i,j) = min([m(i-1,j), m(i,j-1), m(i-1,j-1)])+dissim(i,j);
		end
    end


 	DTW = m(r,1:c);	
%     DTW = normalise(m(r,1:c),2);
end