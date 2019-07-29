%
%
%   Find connected subclusters
%   find their abundances too
%
%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ clu, sabu ] = find_connected( C )


    c = ff(C);
    
    [val, abu] = values2(c);
    
    [sabu, ind] = sort(abu, 'descend');
    
    sval = val(ind);
    clu = c;
    
    for i=1:length(sval)
        clu(find(c==sval(i)))=i;
    end

end

%
%   function cl = ff(C)
%
%	Forest-fire algorithm with the matrix of links
% 
%	v is the vector of values in matrix L
%	c is the vector of cluster numbers
%
%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cl = ff(C)

    C0 = (C>0);
    N = size(C,1);
	cl=zeros(N,1);                  % List of clusters
	cc=1;                           % Current cluster number
	n=0;
    
	for cn=1:N			%
            
        fraction_done = cn/N;
        %disp(char(repmat(8,1,50)))
            
            
        if cl(cn)==0				% Not classified yet
    							% Start new cluster

                v = sparse(cn, 1, 1, N, 1);
                
                
                fv = cn;
                dn=1;
                
                while(dn>0)
                    nv = C0*v+v;
                    fnv = find(nv);
                    dn = length(fnv)-length(fv);
                    v=nv;
                    fv = fnv;
                end

                ind = find(v);
                
                cl(ind)=cc;
                
                cc = cc+1;
                fraction_done = cn/N;
        end
    end
    return
end

%   (c) Alex Koulakov (akula@cshl.edu) 2019 

function [val, abu] = values2(y)

    
    yy=sort(y(:));
    dy=diff(double(yy));
    ind=find([1; dy; 1]);
    val=yy(ind(1:(length(ind)-1)));
    abu=diff(ind);
    
    return
end