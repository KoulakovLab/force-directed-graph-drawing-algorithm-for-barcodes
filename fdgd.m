%
%   Force directed graph drawing algorithm
%   [x,y,z] = fdgd(C, q, U0)
%   graph connectivity, charge, and parabolic confinement strength
%   
%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%
%

function [x,y,z] = fdgd(C, q, U0)

    clear costfunction2
    clear persistent
    
    
    l=1;

    % make sure connection matrix is symmetric

    CC = (C+C')/2;
    
    if 1

    % initial positions of points

    v = sum(CC);
    L = diag(v)-CC;
    
    
    [V,D] = eig(L);
    
    D=diag(D);
    
    [sD, ind]=sort(D, 'descend'); 
    
    P = V(:, ind(2:4) );   % projector
    
    x = P(:,1);
    y = P(:,2);
    z = P(:,3);
    
    dr = pdist([x,y,z]);
    mdr = mean(dr);
    x=x/mdr;
    y=y/mdr;
    z=z/mdr;
    
    hist(dr, 100)
    
    else
        
       x = randn(size(CC,1),1)*100;
       y = randn(size(CC,1),1)*100;
       z = randn(size(CC,1),1)*100;
       
    end

    %scatter3(x0,y0,z0, 30, 'k')
    
    X0=[x;y;z];
    [f, X] = cg2('costfunction2',X0,1e-4,1e-4,10,1,CC,0, U0);    % no charge    first
    
    [f, X] = cg2('costfunction2',X,1e-4,1e-8,10,1,CC,q, U0);    
    
    
    
    N = length(X)/3;
    
    x = X(1 : N);
    y = X((N+1) : (2*N));
    z = X((2*N+1) : end);

    
    
    
    doAnneal = 0;
    
    
    if doAnneal


    E = 0;
    T = 1e1;
    Ebest = E;


scale = 0.01;

Es = ones(1,doAnneal*size(C,1));
Esbest = Es;

step = 1;

for MCStep=1:doAnneal
    %MCStep
    for n=1:size(C,1)

        dx = randn*scale;
        dy = randn*scale;
        dz = randn*scale;
        
        dE = dcostfunction(C, x, y, z, q, l, n, dx, dy, dz);
    
        if rand<exp((-dE)/T)
        
            x(n)=x(n)+dx;
            y(n)=y(n)+dx;
            z(n)=z(n)+dx;
            E=E+dE;
        
            scale = scale*(1+0.01);
        
        else

            scale = scale*(1-0.01);
        
        end
    
        if Ebest>E
            Ebest = E;
            xbest = x;
            ybest = y;
            zbest = z;
        end
    
        %size(E)
        %size(Es(step))
        
        Es(step)=E;
        Esbest(step)=Ebest;
        step=step+1;
        
    end
    
    T = T*(1 - 5/doAnneal); 
    
    
    if ~mod(MCStep,100)
        
        MCStep
    
    %subplot(2,1,1)
    
    plot(Es, 'b')
    hold on
    plot(Esbest, 'r')
    plot(step, T, '.m')
    hold off
    drawnow
    
    scale
    dE
    
    end
end



x = xbest;
y = ybest;
z = zbest;
    
end

% doAnneal








end



function dE = dcostfunction(C, x, y, z, q, l, n, dx, dy, dz)

    r_old = sqrt((x-x(n)).^2+(y-y(n)).^2+(z-z(n)).^2);
    r_new = sqrt((x-x(n)-dx).^2+(y-y(n)-dy).^2+(z-z(n)-dz).^2);

    % first do springs
    
    ind = find(C(:,n));
    Eold = sum((r_old(ind) - l).^2)/2;
    Enew = sum((r_new(ind) - l).^2)/2;
    
    fi_new = -q ./ r_new;
    fi_old = -q ./ r_old;
    
    fi_old(isnan(fi_old)) = 0;
    fi_old(isinf(fi_old)) = 0;

    fi_new(isnan(fi_new)) = 0;
    fi_new(isinf(fi_new)) = 0;
    
    Eold = Eold + sum(fi_old);
    Enew = Enew + sum(fi_new);
    
    dE = Enew-Eold;
end

%
%   Compute cost function
%

function E = costfunction2(In, C, q, U0)

    persistent c
    
    [N1,N2]=size(C);
    if length(c)~=(N1-1)*N1/2
        c=[];
    end
    
    if isempty(c)
        N=size(C,1);
       c = zeros(1, N*(N-1)/2);
       n=1;
       for i=1:N
           for j=(i+1):N
               c(n) = C(i,j);
               n=n+1;
           end
       end
    end

    N = length(In)/3;
    x = In(1 : N);
    y = In((N+1) : (2*N));
    z = In((2*N+1) : end);

    if 0
        
    X = x * ones(1, size(z,1));
    Y = y * ones(1, size(z,1));
    Z = z * ones(1, size(z,1));
    
    dX = X-X';
    dY = Y-Y';
    dZ = Z-Z';

    R = sqrt(dX.^2 + dY.^2 + dZ.^2);
    E = sum(sum((R-1).^2/2 .* C))/2;
    E = E + q * sum(sum( 1 ./ (R+1e-10) ))/2 ;

    else
        %c
        %keyboard
       %R = squareform(pdist([x,y,z])); 
       r = pdist([x,y,z]);
       
       try
            E = sum((r-1) .* (r-1) .* c / 2) + q * sum(1 ./ (r+1e-10));
            mx = mean(x);
            my = mean(y);
            mz = mean(z);
            
            E = E + U0 * sum((x-mx).^2 + (y-my).^2 + (z-mz).^2) / 2;
            
       catch
            size(r)
            size(c)
            keyboard
       end
    end
    

end

function G = grad_costfunction2(In, C, q, U0)

    N = length(In)/3;
    x = In(1 : N);
    y = In((N+1) : (2*N));
    z = In((2*N+1) : end);

    if 1
    X = x * ones(1, size(z,1));
    Y = y * ones(1, size(z,1));
    Z = z * ones(1, size(z,1));
    
    dX = X-X';
    dY = Y-Y';
    dZ = Z-Z';

    R = sqrt( dX.^2 + dY.^2 + dZ.^2);

    else
        R = squareform(pdist([x,y,z]));
        
    end
    
%     Gx = (R-1) .* dX ./ (R+1e-10) .* C - q * dX ./ (R+1e-10).^3; 
%     Gy = (R-1) .* dY ./ (R+1e-10) .* C - q * dY ./ (R+1e-10).^3; 
%     Gz = (R-1) .* dZ ./ (R+1e-10) .* C - q * dZ ./ (R+1e-10).^3; 

    R1 = 1 ./ (R+1e-10);
    R3 = R1 .* R1 .* R1;
    RRR = (C .* (R - 1) .* R1 - q * R3);
    
     Gx = dX .* RRR; 
     Gy = dY .* RRR; 
     Gz = dZ .* RRR; 


    
    gx = sum(Gx,2);
    gy = sum(Gy,2);
    gz = sum(Gz,2);
    
            mx = mean(x);
            my = mean(y);
            mz = mean(z);
            
            gx = gx + U0 * (x-mx);
            gy = gy + U0 * (y-my);
            gz = gz + U0 * (z-mz);
    
    G = [gx; gy; gz];
    
    if ~isempty(find(isnan(G)))
       keyboard 
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% conjugate gradient algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%

%
%       function [f, x] = cg(funfcn,X0,ddx,prec,max_jump,displ,varargin)
%

function [f, x] = cg2(funfcn,X0,ddx,prec,max_jump,displ,varargin)

	tau = (sqrt(5)-1)/2;
	tic
	x = zeros(size(X0));
	x0 = X0;
	g = x;
	g_prev = x;
	p = x;
	p_prev = x;
	
	n = max(size(x));
	
	
	x = x0;
	step = 0;
	
	dddx = ddx;
	MaxPrec=prec
	
	if prec>1
		MaxStep=round(prec);
		MaxPrec=0;
	else
		MaxStep=1e10;
	end
	
	while(1)
	
        for k=1:n
	
		%f0=feval(funfcn,x0,varargin{:});
        %eval(['f0=' funfcn '(x0,varargin{:});']);
        f0=costfunction2(x0, varargin{:});
        
        
		if displ ~= 0, 
            fprintf('\n Step = %i, F(x) = %10.15g', step, f0); 
            toc;
            
            if ~mod(step, 100)
                N=length(x)/3;
                xx = x(1 : N);
                yy = x((N+1) : (2*N));
                zz = x((2*N+1) : end);
               show_graph( varargin{1}, xx, yy, zz, 'r', 'k' ) 
               title(sprintf('step=%g', step))
               drawnow
            end
        end;
		
	
		% Find the gradient
        
        if 0
			
		for i=1:n
			x(i) = x(i) + dddx;
			%f=feval(funfcn,x,varargin{:});
            %eval(['f=' funfcn '(x,varargin{:});']);
            f=costfunction2(x, varargin{:});
            
			g(i) = (f-f0)/dddx;
			x(i) = x(i) - dddx;		
        end
        
        else
            g = grad_costfunction2(x, varargin{:});
        end
	
		% Choose the new search direction
	
		if (k > 1)
			c1 = sum((g-g_prev).*g);
			c2 = sum((g-g_prev).*p_prev);
			if c2==0
				p=-g;
			else
				beta = c1/c2;
				p = -g+beta*p_prev;
			end
		else
			p = -g;		
		end
		
		% Normalize the direction vector
		
		n_p=norm(p);
		if n_p~=0
			p = p/norm(p);
		else
			p=randn(size(p));
			p=p/norm(p);
		end
	
		p_prev = p;
		g_prev = g;
		
		% Interval location
		
		% positive direction
	
		D=max_jump;
	
		a1 = D*(tau);
		%a1 = max_jump/2;
		x = x0 + a1*p;
		%f1 = feval(funfcn,x,varargin{:});
		%eval(['f1=' funfcn '(x,varargin{:});']);
        f1=costfunction2(x, varargin{:});
        
        if 1
        	count2=0;
		while(1)
		
			f = f1;
			a1 = a1*2;
			D=D*2;
			x = x0 + a1*p;
			%f1 = feval(funfcn,x,varargin{:});
            %eval(['f1=' funfcn '(x,varargin{:});']);
            f1=costfunction2(x,varargin{:});
            
			count2=count2+1;
			if ((f1>=f)||(count2>5e1)), break, end;	
		end
        end
    
		% negative direction
		
		a2 = - D*(1-tau);
		%a2 = - max_jump/2;
		x = x0 + a2*p;
		%f2 = feval(funfcn,x,varargin{:});
		%eval(['f2=' funfcn '(x,varargin{:});']);
        f2=costfunction2(x,varargin{:});
        
        	if 0
			while(1)
				
				f = f2;
				a2 = a2*2;
				x = x0 + a2*p;
				%f2 = feval(funfcn,x,varargin{:});
                %eval(['f2=' funfcn '(x,varargin{:});']);
                f2=costfunction2(x,varargin{:});
		
				if (f2>=f), break, end;	
			end, 
		end
		
	
		% Interval reduction
		
		a = a2;
		b = a1;
		c = a + (1-tau)*(b-a);
		d = b - (1-tau)*(b-a);
		
		x = x0+c*p; 
        %Fc = feval(funfcn,x,varargin{:});
        %eval(['Fc=' funfcn '(x,varargin{:});']);
		Fc=costfunction2(x,varargin{:});
		
        x = x0+d*p; 
        %Fd = feval(funfcn,x,varargin{:});
        %eval(['Fd=' funfcn '(x,varargin{:});']);
		Fd=costfunction2(x,varargin{:});
		
		count1=0;
		while (1)
		
			if (Fc<Fd)			
				b=d; d=c;
				c=a+(1-tau)*(b-a);
				Fd=Fc;
				x = x0+c*p; 
				%Fc = feval(funfcn,x,varargin{:});		
                %eval(['Fc=' funfcn '(x,varargin{:});']);
                Fc=costfunction2(x,varargin{:});
			else
				a=c; c=d;
				d=b-(1-tau)*(b-a);
				Fc=Fd;			
				x = x0+d*p; 
				%Fd = feval(funfcn,x,varargin{:});
                %eval(['Fd=' funfcn '(x,varargin{:});']);
                Fd=costfunction2(x,varargin{:});
			end
			
			count1=count1+1;
			
			if (abs(Fc-Fd)<=MaxPrec)||(count1>1e3), break; end;
		end
		
		x0=x;
		%f=feval(funfcn,x0,varargin{:});
        %eval(['f=' funfcn '(x0,varargin{:});']);
		f=costfunction2(x0,varargin{:});
		%dddx = min([ddx abs(c)*1e-3])	
		step = step+1;
		
		%if (step>=MaxStep)
		if ( abs(f-f0) < MaxPrec )||(step>=MaxStep)
			return;
		end;
		
		
	end, end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








