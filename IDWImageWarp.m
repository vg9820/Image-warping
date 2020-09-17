function im2 = IDWImageWarp(im, psrc, pdst)

% input: im, psrc, pdst
%psrc is starting pt ,pdst is ending pt .
%in input is pt not arrow, psrc is equal to pdst.

%% basic image manipulations
% get image (matrix) size
[h, w, dim] = size(im);
% get the number of psrc point
[m, n] = size(psrc);
% P := input 2*m matrix  
% Q := output 2*m matrix
P = psrc';
Q = pdst';
im2 = im*0;
%% Inverse distance-weighted interpolation methods

% compute the D matrix(m*6) of fi
D0 = zeros(m, 4);
for i=1:m
    syms d1 d2 d3 d4 d5 d6;
    % D 3*2
    D = [d1,d2;d3,d4];
    E = 0;
    for j=1:m
        if j==i
            continue;
        end
        % w1:= wij
        w1 = d(i, P(:,j), P);
        % A is 3*1 vector
        fp = Q(:,i)+D*(P(:,j)-P(:,i))-Q(:,j);
        E = E + w1*fp'*fp;
    end
    %to make E(D) minimun,we require dE/d(dij)=0
    E1 = diff(E,d1)==0;
    E2 = diff(E,d2)==0;
    E3 = diff(E,d3)==0;
    E4 = diff(E,d4)==0;
    % why cannot D0(i,:)=solve(E1,E2,E3,E4,E5,E6,d1,d2,d3,d4,d5,d6) here?
    %debug 2!!!
    %class(solve(E1,E2,E3,E4,E5,E6,d1,d2,d3,d4,d5,d6)) is struct
    %D0(i,:) is double,cannot convey value of struct to double!
    [x1,x2,x3,x4] = solve(E1,E2,E3,E4,d1,d2,d3,d4);
    D0(i,:)=[x1,x2,x3,x4];
end
disp(D0);

%% use loop to negate image

for i=1:h
    for j=1:w
        f = 0;
        p = [i;j];
        % omega is the sum
        omega = 0;
        for l=1:m
            omega =omega + d(l,p,P);
        end
        
        for k=1:m
            %D=[D0(k,1),D0(k,2);D0(k,3),D0(k,4)]; 
            D=[1,0;0,1];
            fk=Q(:,k)+D*(p-P(:,k));
            w0= d(k, p, P);
            w1= w0/omega;
            f = f + w1*fk;
        end
        %prevent from overflow of (f1,f2)
        if(f(1)<1 || f(1)>h)
            continue
        end
        if(f(2)<1 || f(2)>w)
            continue
        end
        im2(round(f(1)),round(f(2)),:) = im(i,j,:);
            
    end
end



%% TODO: compute warpped image

