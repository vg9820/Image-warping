function im2 = RBFImageWarp(im, psrc, pdst)

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
%% Radial basis functions methods
% compute r
r=zeros(1,m);
for i =1:m
    r1 = zeros(1,m);
    for j =1:m
        if i==j
            %make sure this r(j) cannot be chosen
            r1(j)=99999;
            continue
        end
        r1(j) = norm(P(:,i)-P(:,j));
    end
    r(i) = min(r1);
end
       
% compute alpha vector
for i=1:m
    a(i) = sym(['a',num2str(i)]);
    b(i) = sym(['b',num2str(i)]);
    E1(i) = sym(['E1',num2str(i)]);
    E2(i) = sym(['E2',num2str(i)]);
end
%set u = 1
u=1;
for i=1:m
    E1(i) = P(1,i);
    E2(i) = P(2,i);
    for j=1:m
        R = power(norm(P(:,i)-P(:,j))^2+r(j)^2,u/2);
        E1(i) = E1(i) + a(j)*R;
        E2(i) = E2(i) + b(j)*R;
    end
    E1(i) = E1(i)==Q(1,i);
    E2(i) = E2(i)==Q(2,i);
end
%solve the equation to get alpha
A1 = solve(E1,a);
A2 = solve(E2,b);
A = zeros(2,m);
%here A1,A2 are struct, we need get value inside
x1 = fieldnames(A1);
x2 = fieldnames(A2);
for i=1:m
    A(1,i) = A1.(x1{i});
    A(2,i) = A2.(x2{i});
end

%% use loop to negate image
for i=1:h
    for j=1:w
        p = [i;j];
        f = p;
        % omega is the sum
        sum = 0;
        for k=1:m
            R = power(norm(p-P(:,k))^2+r(k)^2,u/2);
            sum = sum + A(:,k)*R;
        end
        f = f + sum;
        
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