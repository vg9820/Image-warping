function y = d(i, p, P)
%caculate the distance between p and Xi
d = norm(p - P(:,i));
% choose u
u = 4;
% get di(p)
y = 1/power(d,u);
end

