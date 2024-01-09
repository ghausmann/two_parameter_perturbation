function y=mykron1(x1,x2)
%
% Copyright (C) 2024 Guillermo Hausmann Guil

[s1,~] = size(x1);
[s2,~] = size(x2);

xx1 = kron(x1,ones(s2,1));
xx2 = repmat(x2,s1,1);
%keyboard;
y = xx1.*xx2;