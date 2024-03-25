function y = my_jac_gx(derivs,nxss,nyss,approx,h)
%This function computes the jacobian of g(x) fixing sigma=epsilon=1.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

lx = length(nxss);
x0 = nxss;
x0(lx) = 1;
derivs.hx(lx,lx) = 1;

[aa,bb]=size(derivs.gx);

derivs_c = zeros(aa,bb);

for n=1:lx-1

x_1 = x0;
y_0 = dr_gt(derivs,nyss,approx,(x0-nxss));


x_1(n) = x_1(n) + h;

y_new1 = dr_gt(derivs,nyss,approx,(x_1-nxss));
derivs_c(:,n) = (y_new1 - y_0)/h;

end

derivs_c(:,lx) = derivs.gx(:,lx);

y = derivs_c;