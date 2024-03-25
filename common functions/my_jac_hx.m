function y = my_jac_hx(derivs,nxss,approx,h)
%This function computes the jacobian of g(x) fixing sigma=epsilon=1.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

lx = length(nxss);

x0 = nxss;
x0(lx) = 1;
derivs.hx(lx,lx) = 1;

[aa,bb]=size(derivs.hx);

derivs_c = zeros(aa,bb);

for n=1:lx-1

    x_1 = x0;
    x_1(n) = x_1(n) + h;
    x_new0 = dr_ht(derivs,nxss,approx,(x0-nxss));
    x_new1 = dr_ht(derivs,nxss,approx,(x_1-nxss));
    derivs_c(:,n) = (x_new1 - x_new0)/h;

end

derivs_c(:,lx) = derivs.hx(:,lx);

y = derivs_c;