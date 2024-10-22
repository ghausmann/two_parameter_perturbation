function y = btp_eval(b_bar,St,sol1,sol2,sol3) 
%SOE with PAC as auxiliary model
%Full version with income and interest rate shocks
%One-parameter perturbation
%This function evaluates the polynomial decision rule for bonds

%States
b_1 = St(1);
rt = St(2);
yt = St(3);
s = 1; %perturbation param

%First-order solution
hb = sol1(1);
hr = sol1(2);
hy = sol1(3);

%Second-order solution
hs2 = sol2(1);
hb2 = sol2(2);
hr2 = sol2(3);
hy2 = sol2(4);
hbr = sol2(5);
hby = sol2(6);
hry = sol2(7);

%Third-order solution
hb3 = sol3(1);
hr3 = sol3(2);
hy3 = sol3(3);

hb2r = sol3(4);
hbr2 = sol3(5);
hb2y = sol3(6);
hby2 = sol3(7);
hr2y = sol3(8);
hry2 = sol3(9);
hbry = sol3(10);
hs2b = sol3(11);
hs2r = sol3(12);
hs2y = sol3(13);

%Decision rule
y = b_bar + hb*(b_1-b_bar) + hr*rt + hy*yt ...
    + 0.5*( hb2*(b_1-b_bar)^2 + hr2*rt^2 + hy2*yt^2 + hs2*s^2) ...
    + hbr*(b_1-b_bar)*rt + hby*(b_1-b_bar)*yt+ hry*rt*yt ...
    + (1/6)*( hb3*(b_1-b_bar)^3 + hr3*rt^3 + hy3*yt^3 ) ...
    + 0.5*( hb2r*rt*(b_1-b_bar)^2 + hbr2*(b_1-b_bar)*rt^2 + hb2y*yt*(b_1-b_bar)^2 + hby2*(b_1-b_bar)*yt^2 + hr2y*yt*rt^2 + hry2*rt*yt^2 + hs2b*(b_1-b_bar)*s^2 + hs2r*rt*s^2 + hs2y*yt*s^2 ) ...
    + hbry*(b_1-b_bar)*rt*yt;