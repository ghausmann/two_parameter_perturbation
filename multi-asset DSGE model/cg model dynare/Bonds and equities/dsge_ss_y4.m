function y = dsge_ss_y4(pf,m_params,A,B)
%This function evaluates equation E.26 in Appendix E .

phi = m_params(1);
alpa = m_params(2);
d = m_params(3);
betta = m_params(4);

a = A(1);
as = A(2);

b = B(1);
bs = B(2);


P = ( alpa + (1-alpa)*pf^(1-phi) )^(1/(1-phi));
Ps = ( (1-alpa) + alpa*pf^(1-phi) )^(1/(1-phi));

C = (1/P)*( (1-d) + a*d + as*pf*d + (1-betta)*(b*P + bs*Ps) );
Cs = (1/Ps)*( (1-d)*pf + (1-a)*d + (1-as)*pf*d  -(1-betta)*(b*P + bs*Ps) );


y = alpa*(P^phi)*C + (1-alpa)*(Ps^phi)*Cs - 1;