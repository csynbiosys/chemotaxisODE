function [ dmdt ] = chemosig(t,m,Lt,Lv)
%CHEMOSIG Summary of this function goes here

  %  2.081274522478973   for 100 constant
  %  2.423207034223378   for 200 constant
  %  2.535452158288703   for 250 constant
  %  2.704033304437643   for 350 constant
  %  2.738302210901962   for 375 constant
  %  3.1981   for 1000 constant
  %  3.4692       2000
  
    %Parameters
    KI_TAR_MEASP = 18.2;		%unit: um. Dissociation constant of MeAsp to inactive Tar receptor
    KA_TAR_MEASP = 3000;		%unit: um. Dissociation constant of MeAsp to active Tar receptor
    N_TAR = 6;					%Number of Tar receptors dimers in the complex

    ALFA_TAR_MEASP = 1.7;		%A parmeter in equation (2)
    M0_TAR_MEASP = 1;			%A parmeter in equation (3)

    KR = 0.005;					%unit: 1/s. Linear rate for methylation process
    KB = KR;					%unit: 1/s. Linear rate for demethylaion process

    H = 10.3;					%Hill coefficient of CheY-P response curve
    A_0 = 0.5;	
    
    %Ligand
%      L=100+100*(t>500);
    L=interp1(Lt,Lv,t);
    
    a = 1 ./ (1 + exp(N_TAR .* ((ALFA_TAR_MEASP .* (M0_TAR_MEASP - m))...
        - (log((1 + L./KA_TAR_MEASP)./(1 + L./KI_TAR_MEASP))))));
    %Update the mathylation level (equation 3)
    dmdt = (KR .* (1 - a) - KB .* a);

end

