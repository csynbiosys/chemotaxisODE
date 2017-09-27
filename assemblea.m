a=[];
KI_TAR_MEASP = 18.2;		%unit: um. Dissociation constant of MeAsp to inactive Tar receptor
KA_TAR_MEASP = 3000;		%unit: um. Dissociation constant of MeAsp to active Tar receptor
N_TAR = 6;					%Number of Tar receptors dimers in the complex

ALFA_TAR_MEASP = 1.7;		%A parmeter in equation (2)
M0_TAR_MEASP = 1;			%A parmeter in equation (3)

KR = 0.005;					%unit: 1/s. Linear rate for methylation process
KB = KR;					%unit: 1/s. Linear rate for demethylaion process

H = 10.3;					%Hill coefficient of CheY-P response curve
A_0 = 0.5;
L=interp1(Lt,Lv,t);
atmp=[];
for i=1:length(t)
    atmp= 1 ./ (1 + exp(N_TAR .* ((ALFA_TAR_MEASP .* (M0_TAR_MEASP - y(i)))...
            - (log((1 + interp1(Lt,Lv,t(i))./KA_TAR_MEASP)./(1 + interp1(Lt,Lv,t(i))./KI_TAR_MEASP))))));
    a=[a;atmp];
end