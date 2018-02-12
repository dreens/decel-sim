%% OH ground state in E, B fields without hyperfine
function H = OH_Ham_Simple_SI(Bf,E,beta)

sb = sin(beta);
cb = cos(beta);

uB=9.27401*1e-24; 
g = 0.9355; %1.4032 is 3/2 * the gJ factor for OH, 0.9355
uE=3.33564*1e-30;
dE = 1.67;
h=6.62607*1e-34;
c=2.99792458*1e8;
LD = h*(1.667358e9);


L = [-LD/2*eye(4) zeros(4) ; zeros(4) LD/2*eye(4)];
Z =  diag([-3 -1 1 3 -3 -1 1 3])*g*uB*Bf/2;

S1 = diag([0.6 .2 -.2 -0.6])*uE*E*dE*cb;
S2 = [0 3 0 0 ; 3 0 4 0 ; 0 4 0 3 ; 0 0 3 0];
S2 = -0.2*sqrt(S2)*uE*E*dE*sb;
S3 = S1 + S2;
S = [zeros(4) S3 ; S3 zeros(4)];

H = L + Z + S;

end