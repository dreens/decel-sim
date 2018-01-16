%% Test Zero Lift
% I want to examine how hard it is to lift the zero with large electric
% fields.
Bp = 400; %T/m
B0 = 5000e-4;
E0 = 10e6; %V/m, 100kV/cm
B = @(x,y,z) [x -y B0/Bp]*Bp;
E = [0 E0 0];
Tnan = @(x,y,z) acos(sum(B(x,y,z).*E)/(norm(B(x,y,z))*norm(E0)));
iif = @(varargin) varargin{3-varargin{1}};
check = @(x) iif(isnan(x),0,x);
T = @(x,y,z) check(Tnan(x,y,z));
x = -1:.01:1;
y = x;
[xx, yy] = meshgrid(x,y);
V1 = zeros(size(xx));
V2 = V1;
for i=1:length(x)
for j=1:length(y)
    xs = x(i)*1e-3;
    ys = y(j)*1e-3;
    BB = norm(B(xs,ys,0));
    TT = T(xs,ys,0);
    ham = OH_Ham_Simple_SI(BB,norm(E),TT);
    energy = sort(eig(ham));
    V1(i,j) = energy(end);
    V2(i,j) = energy(end-1);
end
end

figure;
surf(xx,yy,V1)
figure;
surf(xx,yy,(V1 - V2)/6.626e-28);

%answer: its brutally hard- in fact impossible- to get rid of LZ loss when E-field is really big.
