xl = 0;
xr = 1;
tb= 0;
tt= 1;
f = @(x) 0*x;
g = @(x) 2*pi*sin(pi*x);
l = @(t) 0*t;
r = @(t) 0*t;
c = 2; % = sqrt(4)
u_real = @(x,t) sin(pi*x).*sin(2*pi*t);
h = 0.05;
k = 0.9*h/c;
% h=(xr-xl)/M; k=(yt-yb)/N; m=M-1; n=N;
M = round((xr-xl)/h);
N = round((tt-tb)/k);
%M = 10; %h=0.1
%N = 500; %k=0.002
%N = 200; %k>0.003
w = finiteDifferenceWaveFunction(f, g, l, r, c, xl, xr, tb, tt, M, N);
[x, t] = meshgrid(1:(M+1), 1:(N+1));
x = x / (M+1);
t = t / (N+1);
w_real = u_real(x,t);

mesh(x,t,w'); % 3-D plot of solution w
axis([xl xr tb tt 0 1]);
%mesh(x,t,w_real); % 3-D plot of real w
%mesh(x, t, abs(w' - w_real))