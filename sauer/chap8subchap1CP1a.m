f=@(x) 2*cosh(x);
l=@(t) 2*exp(2*t);
r=@(t) (exp(2) + 1)*exp(2*t-1);
D = 2;
xl = 0; xr = 1; tb = 0; tt = 1.1;
M = 10; %h=0.1
N = 500; %k=0.002
%N = 400; %k>0.003
w = forwardDiffHeat(f, l, r, D, xl, xr, tb, tt, M, N);
w_real_func = @(x,t) exp(2*t+x) + exp(2*t-x);
[x, t] = meshgrid(1:(M+1), 1:(N+1));
x = x / (M+1);
t = t / (N+1);
w_real = w_real_func(x,t);

mesh(x,t,w'); % 3-D plot of solution w
view(60,30);
axis([xl xr tb tt 0 30]);
%mesh(x,t,w_real); % 3-D plot of solution w
%mesh(x, t, abs(w' - w_real))