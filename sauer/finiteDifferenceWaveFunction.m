% Program 8.1 Forward difference method for heat equation
% input: space interval [xl,xr], time interval [yb,yt],
% number of space steps M, number of time steps N
% output: solution w
% Example usage: w=heatfd(0,1,0,1,10,250)
function w=finiteDifferenceWaveFunction(f, g, l, r, c, xl,xr,yb,yt,M,N)
h=(xr-xl)/M; k=(yt-yb)/N; m=M-1; n=N;
sigma=c*k/(h);

a = diag(2-2*sigma*sigma*ones(m,1));
a = a + diag(sigma*sigma*ones(m-1,1), 1);
a = a + diag(sigma*sigma*ones(m-1,1),-1); % define matrix a

lside=l(yb+(0:n)*k); rside=r(yb+(0:n)*k);

xes = (1:m)*h;
w(:,1) = f(xl+xes);
w(:,2) = 0.5*a*f(xl+xes)' + k*g(xes)';
w(:,2) = w(:,2) + 0.5*sigma*sigma*[lside(1);zeros(m-2,1);rside(1)]; % initial conditions
for j=2:n
    w(:,j+1) = a*w(:,j) - w(:,j-1) + sigma*sigma*[lside(j);zeros(m-2,1);rside(j)];
end
w=[lside;w;rside]; % attach boundary conds
