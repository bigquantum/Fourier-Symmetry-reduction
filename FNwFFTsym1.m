
%------------------------------------------------------------------------
%--- FN system in 1D. SR with FOURIER FFT-------------------------------------
%------------------------------------------------------------------------
clc
close all;
clear all;

%%

% x -> N,i

%--dimensions...........................................................
p.N = 500;
p.Lx = 13;
p.T = 3000;
p.dx = p.Lx/(p.N-1); % step size
p.diff = 0.001;

x = linspace(-p.Lx/2,p.Lx/2,p.N);
kx = (2*pi/p.Lx) * [-p.N/2:p.N/2-1]';
kx = fftshift(kx);

%prameters from ODE
p.alpha = 0.1;
p.beta = 0.5;
p.gama = 1;
p.epsilon = 0.001;
p.delta = 0.0;

moving = 0;
method = 1;

%% Reduced symmetry

addpath('../getRoots')

file_name = './DATA/icondFFT01.mat';
dataFN = load(file_name);

u = dataFN.uout;
v = dataFN.vout;
p = dataFN.p;
x = dataFN.x;
kx = dataFN.kx;
p.dt = 0.01;
scale = 100;
p.T = 20000*scale;
p.Tau = 2000*scale; % Recording time
p.M = 20000;
p.rate = p.Tau/(p.M/50); % Sample rate

%Initial condition
y0 = [u' ; v'];
u0 = y0(1:p.N);
v0 = y0((p.N+1):(2*p.N));
    
% Symmetry stuff
k0 = 2*pi/p.Lx; % kx(2);
uhat = fft(u);
uhat0 = uhat(2);
phi0 = angle(uhat0);

t = 0;
s = 0;
d = 0;
distance = [];
cdot = [];
intu = [];
intv = [];
APL = [];

while t <= p.T
    
    [y0,s] = RK4fftShift(y0,p,kx,s,phi0,moving,method);
    
    u = y0(1:p.N);
    v = y0((p.N+1):2*p.N);
    
    if (mod(t,p.rate) == 0)
        plot1D(x,u,v,p);
        disp(t)
    end
    
    t = t + p.dt*scale;
    
    % Save data
    d = d + s*kx(2)*p.Lx/(2*pi);
    distance = [distance d];
    cdot = [cdot (s*kx(2)*p.Lx/(2*pi))/p.dt];

end

uout = y0(1:p.N);
vout = y0((p.N+1):2*p.N);

% uout = ifft(y(end,1:p.N));
% vout = ifft(y(end,(p.N+1):2*p.N));

figure(607)
plot(x,u0,'go-')
hold on
plot(x,v0,'ko-')
hold on
plot(x,uout','bo-')
hold on
plot(x,vout','rx-')
grid on



