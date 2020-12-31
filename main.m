close all
clear all
clc

% use fmincon to do minimization - constant volume is a linear constraint

% problem parameters
M = 2;
a = 340;
L = 100;
V = 500;
VA = a*M;
H = 20000;
T = sqrt(2*H) - sqrt(L);

% required frequency components
fL = 10;
fH = 10000;

% solution parameters
N = roundn( fH*2*(3*L)/VA, 3 );
dx = (3*L)/N;
x = -L+dx/2:dx:2*L-dx/2;
dt = 5*dx*(500/V)/2;

S0 = zeros(size(x));
ind = find(x>=0 & x<=L);
S0(ind) = 16*V/3/L/pi*(4*x(ind)/L - 4*x(ind).^2/L/L).^(3/2);

u0 = init_u( S0, x, dx, L, V );
u = burgers(u0, dx, T, dt);
[P0, f] = calc_ft( u(end,:), x, N, VA, fL, fH );
I0 = calc_I( P0, f );

% uT = u(end,:);
% % duT = 0.0025*exp(-15*(x-x(2320)).^2);
% ind = find(abs(uT) > 1e-6);
% duT = 0.0025*exp(-15*(x-x(ind(end)-15)).^2);
% duT(ind(end):end) = 0.0;
% [P1, f] = calc_ft( uT+duT, x, N, VA, fL, fH );
% I = calc_I( P1, f );
% dI1 = I-I0
% 
% dx = x(2)-x(1);
% int_vec = ones(size(x))*dx/VA;
% int_vec(1) = int_vec(1)/2;
% int_vec(end) = int_vec(end)/2;
% psi0 = init_psi( P0, x, f, VA );
% dI2 = sum(int_vec.*psi0.*duT)

% I0 = calc_I( P, f, fL, fH );
% sens = zeros(size(x));
% for s = 3001:30:6000
%     tic
%     u0 = init_u( S0, x, dx, L, V );
%     S0(s) = 1.01*S0(s);
%     u = burgers(u0, dx, T, dt);
%     [P, f] = calc_ft( u(end,:), x, N, VA, L );
%     I = calc_I( P, f, fL, fH );
%     sens(s) = (I-I0)/0.01;
%     toc
%     s
% end
% plot(x(3001:30:6000),sens(3001:30:6000))

psi0 = init_psi( P0, x, f, VA );
disp('Init Done')

psi = adj_burgers( psi0, u, x, dx, T, dt );
disp('Adjoint Done')

grad = calc_grad( psi, x, dx, L );

% figure
% semilogy(f,abs(P(1:N/2+1)),'r*')
% ylabel('|P(f)|')
% xlabel('f')
% grid on
% xlim([10 10000])