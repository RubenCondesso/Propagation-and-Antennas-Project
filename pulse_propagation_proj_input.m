%--------------------------------------------------------------------------
% INPUT FILE containing data for the script PULSE PROPAGATION
%
% Carlos R. Paiva
%
% 27 April 2017
%
% Edited for PA project by Miguel Carreiro and Ruben Condesso
%
% 5 December 2017

%--------------------------------------------------------------------------
% DEFAULT CHIRP PARAMETER AND DISTANCE
%--------------------------------------------------------------------------
C = 0;                         % Chirp parameter (C=0, C=-6, C=6)
L_km = 50
%--------------------------------------------------------------------------

pulse = input('Please enter the integer value for PULSE between 1 and 3:     ')

% There are 3 different classes of input pulses:
%
% [PROJETO] pulse = 1   --->   chirped Gaussian pulse
% [PROJETO] pulse = 2   --->   chirped super-Gaussian pulse
% [PROJETO] pulse = 3   --->   chirped hyperbolic secant pulse

if pulse == 1
   C = input('Insert chirp parameter:   ')
   L_km = input('Insert distance in km:   ')
   input_pulse = ['Gaussian pulse with chirp ' num2str(C) ' and distance ' num2str(L_km)]
elseif pulse == 2
   m_Input = input('Insert m value:  ')
   C = input('Insert chirp parameter:   ')
   L_km = input('Insert distance in km:   ')
   input_pulse = ['super-Gaussian pulse with m = ' num2str(m_Input) ', chirp ' num2str(C) ' and distance ' num2str(L_km)]
elseif pulse == 3
   C = input('Insert chirp parameter:   ')
   L_km = input('Insert distance in km:   ')
   input_pulse = ['hyperbolic secant pulse with chirp ' num2str(C) ' and distance ' num2str(L_km)]
else
   input_pulse = 'Error! No signal encountered! Closing...'
   return;
end
 
%
%--------------------------------------------------------------------------
%
% Input data
%
T_0_ps = 10
T_0 = T_0_ps*(1e-12);
c = 299792458;
lambda = 1550e-9;
lambda_nm = 1550

%
%--------------------------------------------------------------------------
%
% Input data: dispersion
%
beta_2_usual_units = -1;
beta_3_usual_units = 0.08;

beta_2 = beta_2_usual_units*1e-27;
beta_3 = beta_3_usual_units*1e-39;
%
%--------------------------------------------------------------------------

if beta_2 == 0
   % Em principio nao sera utilizado no projeto
   zeta_prime = L_km/L_D_prime_km
   L_D_prime_km = 1e-3*T_0^3/abs(beta_3)
   ztotal = zeta_prime;
else
   % unidades corretas?
   L_D_km = 1e-3*(T_0^2)/abs(beta_2)
   zeta = L_km/L_D_km
   ztotal = zeta;
   kappa = abs(beta_3)/(abs(beta_2)*T_0)
end

%
%--------------------------------------------------------------------------
%
% Input data for computations
%

if pulse == 1 
   t0 = 100;
   t1 = 10;
elseif pulse == 2
   t0 = 500;
   t1 = 15;
elseif pulse == 3
   t0 = 500;
   t1 = 10;
elseif pulse == 4
   t0 = 100;
   t1 = 10;
elseif pulse == 5
   t0 = 100;
   t1 = 10;
else
   t0 = 50;
   t1 = 4;
end

W_max = 20;
W_min = -W_max;

p = 14;  N_t = 2^p;          % Number of points for the FFT        
t = linspace(-t0,t0,N_t);    % Vector defining the TIME window

M = 2*2^10;
for m = 1:M
    T(m) = t(7*2^10+m);
end

%--------------------------------------------------------------------------
% INPUT PULSE
%--------------------------------------------------------------------------
% VALOR DE C AQUI!!

a_max = 1;                     % Maximum amplitude

%--------------------------------------------------------------------------

if pulse == 1

% Gaussian pulse (to = 100, t1 = 10, p = 14)
% m=1
a = exp(-0.5*(1+1i*C)*t.^2);
%
f_D = @(t) exp(-t.^2);
D = integral(f_D,-inf,inf);
f_N = @(t) t.^2.*exp(-t.^2);
N = integral(f_N,-inf,inf);
sigma_0 = sqrt(N/D);
sigma_0_ps = sigma_0*T_0_ps

%--------------------------------------------------------------------------
% nao usar para impulso gaussiano simples
elseif pulse == 2

% Super-Gaussian pulse (to = 500, t1 = 20, p = 14)
% Project pulse (T_0 = 10 ps, ...)
% m=3, m=10
%
m = m_Input;
a = exp(-0.5*(1+i*C)*(t).^(2*m));
f_D = @(t) exp(-t.^(2*m));
D = integral(f_D,-inf,inf);
f_N = @(t) t.^2.*exp(-t.^(2*m));
N = integral(f_N,-inf,inf);
sigma_0 = sqrt(N/D);
sigma_0_ps = sigma_0*T_0_ps

%--------------------------------------------------------------------------

elseif pulse == 3

% Hyperbolic secant pulse (to = 100, t1 = 15, p = 14)
%
a = sech(t).*exp(-0.5*1i*C*t.^2);
%
f_D = @(t) sech(t).^2;
D = integral(f_D,-inf,inf);
f_N = @(t) t.^2.*sech(t).^2;
N = integral(f_N,-inf,inf);
sigma_0 = sqrt(N/D);
sigma_0_ps = sigma_0*T_0_ps

%--------------------------------------------------------------------------

end
