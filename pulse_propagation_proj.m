% Pulse propagation in linear SMF (single-mode fibers):
% Numerical simulation using FFT
%
% Carlos R. Paiva
%
% 21 April 2017
%
% Edited for PA project by Miguel Carreiro and Ruben Condesso
%
% 5 December 2017

clear all            
close all         

%--------------------------------------------------------------------------
% Read the input DATA in file PULSE_PROPAGATION_DATA
%--------------------------------------------------------------------------
%
pulse_propagation_proj_input;
%
%--------------------------------------------------------------------------

% N_z = Número de passos na distância
%
N_z = 500;                     
%
z = linspace(0,ztotal,N_z);        % Vector distância com 'Nz' posições


Ts = t(2)-t(1);                    % Separação entre amostras
Ws = 2*pi/Ts;                      % Largura total da janela espectral [rad/s]
W = Ws*[0:N_t/2-1 -N_t/2:-1]/N_t;  % Vector de frequências [rad/s]
WE = fftshift(W);

U(1,:) = a;                       % Armazena na matriz U o impulso inicial
A = fft(a);                       % Coefs. de Fourier na posição z = 0
Aux = fftshift(A)*Ts;             % Espectro invertido
F(1,:) = Aux([N_t:-1:1]);         % Armazena na matriz F o espectro ordenado

for i = 2:N_z                                     % Ciclo para cada posição de z
    coeff_1 = sign(beta_2)*1i/2*W.^2;
    coeff_2 = sign(beta_3)*kappa*1i/6*W.^3;
    coeff_exp = (coeff_1+coeff_2).*z(i);
    Az = A.*exp(coeff_exp);
    a = ifft(Az);                                 % Impulso no domínio temporal:  a(z,t)
    U(i,:) = a;                                   % Armazena na matriz U o impulso temporal
    Aux = fftshift(Az)*Ts;                        % Espectro invertido
    F(i,:) = Aux([N_t:-1:1]);                     % Armazena na matriz F o espectro ordenando
end

for m = 1:M
    for i = 1:N_z
        U_prime(i,m) = U(i,7*2^10+m);
    end
end

I_max = max(abs(U(end,:).^2));
S_max = max(abs(F(1,:).^2));

denominator = trapz(t,abs(U(end,:).^2));
mom_2 = trapz(t,t.^2.*abs(U(end,:).^2))/denominator;
mom_1 = trapz(t,t.*abs(U(end,:).^2))/denominator;
sigma_p = sqrt(mom_2-mom_1^2);
sigma_p_ps = sigma_p*T_0_ps
eta = sigma_p/sigma_0

figure(1);         
plot(t,abs(U(1,:)).^2,'b',t,abs(U(end,:).^2),'r','LineWidth',2);        % Gráfico 2D do impulso inicial/final
ax = gca;
ax.FontSize = 20;
ax.Color = [255 255 204]/255;
ax.FontName = 'Times New Roman';
grid on
h = legend('input pulse','output pulse');
c = h.TextColor;
h.TextColor = [1 0 1];
xlabel('Time');
ylabel('Optical Intensity');
axis([-t1 t1 0 a_max]); 
hold on
plot([-t1 t1],[1 1]*I_max,'--r','LineWidth',1)

figure(2);         
plot(WE,abs(F(1,:).^2),'b','LineWidth',2);        % Gráfico 2D do espectro do impulso inicial/final
ax = gca;
ax.FontSize = 20;
ax.Color = [204 204 255]/255;
ax.FontName = 'Times New Roman';
grid on
xlabel('Angular Frequency');
ylabel('Energy Spectral Density');
axis([W_min W_max 0 S_max])

[X,Y] = meshgrid(T,z);                              % Cria uma grelha rectangular de dimensões (Nt,Nz)
                                                    % usada por 'mesh' 
figure(3);              
mesh(X,Y,abs(U_prime).^2);                                % Gráfico da evolução do impulso
xlabel('Time');
ylabel('Distance');
zlabel('Optical Intensity');

figure(4);
mesh(X,Y,abs(U_prime).^2);
xlabel('Time');
ylabel('Distance');
zlabel('Optical Intensity');
view(30,45);                                         % view(azimute,elevação);   por defeito azimute=-37.5
                                                     % elevação=30
axis tight;

%[XE,YE] = meshgrid(WE,z);
%figure(5);
%mesh(XE,YE,abs(F).^2);
%xlabel('Angular Frequency');
%ylabel('Distance');
%zlabel('Energy Spectral Density');
%axis tight;
                                