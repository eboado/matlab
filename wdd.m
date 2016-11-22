clc;
clear;
alpha = 0.5;
mu = [0.2 0.5 0.7];
Sigma_T = [0.5; 0.8; 1];
Sigma_S = [0.1 0.0 0.0; 0.3 0.1 0.1; 0.1 0.3 0.3];
Qe = [1.5; 0; 0.2]; 
delta = .1;

x_i = 0;
x_f = 2;
n = (x_f-x_i)/delta;

x_forward = linspace(x_i,x_f,n);

Psi_half = zeros(1,n);

Psi = zeros(1,n);
S = zeros(1,n);
S_prev = zeros(1,n);
Psi_reverse = zeros(1,n);

RelativeError = 1;

for k = 1:3
PS = zeros(1,n);
PSR = zeros(1,n);
for h =1:numel(mu);
    
while RelativeError > 1e-4
    
    for i = 1:length(Psi)
        if k == 1
        Psi_half(1) = 0.5;
        else
        Psi_half(1) = 0;
        end
        S(i) = Sigma_S(k,1) * Psi(i) + Sigma_S(k,2) * Psi(i) + Sigma_S(k,3) * Psi(i) + Qe(k);
        Psi(i) = ( S(i) + (2*abs(mu(h))*Psi_half(i))/((1+alpha)*delta) ) / ( Sigma_T(k) + ( (2*abs(mu(h))) / ((1+alpha)*delta) ) );
        Psi_half(i+1) = ( (2*Psi(i)) / (1+alpha) ) - ( ( (1-alpha)/(1+alpha) ) * Psi_half(i));
    end
    
    Psi_half_reverse(1) = Psi_half(i+1);    
    
    
    for i = 1:length(Psi_reverse)
        S_reverse(i) = Sigma_S(k,1) * Psi(i) + Sigma_S(k,2) * Psi(i) + Sigma_S(k,3) * Psi(i) + Qe(k);
        Psi_reverse(i) = ( S_reverse(i) + (2*abs(mu(h))*Psi_half_reverse(i))/((1-alpha)*delta) ) / ( Sigma_T(k) + ( (2*abs(mu(h))) / ((1-alpha)*delta) ) );
        Psi_half_reverse(i+1) = ( (2*Psi_reverse(i)) / (1-alpha) ) - ( ( (1+alpha)/(1-alpha) ) * Psi_half_reverse(i));
    end
    
    
    Diff = abs(S_prev - S);
    AE = 0;
    for i = 1:length(Diff)
        AE = AE + (Diff(i))^2;
    end
    AbsoluteError = sqrt(AE);
    
    Denominator = 0;
    for i = 1:length(S)
        Denominator = Denominator + (S(i))^2;
    end
    Den = sqrt(Denominator);
    
    RelativeError = AbsoluteError/Denominator;
    S_prev = S;
end

PS = PS + Psi;
PSR = PSR + Psi_reverse;
Psi = PS;
Psi_half = zeros(1,n);
Psi_half_reverse = zeros(1,n);
end

Psi = PS + PSR;
plot(x_forward, Psi)
hold on
str = sprintf('Graph at \\Delta %f and \\alpha %f', delta, alpha);
grid on
title(str);
legend({'Group 1','Group 2','Group 3'});
end
