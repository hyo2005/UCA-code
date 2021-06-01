clear all;
close all;
clc; 


%%
tic
Nt = 2;    % number of transmit antennas
Nr = 32;    % number of receive antennas
L   = 4;    % channel order
M   = 2;    % Number of multipaths (assumption: M  = L)   
K  = 64;        % OFDM subcarriers
sigmax2=1;
Pxp=0.15*sigmax2;
Np      = 4;      % Pilot symbols
F  = dftmtx(K);
FL  = F(:,1:L);


%% Signal Generation
% we use the Zadoff-Chu sequences
U = 1:2:7;
ZC_p = [];
for u = 1 : Nt
    for k = 1 : K
        ZC(k,u) = sqrt(Pxp) * exp( ( -1i * pi * U(u) * (k-1)^2 ) / K );
    end
    ZC_p = [ZC_p; ZC(:,u)];
end

%% Channel generation 
fading=[0.8,0.6,0.4,0.2;0.9,0.7,0.5,0.3];
delay=[110,240,180,900;570,350,742,236]*10^(-8);%s
DOA=[pi/2,pi/4,pi/6,pi/8;pi/3,pi/5,pi/7,pi/9];
d_nor=0.5;
R_nor=0.5*d_nor/sin(pi/Nr);

%% Derivative

dev_h_fading=[];
dev_h_delay=[];
dev_h_angle=[];

dev_h_fading_UCA=[];
dev_h_delay_UCA=[];
dev_h_angle_UCA=[];

for Nr_index=1:Nr;
%ULA
Br_fading = spec_chan_derive_fading_ULA(fading,delay,DOA,d_nor,Nr_index,L,M,Nt);
dev_h_fading=[dev_h_fading; transpose(Br_fading)];

Br_delay = spec_chan_derive_delay_ULA(fading,delay,DOA,d_nor,Nr_index,L,M,Nt);
dev_h_delay=[dev_h_delay; transpose(Br_delay)];

Br_angle = spec_chan_derive_angle_ULA(fading,delay,DOA,d_nor,Nr_index,L,M,Nt);
dev_h_angle=[dev_h_angle; transpose(Br_angle)];
    
%UCA    
Br_fading_UCA = spec_chan_derive_fading_UCA(fading,delay,DOA,R_nor,Nr_index,Nr,L,M,Nt);
dev_h_fading_UCA=[dev_h_fading_UCA; transpose(Br_fading_UCA)];

Br_delay_UCA = spec_chan_derive_delay_UCA(fading,delay,DOA,R_nor,Nr_index,Nr,L,M,Nt);
dev_h_delay_UCA=[dev_h_delay_UCA; transpose(Br_delay_UCA)];

Br_angle_UCA = spec_chan_derive_angle_UCA(fading,delay,DOA,R_nor,Nr_index,Nr,L,M,Nt);
dev_h_angle_UCA=[dev_h_angle_UCA; transpose(Br_angle_UCA)];

end

%% Derivation of $h$ w.r.t. (bar{h},tau,alpha) %% channel specular parameters
%G_ULA
G = [dev_h_fading,dev_h_delay,dev_h_angle]; 

%G_UCA
G_UCA = [dev_h_fading_UCA,dev_h_delay_UCA,dev_h_angle_UCA]; 
%% ------------------------------------------------------------------------
 
X = [];
for ii = 1 : Nt
    %X        = [X (F'/sqrt(K))*diag(ZC(:,ii))*FL];
    X        = [X diag(ZC(:,ii))*FL];
end

%% CRB 
%Loop SNR
SNR = -10:5:30;
for snr_i = 1 : length(SNR)
    
    sigmav2 = 10^(-SNR(snr_i)/10);

  %usual
    X_nga=kron(eye(Nr),X);
    Iop      = Np*X_nga'*X_nga / sigmav2;
    CRB_op(snr_i) = abs(trace(pinv(Iop)));
    
  %spec_ULA  
    Iop_spec = Np*((-1)/(sigmav2)^2)*G*G'* X_nga'*sigmav2*eye(Nr*K)*X_nga*G*G';
    CRB_op_spec(snr_i) = abs(trace(pinv(Iop_spec)));
      
  %spec_UCA 
    Iop_spec_UCA = Np*((-1)/(sigmav2)^2)*G_UCA*G_UCA'* X_nga'*sigmav2*eye(Nr*K)*X_nga*G_UCA*G_UCA';
    CRB_op_spec_UCA(snr_i) = abs(trace(pinv(Iop_spec_UCA)));      
end

%figure

semilogy(SNR,CRB_op,'-b*','LineWidth',1)
hold on; semilogy(SNR,CRB_op_spec,'-dr','LineWidth',1);
hold on; semilogy(SNR,CRB_op_spec_UCA,'-g>','LineWidth',1);

grid on
ylabel('CRB')
xlabel('SNR(dB) - UCA')
legend('usual OP','spec OP UCA', 'spec OP UCA')
title(' ')
%axis([-10 20 1e-4 1e2])
hold on;