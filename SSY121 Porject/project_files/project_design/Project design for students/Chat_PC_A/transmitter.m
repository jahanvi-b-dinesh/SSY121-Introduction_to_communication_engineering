%clc; clear all; close all;
%Random binary string generated of length 432 for test purpose 
% The Message can allow 50 Character with in binary come up as 50 * 8 = 400
% so maybe 32 remain ??
%        info_bits = randi([0,1],1,432);
%fc = 4000;    
function transmitter(packet,fc)
info_bits=packet';
fsamp = 35000; % sampeling frequency
Tsymb = 0.006; % Symbol period
Tsamp = 1/fsamp; %Sampling period
M = 2; % Bits per symbol(QPSK)
fsymb = 1/Tsymb; % symbol frequency

preamble=[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1];% bits added to detect targeted signal
bits=[preamble,info_bits];

message=buffer(bits,2)';%transform bits into message, transpose and invert (buffer)
% Map constellation to message
messageDec = bi2de(message,'left-msb') +1;
constellation = [1+1i, -1+1i, -1-1i,1-1i]/sqrt(2);
x = constellation(messageDec);

% Sample up constellation points to aid(-kn)
fsfd = (1)*fsamp/fsymb; % Samples per symbol
x_up = upsample(x,fsfd);

% Create Pulse shape
alpha = 0.6; % variable for pulse shape
G = 0.006; % Variable for pulse shape, which determines the symbol time, bandwidth is (1+alpha)(2G)
span=4;% which determines number of peaks in the rrc_pulse
[rrc_p, t] = rtrcpuls(alpha,G,fsamp,span); % RRC with alpha and G factors

% Convolve x_up with the pulse shape to create a train of pulses
pulse_train = conv(rrc_p,x_up);

% Carrier signal
carr_cos =@(t) sqrt(2)*cos(2*pi*fc*t);
carr_sin =@(t) sqrt(2)*sin(2*pi*fc*t);
carrier_len = length(pulse_train);
t_new=0:Tsamp:(carrier_len-1)*Tsamp;%define sampling time of carrier

%put baseband signal on carrier
sig_send = real(pulse_train).*carr_cos(t_new) + imag(pulse_train).*carr_sin(t_new);

figure;
subplot(2,1,1); plot(sig_send);

soundsc(sig_send,fsamp)%send out the signal by soundcard
audiowrite("trans_sound.wav", sig_send,fsamp)
end

