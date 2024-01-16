function [audio_recorder] = receiver(fc)

%-----------------------------------------------------------
% Veriable Defination
%-----------------------------------------------------------

fs = 35000;         % Sampling frequency
Tsymb = 0.006;      % Symbol period
Ts = 1/fs;          % Sampling period
M = 2;              % Bits per symbol(QPSK)
fsymb = 1/Tsymb;    % Symbol frequency

%-----------------------------------------------------------
% Veriable for Pluse Shaping
%-----------------------------------------------------------

alpha = 0.6;        % variable for pulse shape
G = 0.006;          % Variable for pulse shape, which determines the symbol time, bandwidth is (1+alpha)/(2*G)
span=4;             % which determines number of peaks in the rrc_pulse


%-----------------------------------------------------------
% use rtrcpuls function to generate RRC
%-----------------------------------------------------------
% RRC with alpha and G factors
[rrc_p, t] = rtrcpuls(alpha,G,fs,span); 

%-----------------------------------------------------------
% Preamble Sequence of bits
%-----------------------------------------------------------
% bits added before packet for detection
preamble=[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1];
% transform preamble into message(2 bits per symble)
p_message=buffer(preamble,2)'; 

%-----------------------------------------------------------
% Map constellation to message
%-----------------------------------------------------------
mess_Dec = bi2de(p_message,'left-msb') +1;
constellation = [1+1i, -1+1i, -1-1i,1-1i]/sqrt(2);
pre_map_c = constellation(mess_Dec);

%-----------------------------------------------------------
% Upsample the presample-message and Convolved with RRC
%-----------------------------------------------------------

up_preamble=upsample(pre_map_c,round(fs/fsymb));
pulse_preamble=conv(up_preamble,rrc_p);

%-----------------------------------------------------------
% Efficient Synchronization by taking 60% of preamble for detection
%-----------------------------------------------------------

len_preamble=length(pulse_preamble);
len_pre_used=floor(3*len_preamble/5);
len_pre_left=ceil(2*len_preamble/5);
pulse_preamble=pulse_preamble(1:len_pre_used);

%-----------------------------------------------------------
% The Preamble pluse train modulated onto a carrier signal 
% by multiplying with carrier Cosine and Sine waveform
%-----------------------------------------------------------

carrier_len = length(pulse_preamble);
% Define sampling time of carrier signnal
tt_new=0:Ts:(carrier_len-1)*Ts;
carr_cos =@(t) sqrt(2)*cos(2*pi*fc*t);
carr_sin =@(t) sqrt(2)*sin(2*pi*fc*t);
preamble_carrier = real(pulse_preamble).*carr_cos(tt_new) + imag(pulse_preamble).*carr_sin(tt_new);


%-----------------------------------------------------------
% Store variables in function so that it can be used in audioTimerFcn(recObj, event, ~)
%-----------------------------------------------------------

inputval = inputValues(fs,fc, fsymb, rrc_p, constellation, preamble_carrier);

% Create the recorder
audio_recorder = audiorecorder(fs,24,1);

%-----------------------------------------------------------
% ADD USER DATA FOR CALLBACK FUNCTION
%-----------------------------------------------------------

audio_recorder.UserData.receive_complete = 0;
audio_recorder.UserData.pack  = [];         % Allocate for data package
audio_recorder.UserData.pwr_spect = [];     % Allocate for PSD
audio_recorder.UserData.const = [];         % Allocate for constellation
audio_recorder.UserData.eyed  = [];         % Allocate for eye diagram
audio_recorder.UserData.inputval= inputval; % Input vars for the callback function

%-----------------------------------------------------------
% Attach callback function
%-----------------------------------------------------------

time_value = 0.1;
% Use cells to specify function and its arguments
set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn); 
% Start recording
record(audio_recorder); 
end


%-----------------------------
%-----------------------------
% % CALLBACK FUNCTION %%
%-----------------------------
%-----------------------------

function audioTimerFcn(recObj, event, ~)

% Store recieved audio signal in "pshape" 
pshape= (getaudiodata(recObj)).'; 

Ts = 1/recObj.UserData.inputval.fs; % Sampling period
fc = recObj.UserData.inputval.fc;   % Carrier frequency 

%-----------------------------
% % Signal Detection %%
%-----------------------------
% Looking for a convolution that passes threshhold --> data being received
% Compute the Energy of Preamble 

Ecorr = sum(recObj.UserData.inputval.preamble_carrier.^2);

% Correlation output - Received signal with flip version preamble divided
% by the energy of preamble to normalize the correlation output
corr = conv(pshape, fliplr(recObj.UserData.inputval.preamble_carrier))/Ecorr;

% Find the maximum value and index of maximum value in the correlation
[maxval,index] = max(corr);

% Initialization of data signal array.
dataSig = [];

% Comparing (maximum value / Energy of received signal) with threshold
% 0.05, If data is received, wait for package and save. wait for 0.5 sec

if maxval/Ecorr > 0.05   
    disp('found data :)')
    pause(0.5)
    
    % Store recieved audio signal in "pshape"
    pshape = getaudiodata(recObj);
    corr = conv(pshape, fliplr(recObj.UserData.inputval.preamble_carrier))/Ecorr;

    % Find the preamble
    [maxval,index] = max(corr); 
    
    % Remove the preamble and find the starting point of targeted signal
    startData = index -length(recObj.UserData.inputval.preamble_carrier)+1;
    % End point of targeted signal
    endData = startData+(216+13)*(recObj.UserData.inputval.fs/recObj.UserData.inputval.fsymb)+ length(recObj.UserData.inputval.pulseShape)-1;
    % Store the targetd signal
    dataSig = pshape(startData:endData);
    

    %-----------------------------
    % Passband -TO- baseband
    %-----------------------------
    % Define sampling time
    tt_new=0:Ts:(length(dataSig)-1)*Ts;
    
    carr_cos =@(t) sqrt(2)*cos(2*pi*fc*t);
    carr_sin =@(t) sqrt(2)*sin(2*pi*fc*t);
    data_bb = carr_cos(tt_new).*dataSig.' + 1i*carr_sin(tt_new).*dataSig.';
     
    %-----------------------------
    %  Matched filter
    %-----------------------------        
    
    Errc = sum(recObj.UserData.inputval.pulseShape.^2);
    data_mf = conv(data_bb, conj(fliplr(recObj.UserData.inputval.pulseShape)))/Errc;
    % Signal processed by match filter, but still have tail trancient
    data_mf = data_mf(length(recObj.UserData.inputval.pulseShape):length(data_mf));
    
    fsfd = round(recObj.UserData.inputval.fs/recObj.UserData.inputval.fsymb);

    %-----------------------------
    % %Symbol --> samples
    %-----------------------------       
   
    data_samp = downsample(data_mf,fsfd);
    firstPar = data_samp(1);
    % Take the targeted message
    data_samp = data_samp(14:216+13);

    %-----------------------------
    %  Fix attenuation and angle offset
    %-----------------------------   
    
    angOffset = angle(firstPar) - angle(1 -1i);
    ampOffset = abs(firstPar);
    data_samp = (exp(-1i*angOffset)/ampOffset).*data_samp;

    %-----------------------------
    % % Samples --> bits
    %----------------------------- 
    
    recovered_message = data_samp;
    d_message=[];
    for ii=1:length(recovered_message)
        if real(recovered_message(ii))>=0&imag(recovered_message(ii))>0
            d_message(ii)=1;
        elseif real(recovered_message(ii))<0&imag(recovered_message(ii))>=0
            d_message(ii)=2;
        elseif real(recovered_message(ii))<=0&imag(recovered_message(ii))<0
            d_message(ii)=3;
        else
            d_message(ii)=4;
        end
    end
    d_message=d_message-1;
    d_bits=de2bi(d_message,'left-msb')';
    d_bits=d_bits(:)';
        
    %-----------------------------
    %  Record the estimated bits
    %----------------------------- 
    recObj.UserData.pack = d_bits;
    
    %-----------------------------
    %  Save the sampled symbols
    %----------------------------- 
    recObj.UserData.const = data_samp;

    %-----------------------------
    %  Provide the matched filter output for the eye diagram
    %----------------------------- 
    recObj.UserData.eyed.fsfd = fsfd;
    % Remove the transcient caused by match filter
    recObj.UserData.eyed.r = data_mf(1:end-9*fsfd+1);   
    
    %-----------------------------
    %  Compute the PSD and save it. 
    %  Note that it has to be computed on 
    %  the BASE BAND signal BEFORE matched filtering
    %----------------------------- 
    % pwr_spect.f will be normalized frequencies
    [pxx, f] = pwelch(data_mf,1024,768,1024, recObj.UserData.inputval.fs); 
    
    % Shift to be centered around fs
    f = fftshift(f); 
    
    % Center to be around zero
    f(1:length(f)/2) = f(1:length(f)/2) - recObj.UserData.inputval.fs; 
    
    % Shift, normalize and convert PSD to dB
    p = fftshift(10*log10(pxx/max(pxx))); 
    recObj.UserData.pwr_spect.f = f;
    recObj.UserData.pwr_spect.p = p;
      
    recObj.UserData.receive_complete = 1;
    stop(recObj)       
end
   
end