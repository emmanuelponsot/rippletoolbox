%% Examples of how to use the functions contained in the ripple toolbox
% E. Ponsot, Jan. 2020

close all;
clear;
clc;

%% Main parameters
omega_target=1;
w_target=8;

omega_temp=omega_target;
w_temp=w_target; 
phi_start_target_temp=0;

fmin=125;
f0=fmin;
fs=20000;
fmax=8000;
ntones = 400;
Targetlevel_dB=10;

duration_sound=0.5;
shape_target='squircle';
coef_smoothing_target=0;
plot_01=0;

%% 1 - Create and play a Downward ripple sound 

w_temp=8; % negative sign
d_temp=1; % modulation depth, between 0 and 1
plot_01=1; % plot the spectrogram (y=1, no=0)

[coef_m,target_wav_v]= fct_toolbox_Generate_Ripple(fmin, fmax, fs, duration_sound, ntones, f0, d_temp, omega_temp, w_temp, phi_start_target_temp, coef_smoothing_target, Targetlevel_dB,shape_target, plot_01);
soundsc(target_wav_v,fs);
imagesc(coef_m)
axis xy


%% 2 - Create and play an Upward ripple sound 

w_temp=-8; % negative sign
d_temp=1;

[coef_m,target_wav_v]= fct_toolbox_Generate_Ripple(fmin, fmax, fs, duration_sound, ntones, f0, d_temp, omega_temp, w_temp, phi_start_target_temp, coef_smoothing_target, Targetlevel_dB,shape_target, plot_01);
soundsc(target_wav_v,fs);



%% 3 - Example of a 3-interval ripple detection task
% you hear 3 noises, and only one is modulated : you have to detect which
% one it is !

ISI_n=0.3;
d_offset=0.5; % changing the modulation will result in an easier/more difficult task
order_v=[0 0 d_offset];
order_v=order_v(randperm(length(order_v)));
plot_01=1;
blank_v=zeros(1,floor(ISI_n*fs));
seq_v=[];

for ii=1:length(order_v)
    omega_temp=omega_target;
    w_temp=w_target;
    d_temp=order_v(ii)*d_offset;
    phi_start_target_temp=2*pi/2*randn(1);
    Targetlevel_dB=10;
    [~,target_wav_v]= fct_toolbox_Generate_Ripple(fmin, fmax, fs, duration_sound, ntones, f0, d_temp, omega_temp, w_temp, phi_start_target_temp, coef_smoothing_target, Targetlevel_dB,shape_target, plot_01);
    seq_v=[seq_v target_wav_v blank_v];
end

soundsc(seq_v,fs);
order_v

%% 3 - Example of a 3-interval ripple discrim task
% you hear 3 ripple noises, you have to detect which one is different from
% the other two !

% default params of the 2 identical stimuli
omega_target=1;
w_target=8;
d_temp=1;

% offset values on these params of the odd stimulus
omega_offset=0.5;
w_offset=+4;

order_v=[0 0 1];
order_v=order_v(randperm(length(order_v)));

blank_v=zeros(1,floor(ISI_n*fs));
seq_v=[];

for ii=1:length(order_v)
    omega_temp=omega_target+order_v(ii)*omega_offset;
    w_temp=w_target+order_v(ii)*w_offset;
    phi_start_target_temp=2*pi/2*randn(1);
    Targetlevel_dB=10;
    [coef_target_m,target_wav_v]= fct_toolbox_Generate_Ripple(fmin, fmax, fs, duration_sound, ntones, f0, d_temp, omega_temp, w_temp, phi_start_target_temp, coef_smoothing_target, Targetlevel_dB,shape_target, plot_01);
    seq_v=[seq_v target_wav_v blank_v];
end

soundsc(seq_v,fs);
order_v

%% 4 - Example of a 2-interval ripple target in ripple noise detection task
% Generate the Orientation noise and add it to the target
% % the parameters omega and w are not directly related to rate (Hz) and scale
% (ripples/oct), please see 
% for details about the procedure see E. Ponsot, L. Varnet, E. Daoud, N. Wallaert, S. Shamma, C. Lorenzi, & P. Neri (under review). 
% Mechanisms of spectrotemporal modulation detection for normal- and hearing- impaired listeners. (bioRxiv 2020/894667)
% here your task is to detect the interval that contains the upward ripple (as the one in example 1
% above) embedded in the ripple noise 


% main parameters you can play with

%%%%%%%%%
Targetlevel_play_dB=15; %  This will change the SNR between the ripple target and the ripple noise, i.e. the difficulty of the task!!
%%%%%%%%%

%%%%%%%%%
% params of the target ripple
spacing_target=10;
angle_degree_target=-45;
% equivalence in rate and scale values
[rate,scale] = fct_toolbox_conversion_AngleTORateScale(angle_degree_target,spacing_target,fmin,fmax,duration_sound)
%%%%%%%%%

%%%%%%%
% nuber of rotated components to create the ripple-noise
nb_angles=12;
%%%%%%%%




% other parameters
phi_start_list_v=[0 pi/2 pi 3*pi/2];
phi_start_target = phi_start_list_v(1); % value chosen by hand to center the target
sd_noise_dB=3;
fs=20000;
ntones = 400;
duration_sound=0.5;
ISI_n=0.1;
fmin=250;
f0=fmin;
fmax=8000;
d=1;
angle_degree_surround=angle_degree_target;
spacing_surround=spacing_target;
Surroundlevel_dB=NaN; 

nb_noises_n=1;
plot_01=1;
spacing_nontarget=spacing_target;
angle_degree_nontarget=angle_degree_target;

shape_target='squircle';
shape_noise='squircle';
shape_surround='shortsurround'; %'squircle'
coef_smoothing_target=1;
coef_smoothing_noise=1;

% create the coefficients of the spectrotemporal envelopes of the ripple
% noise first, because it may take some time
[coef_noise1_generated_m, scale_rate_values1_m, level_noise1_generated_m, phi_noise1_generated_m, phi_start1_target_v]= fct_toolbox_generateCoefMatrix_ORIENTnoise_final(fmin, fmax, fs, duration_sound, ntones, f0, d, spacing_target,phi_start_list_v, angle_degree_target,nb_angles, sd_noise_dB, nb_noises_n);
[coef_noise2_generated_m, scale_rate_values2_m, level_noise2_generated_m, phi_noise2_generated_m, phi_start2_target_v]= fct_toolbox_generateCoefMatrix_ORIENTnoise_final(fmin, fmax, fs, duration_sound, ntones, f0, d, spacing_nontarget,phi_start_list_v, angle_degree_nontarget,nb_angles , sd_noise_dB, nb_noises_n);

phi_start_surround=phi_start1_target_v;


% Create the sounds and play
for ss=1:nb_noises_n
    
    % generate the noise only
    coef_noise_play1_m=zeros(size(coef_noise1_generated_m,2),size(coef_noise1_generated_m,3));
    coef_noise_play1_m(:,:)=coef_noise1_generated_m(ss,:,:);
    phi_start_target1=phi_start1_target_v(ss);
   
    level_noise1_m=zeros(1,size(level_noise1_generated_m,2));
    phi_noise1_m=zeros(1,size(phi_noise1_generated_m,2));
    level_noise1_m(:)=level_noise1_generated_m(ss,:);
    phi_noise1_m(:)=phi_noise1_generated_m(ss,:);
    
    % generate the noise + target
    coef_noise_play2_m=zeros(size(coef_noise2_generated_m,2),size(coef_noise2_generated_m,3));
    coef_noise_play2_m(:,:)=coef_noise2_generated_m(ss,:,:);
    phi_start_target2=phi_start2_target_v(ss);
    
    level_noise2_m=zeros(1,size(level_noise2_generated_m,2));
    phi_noise2_m=zeros(1,size(phi_noise2_generated_m,2));
    level_noise2_m(:)=level_noise2_generated_m(ss,:);
    phi_noise2_m(:)=phi_noise2_generated_m(ss,:);
    
    % create the sequences
    seq_wav=[];
    rand_order_n=randi(2)-1;
    spacing_surround=spacing_target;
    
    % Generate the ORIENT sounds
    % and plot the construction of the stimuli in the two intervals
    [targetPlusNoise1_wav_v, coef_noise1_m, coef_target1_m, coef_total1_m, target1_wav_v, Noise1_wav_v, realSNR1_dB] = fct_toolbox_RippleInORIENTnoise_Surround(fmin, fmax, fs, duration_sound, ntones, f0, d, angle_degree_target, spacing_target,phi_start_target1, coef_smoothing_target, coef_noise_play1_m, coef_smoothing_noise, Targetlevel_play_dB,shape_target,shape_noise,shape_surround, angle_degree_surround, spacing_surround, phi_start_surround, Surroundlevel_dB,plot_01);
    [nontargetPlusNoise2_wav_v, coef_noise2_m, coef_target2_m, coef_total2_m, target2_wav_v, Noise2_wav_v, realSNR2_dB] = fct_toolbox_RippleInORIENTnoise_Surround(fmin, fmax, fs, duration_sound, ntones, f0, d, angle_degree_nontarget, spacing_nontarget,phi_start_target2, coef_smoothing_target, coef_noise_play2_m, coef_smoothing_noise, Targetlevel_play_dB,shape_target,shape_noise,shape_surround, angle_degree_surround, spacing_surround, phi_start_surround, Surroundlevel_dB,plot_01);
    
   % randomize the order of target and non-target
    if rand_order_n == 0
        seq_wav=[seq_wav targetPlusNoise1_wav_v zeros(1,ISI_n*fs) Noise2_wav_v]; % target is in the 1st interval
        disp('target was first');
        soundsc(seq_wav,fs);
    else
        seq_wav=[seq_wav Noise2_wav_v zeros(1,ISI_n*fs) targetPlusNoise1_wav_v]; % target is in the 2nd interval
        disp('target was second');
        soundsc(seq_wav,fs);
    end
end


