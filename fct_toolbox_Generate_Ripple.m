function [coef_target_m,modtarget_wav_v,carrier_wav_v]= fct_toolbox_Generate_Ripple(fmin, fmax, fs, duration_sound, ntones, f0, d, omega_target, w_target, phi_start_target, coef_smoothing_target, Targetlevel_dB,shape_target, plot_01)

%%
% This function creates a Ripple audio signal 
%
% INPUTS:
% fmin, fmax: spanning the fmin to fmax region (in Hz)
% fs : sampling frequency in Hz
% duration_sound : in seconds
% ntones : number of tones (log-spaced across the [fmin-fmax] region) for
% the carrier
% f0 : the first tone freq.
% d : amplitude modulation depth of the ripple (between 0 and 1)
% omega : ripple density in cycle/octave - no sense here
% w_target : temporal rate in Hz
% phi_start_target : starting phase of the ripple
% coef_smoothing target:  % if == 0 => normal ripples; if >0, it corresponds to the smoothing coefficient for the gaussian window applied to the disc centered in the image (in units of radius) to make the gabor shape
% Targetlevel_dB: the level of the target
% plot_01: 0 <=> don't plot / 1 <=> plot the targeted spectrogram (freq on
% a log axis)
%
%
% OUTPUTs:

% coef_target_m: 2D matrix of amp coefs for the enveloppe of the target
% (theoretical log-freq spectrogram)
% modtarget_wav_v: wav file of the ripple signal = enveloppe*carrier
% carrier_wav_v: wav file of the carrier signal alone
%
%%
%  E. Ponsot, 
% last seen on 07/01/2020

%% Parameters

% create the frequenciy vector and the amplitudes/phases of the different components of the carrier
t_v = linspace(0, duration_sound, duration_sound * fs);
f_v=zeros(1,ntones);
for ii=1:ntones
    f_v(ii) = fmin * (fmax / fmin).^((ii - 1)/(ntones - 1));
end

gamma_v = rand(1, ntones); % random amplitudes
phi_v = 2*pi*rand(1,ntones); % random phases
x_v=log(f_v/f0);


%% Windowing
%  Select only a disc (or a square) with smoothed contours in the middle
radius_small=floor(ntones/8); % the radius of the disc is not in argument of the function but it can me manually changed here
[mask_rescaled_small_v] = fct_toolbox_mask(ntones,duration_sound,fs,radius_small,coef_smoothing_target,shape_target); %  for the target

% to make a rotation of a gabor
coef_target_m=zeros(ntones, duration_sound * fs);
for ii=1:ntones
    coef_target_m(ii,:)=10^((Targetlevel_dB)/20)*(1+d*sin(2*pi*(w_target*t_v + omega_target*x_v(ii))+phi_start_target));
end


% windowing the target
if coef_smoothing_target>0
    coef_target_m(:,:)=  coef_target_m(:,:).*mask_rescaled_small_v;
end


%% Plot
% make the plots of the spectrograms
if plot_01 ~= 0
    figure
    imagesc(abs(coef_target_m))
    axis xy
    colorbar
    title('ripple spectrogram')
end
colormap('gray')

%% apply the overall matrix of coefficients (the enveloppe) to the carrier
coef_ii=zeros(1,duration_sound * fs);
modtarget_wav_v=zeros(1,duration_sound * fs);
carrier_wav_v=zeros(1,duration_sound * fs);
for ii=1:ntones
    coef_ii(:)=coef_target_m(ii,:);
    s_ii=gamma_v(ii)*1./sqrt(f_v(ii))*sin(2*pi*f_v(ii)*t_v+phi_v(ii));
    modtarget_wav_v=modtarget_wav_v+coef_ii.*s_ii;
    carrier_wav_v=carrier_wav_v+s_ii;
end


%%
% normalize rms of the target_wav_v output
modtarget_wav_v=modtarget_wav_v./sqrt(sum(abs(modtarget_wav_v).^2)/length(modtarget_wav_v));
carrier_wav_v=carrier_wav_v./sqrt(sum(abs(modtarget_wav_v).^2)/length(modtarget_wav_v));

end
