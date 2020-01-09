function [targetPlusNoise_wav_v, coef_noise_m, coef_target_m, coef_total_m, target_wav_v, Noise_wav_v, realSNR_dB]= fct_toolbox_RippleInORIENTnoise_Surround(fmin, fmax, fs, duration_sound, ntones, f0, d,angle_degree_target, spacing_target,phi_start_target, coef_smoothing_target, coef_noise_generated_m, coef_smoothing_noise, Targetlevel_dB,shape_target,shape_noise,shape_surround, angle_degree_surround, spacing_surround, phi_start_surround, Surroundlevel_dB, plot_01)

%% fct_toolbox_RippleInORIENTnoise_Surround
%%
% This function creates a signal made of a Gabor/Ripple target (defined by its orientation) embedded in
% a Spectro-Temporal Modulation noise made of gabor/ripple signals (created
% earlier) using pre-defined noise envelopes saved in
% coef_noise_generated_m. You can add a spectrotemporal surround as well.
%
% INPUTS:
% fmin, fmax: spanning the fmin to fmax region (in Hz)
% fs : sampling frequency in Hz
% duration_sound : in seconds
% ntones : number of tones (log-spaced across the [fmin-fmax] region) for
% the carrier
% f0 : the first tone freq.
% d : amplitude modulation depth of the ripple (between 0 and 1)
% angle_degree_target : orientation of the target (in degrees)
% spacing_target : distance parameter controling the distance betweem 2 flanks of the ripples
% phi_start_target : starting phase of the target
% on the rate axis (must be odd numbers to match the target parameters)
% coef_noise_generated_m: coeffient matrix of the noise enveloppe generated offline
% coef_smoothing_target:  % if == 0 => normal ripples; if >0, it corresponds to the smoothing coefficient for the gaussian window applied to the disc centered in the image (in units of radius) to make the ripple shape
% coef_smoothing_noise: same thing for the orientation noise
% Targetlevel_dB: level of the target (dB)
% shape_target, shape_noise: shape of the window where the target is in - disc, square or squircle
% angle_degree_surround: orientation of the surround (in degrees)
% spacing_surround :  distance parameter 
% phi_start_surround : starting phase of the surround
% Surroundlevel_dB : level of the surround (dB)
% plot_01: 0 <=> don't plot / 1 <=> plot the targeted spectrogram (freq on
% a log axis)
%
%
% OUTPUTs:
% targetPlusNoise_wav_v: wav file of the total [target + orientation noise] signal created
% coef_noise_m: matrix of amp coefs for the enveloppe of the noise
% coef_target_m: matrix of amp coefs for the enveloppe of the target
% coef_total_m: matrix of amp coefs for the enveloppe of the target+noise
% These are optional:
% target_wav_v: wav file of the target (without noise)
% Noise_wav_v: wav file of the orientation noise (without target) created
% realSNR_dB: the target/noise SNR (in dB)
%
%%
%  E. Ponsot, last viewed on 07/01/2020
%%

%% Parameters

targetPlusNoise_wav_v=zeros(1,duration_sound*fs);
targetPlusNoise_noSurround_wav_v=zeros(1,duration_sound*fs);

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
radius_big=floor(ntones/8); % the radius of the disc is not in argument of the function but it can me manually changed here
radius_surround=floor(ntones/8);

[mask_rescaled_small_v] = fct_toolbox_mask(ntones,duration_sound,fs,radius_small,coef_smoothing_target,shape_target); %  for the target
[mask_rescaled_big_v] = fct_toolbox_mask(ntones,duration_sound,fs,radius_big,coef_smoothing_noise, shape_noise); % for the noise

if strcmp(shape_surround,'squircle')
     mask_rescaled_surround_v=1-mask_rescaled_big_v;
else
    [mask_rescaled_surround_v] = fct_toolbox_mask(ntones,duration_sound,fs,radius_surround,coef_smoothing_noise, shape_surround); % for the noise
    mask_rescaled_surround_v=mask_rescaled_surround_v-mask_rescaled_big_v;
end

size_edge= ntones/64; % you have to change here if you want
radius_smoothing=ntones/16; % you have to change here if you want
smooth_edges_m=fct_toolbox_smooth(ntones, duration_sound, fs, size_edge, radius_smoothing, coef_smoothing_noise);
mask_rescaled_surround_v=mask_rescaled_surround_v.*smooth_edges_m;

%% Loop
% to make a rotation of a ripple
coef_noise_m=zeros(ntones, duration_sound * fs);
coef_target_m=zeros(ntones, duration_sound * fs);
coef_surround_m=zeros(ntones, duration_sound * fs);

halfSize_x = duration_sound * fs/2;
halfSize_y = ntones/2;
[xx,yy] = meshgrid(linspace(-halfSize_x,halfSize_x,2*halfSize_x), linspace(-halfSize_y,halfSize_y,2*halfSize_y));


% reading the noise
coef_noise_m(:,:) = coef_noise_generated_m(:,:);

% windowing the noise
if coef_smoothing_noise>0
    coef_noise_m(:,:)=  coef_noise_m(:,:).*mask_rescaled_big_v;
end

%create the target
coef_target_m(:,:) = coef_target_m(:,:)  +  10^((Targetlevel_dB)/20).*(  1 + d*sin(2*pi*spacing_target*duration_sound/(ntones)*(xx*cosd(angle_degree_target)*ntones/(duration_sound * fs) + yy*sind(angle_degree_target))+phi_start_target));

% windowing the target
if coef_smoothing_target>0
    coef_target_m(:,:)=  coef_target_m(:,:).*mask_rescaled_small_v;
end

% add noise and target
coef_total_m(:,:) =  coef_noise_m(:,:) + coef_target_m;


% create the surround 
if isnan(Surroundlevel_dB)==0
    
    coef_surround_m(:,:) = 10^((Surroundlevel_dB)/20)*(  1 +  d*sin(2*pi*spacing_surround*duration_sound/(ntones)*(xx*cosd(angle_degree_surround)*ntones/(duration_sound * fs) + yy*sind(angle_degree_surround))+ phi_start_surround)  );
    
    if Surroundlevel_dB == 666 % constant background
        coef_surround_m(:,:) = 10^((Targetlevel_dB)/20)*ones(size(coef_surround_m,1),size(coef_surround_m,2));
    end
    
    if coef_smoothing_target>0
        coef_surround_m(:,:)=  coef_surround_m(:,:).*mask_rescaled_surround_v;
    end
end

% add surround to noise only and to target+noise
coef_noise_plussurround_m(:,:) =  coef_noise_m(:,:) +  coef_surround_m;
coef_total_plussurround_m(:,:) =  coef_total_m(:,:) +  coef_surround_m;


%% Plot
% make the plots of the theoretical spectrograms
if plot_01 ~= 0
    figure
    subplot(2,2,1)
    imagesc(abs(  coef_target_m))
    axis xy
    colorbar
    title('target')
    subplot(2,2,2)
    imagesc(abs(  coef_surround_m))
    axis xy
    colorbar
    title('(surround: optional)')
    subplot(2,2,3)
    imagesc(abs(coef_noise_plussurround_m))
    axis xy
    colorbar
    %caxis([0 30])
    title('noise')
    subplot(2,2,4)
    imagesc(abs(coef_total_plussurround_m))
    axis xy
    colorbar
    title('target +noise + surround')
    %caxis([0 30])
end
colormap('gray')


%% apply the overall matrix of coefficients (the enveloppe) to the carrier
coef_ii=zeros(1,duration_sound * fs);
for ii=1:ntones
    coef_ii(:)=coef_total_plussurround_m(ii,:);
    s_ii=gamma_v(ii)*1./sqrt(f_v(ii))*sin(2*pi*f_v(ii)*t_v+phi_v(ii));
    targetPlusNoise_wav_v=targetPlusNoise_wav_v+coef_ii.*s_ii;
end


%% apply the overall matrix of coefficients (the enveloppe WITHOUT THE SURROUND) to the carrier
coef_ii=zeros(1,duration_sound * fs);
for ii=1:ntones
    coef_ii(:)=coef_total_m(ii,:);
    s_ii=gamma_v(ii)*1./sqrt(f_v(ii))*sin(2*pi*f_v(ii)*t_v+phi_v(ii));
    targetPlusNoise_noSurround_wav_v=targetPlusNoise_noSurround_wav_v+coef_ii.*s_ii;
end

%% normalize rms of the targetPlusNoise_v output %%%% without potential spectrotemporal surround 
if isnan(Surroundlevel_dB)==0
targetPlusNoise_wav_v=targetPlusNoise_wav_v./sqrt(sum(abs(targetPlusNoise_noSurround_wav_v).^2)/length(targetPlusNoise_noSurround_wav_v)); % 
else
targetPlusNoise_wav_v=targetPlusNoise_wav_v./sqrt(sum(abs(targetPlusNoise_wav_v).^2)/length(targetPlusNoise_wav_v));
end


%% if you need the wav file of the target and the noise only, and also to calculate the real signal-to-noise ratio
if nargout>5
    target_wav_v=zeros(1,duration_sound*fs);
    Noise_wav_v=zeros(1,duration_sound*fs);
    coef_target_ii=zeros(1,duration_sound * fs);
    coef_noise_ii=zeros(1,duration_sound * fs);
    
    for ii=1:ntones
        % create the noise only
        coef_noise_ii(:)=coef_noise_m(ii,:);
        s_ii=gamma_v(ii)*1./sqrt(f_v(ii))*sin(2*pi*f_v(ii)*t_v+phi_v(ii));
        Noise_wav_v=Noise_wav_v+coef_noise_ii.*s_ii;
        
        % create the target only
        coef_target_ii(:)=coef_target_m(ii,:);
        s_ii=gamma_v(ii)*1./sqrt(f_v(ii))*sin(2*pi*f_v(ii)*t_v+phi_v(ii));
        target_wav_v=target_wav_v+coef_target_ii.*s_ii;
    end
    
    % normalize rms
    rms_noise=sqrt(sum(abs(Noise_wav_v).^2)/length(Noise_wav_v));
    rms_target=sqrt(sum(abs(target_wav_v).^2)/length(target_wav_v));
    Noise_wav_v=Noise_wav_v./rms_noise;
    target_wav_v=target_wav_v./rms_target;
    realSNR_dB=20*log10(rms_target/rms_noise);
end


end
