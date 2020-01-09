function [newsignal_v] = fct_toolbox_fade_onoff(signal_v,duration_n,fs)

% To  smooth a mono signal onset and offset with a linear ramp on amplitude
% INPUTS: signal_v, duration_n (duration of the smoothing ramp in seconds),
% fs (in Hz)
% OUTPUT: new signal
%%
% E. Ponsot, 2012
%%
nbpoints=floor(duration_n*fs);
onset_profile=linspace(0,1,nbpoints);
offset_profile= 1 - onset_profile;

if size(signal_v,1)<size(signal_v,2)
    newsignal_v=signal_v;
else
    newsignal_v=signal_v';
end

newsignal_v(1:nbpoints)=newsignal_v(1:nbpoints).*onset_profile;
newsignal_v(end-nbpoints+1:end)=newsignal_v(end-nbpoints+1:end).*offset_profile;