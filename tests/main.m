% by andrew behling
% 
% a multiband convolution reverb algorithm main loop using a timestretched 
% impulse signal to modify the reverb length
% 
% references:
% 
% [1] OLAW.m, Gary P. Scavone, McGill University, 2007.
% [2] pvoc.m, D. P. W. Ellis, Columbia University, 2002.

%% definitions
fs = 44100;
ts = 1/fs;

% and audio processing respectively
blockSize = 64;

% desired T60 in s
T60_3b = [5 3 1.5];

% desired frequency bands for reverb
freqBands = [100, 2000, 20000];

% buffer for where the IR should start
startBuffer = 1024;

%% importing an audio file for the input

% import each audio file
[in, fs_in] = audioread('glasslid5.mp3');
% [in, fs_in] = audioread('guitar.wav');

in = transpose(in);

% handinling a mono input
if size(in, 1) == 1

    placeholder = in;
    in = zeros(2, size(placeholder, 2));
    in(1,:) = placeholder;
    in(2,:) = placeholder;

end

%% importing an IR I made
% [impulse, fs_impulse] = audioread('FDA.mp3');
[impulse, fs_impulse] = audioread('creepy.mp3');

%% importing an opensource IR
[impulse_L, fs_impulse] = ...
    audioread('AIR_1_4/AIR_wav_files/air_stairway_1_1_1_165_mls.wav');
[impulse_R, ~] = ...
    audioread('AIR_1_4/AIR_wav_files/air_stairway_1_1_2_165_mls.wav');


impulse_L = transpose(impulse_L);
impulse_R = transpose(impulse_R);

L_impulse = size(impulse_L,2);
impulse = zeros(2, L_impulse);

impulse(1,:) = impulse_L;
impulse(2,:) = impulse_R;


%%

impulse = transpose(impulse);

% resampling if the sample rate is not the same as defined above
if fs_in ~= fs
    
    in = audioresample(in, 'InputRate', fs_in, 'OutputRate', fs);

end

if fs_impulse ~= fs
 
     impulse = audioresample(impulse, 'InputRate', fs_impulse, ...
         'OutputRate', fs);
 
end

% defining lengths of each file
N_in = size(in);
N_in = N_in(2);

% impulse(1,:) = impulse(1,:) ./ max(abs(impulse(1,:)));
% impulse(2,:) = impulse(2,:) ./ max(abs(impulse(2,:)));

N_impulse = size(impulse, 2);

if N_impulse > N_in

    in = [in zeros(2,N_impulse - N_in)];
    N = N_impulse;

else

    impulse = [impulse zeros(2,N_in - N_impulse)];
    N = N_in;

end

%% finding the t60 repsonse of the impulse

% see getImpulseT60.m
T60 = getImpulseT60(impulse, 60, freqBands, fs);

%% stretching the impulses to the desired lengths

% see get3bandImpulse.m
impulse_3b = get3bandImpulse(impulse, T60, T60_3b, 1024, fs);

soundsc(squeeze(impulse_3b(:, 3,:)), fs);

%% convolving the signal
outWet = NRTConvolve(in, impulse_3b, freqBands, blockSize, 4096, fs);


%%
% soundsc(outWet, fs);
mix = 1;
mixWet = mix;
mixDry = 1 - mixWet;

NOut = size(outWet,2);
out = zeros(2,NOut);

for i = 1:NOut

    out(:,i) = mixDry * in(:,i) + mixWet * outWet(:,i);

end

soundsc(out, fs);

