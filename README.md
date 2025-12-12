# mumt307-finalproject-convolution-reverb
Author: Andrew Behling\
Date Created: 2025-12-12

A convolution reverb in MATLAB for the MUMT 307 final project.

  This repo contains the code that I wrote for the final project in MUMT 307. It contains a folder called tests with the applicable files both written and taken from sources cited in the references section. The `main.m` file contains a main loop that calls pertinent functions for determining the impulse response of a sample on frequency bands, stretches the impulse response to a desired length, and produces a convolved signal that can be mixed according to user specifications. Each distinct section apart from initial audio file loading and user specification uses a matlab function contained in a separate file. References are listed in each file where material is taken from external sources as well as at the end of this markdown. This markdown will serve to offer insights to the workings of each function as well as a number of issues faced and improvements. The audio files used in testing are all in a directory called audio. The AIR_1_4 directory contains a number of open source impulse responses available at the [AIR Database](https://www.iks.rwth-aachen.de/en/research/tools-downloads/databases/aachen-impulse-response-database/). The rest of the audio was either recorded or synthesized by me.

The `main.m` program begins with user specfications as follows:
```
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
```
The adjustable parameters are the desired T60 response time for reverberation (`T60_3b`), the desired frequency bands (`freqBands`), and an arbitrary start delay to remove the excitation of the impulse response.

Next, the desired audio files are imported and processed based on a few edge cases: mono impulse response, and sample versus impulse size. Please see the `main.m` for details on how this was completed. Next the `getImpulseT60()` function is called with user defined parameters: impulse file, dB threshold, frequency bands, and sample rate. The dB threshold can be modified as recordings made with less capable microphones or in noisier environments can cause a 60 dB to be impossible. The function works as follows:
* The impulse signal is blocked and windowed with a flat top window for amplitude accuracy.
* For each block, an FFT is taken with zero padding.
```
 YT60 = 20*log10(abs(fft([ impulse(:, iStart+1:iStart+blockSizeT60), zeros(2,zeroPadT60)], [], 2)));
```
* In each frequency band, the peaks are found using `findpeaks.m`, and the maximum value is stored in a previously preallocated matrix.
* In the case that a peak was not found in that frequncy band on the specific block, the peaks for that block are set to the initial values.
```
% low frequency band
[l_Pks_L, ~] = findpeaks(YT60(1,1:i_ml));
[l_Pks_R, ~] = findpeaks(YT60(2,1:i_ml));
    
[l_Pks_L, ~] = max(l_Pks_L);
[l_Pks_R, ~] = max(l_Pks_R);

% mid frequency band
[m_Pks_L, ~] = findpeaks(YT60(1,i_ml+1:i_mh));
[m_Pks_R, ~] = findpeaks(YT60(2,i_ml+1:i_mh));
        
[m_Pks_L, ~] = max(m_Pks_L);
[m_Pks_R, ~] = max(m_Pks_R);
        
% high frequency band
[h_Pks_L, ~] = findpeaks(YT60(1,i_mh+1:i_h));
[h_Pks_R, ~] = findpeaks(YT60(2,i_mh+1:i_h));
        
[h_Pks_L, ~] = max(h_Pks_L);
[h_Pks_R, ~] = max(h_Pks_R);

    
if isempty(l_Pks_L), l_Pks_L = Pks_3b(1,1,1); end
if isempty(m_Pks_L),  m_Pks_L = Pks_3b(1,2,1); end
if isempty(h_Pks_L), h_Pks_L = Pks_3b(1,3,1); end
        
if isempty(l_Pks_R), l_Pks_R = Pks_3b(1,1,2); end
if isempty(m_Pks_R), m_Pks_R = Pks_3b(1,2,2); end
if isempty(h_Pks_R),  h_Pks_R = Pks_3b(1,3,2); end
    
Pks_3b(n,:,1) = [l_Pks_L m_Pks_L h_Pks_L];
Pks_3b(n,:,2) = [l_Pks_R m_Pks_R h_Pks_R];
```
* The peak value matrix is iterated on for each band and channel (6 in total), and for the first value that is below the initial maximum minus, the block number is stored and multiplied by the block size and sample period to get a decay time.
```
try T60(1,1) = find(Pks_3b(:,1,1) < Pks_3b(1,1,1) - threshold, 1) * hop;
catch
    T60(1,1) = N_impulse;
end

% left mid band 'T60' response
try T60(1,2) = find(Pks_3b(:,2,1) < Pks_3b(1,2,1) - threshold, 1) * hop;
catch
    T60(1,2) = N_impulse;
end

% left high band 'T60' repsonse
try T60(1,3) = find(Pks_3b(:,3,1) < Pks_3b(1,3,1) - threshold, 1) * hop;
catch
    T60(1,3) = N_impulse;
end

% right low band 'T60' response
try T60(2,1) = find(Pks_3b(:,1,2) < Pks_3b(1,1,2) - threshold, 1) * hop;
catch
    T60(2,1) = N_impulse;
end

% right mid band 'T60' response
try T60(2,2) = find(Pks_3b(:,2,2) < Pks_3b(1,2,2) - threshold, 1) * hop;
catch
    T60(2,2) = N_impulse;
end

% right high band 'T60' repsonse
try T60(2,3) = find(Pks_3b(:,3,2) < Pks_3b(1,3,2) - threshold, 1) * hop;
catch
    T60(2,3) = N_impulse;
end

T60 = T60 / fs;
```
* In the case the threshold is never met, the exception is handled by setting the decay time to the length of the impulse sample submitted by the user.

Once the 6 response values are created, the impulse is stretched according to the users desired decay time. This is done using a [phase vocoder](https://www.ee.columbia.edu/~dpwe/resources/matlab/pvoc/) from Dan Ellis from Columbia. I opted not to develop a phase vocoder myself as it seemed outside the scope of the project. The function `get3bImpulse()` takes the true and desired T60 values for each band and channel and returns 6 different impulse signals corresponding to the 3 frequency bands and L/R channels. The `pvoc.m` code will not be discussed as Ellis has ample and detailed documentation at the link above.
```
% see pvoc.m from [2]
r = T60 ./ T60_3b;

L_impulse = max(T60_3b) * fs + buffer;

% impulse_3b = impulse_3b(sample,band,channel)
impulse_3b = zeros(L_impulse,3,2);

for i = 1:3

    L_impulse = pvoc(impulse(1,:), r(1,i), 256);
    R_impulse = pvoc(impulse(2,:), r(2,i), 256);

    LL = size(L_impulse, 1);
    LR = size(R_impulse, 1);

    impulse_3b(1:LL,i,1) = L_impulse;
    impulse_3b(1:LR,i,2) = R_impulse;
end
```
The function `NRTConvolve()` is used to convolve the two signals (input and three band impulse). It operates as follows:
* An overlapping hanning windowed STFT is done for the impulse and input signals.
```
iStart = (n-1) * hop;

YOut = zeros(2,blockSize+zeroPad);

l_yImpulse = squeeze(impulse_3b(iStart+1:iStart+blockSize,1,:));
m_yImpulse = squeeze(impulse_3b(iStart+1:iStart+blockSize,2,:));
h_yImpulse = squeeze(impulse_3b(iStart+1:iStart+blockSize,3,:));

l_yImpulse = transpose(l_yImpulse) .* w;
m_yImpulse = transpose(m_yImpulse) .* w;
h_yImpulse = transpose(h_yImpulse) .* w;

yIn = in(:,iStart+1:iStart+blockSize) .* w;

% taking the STFT
l_YImpulse = fft([ l_yImpulse zeros(2,zeroPad)], [], 2);
m_YImpulse = fft([ m_yImpulse zeros(2,zeroPad)], [], 2);
h_YImpulse = fft([ h_yImpulse zeros(2,zeroPad)], [], 2);
YIn = fft([yIn zeros(2,zeroPad)], [], 2);
```
* The signals are multiplied at each sample to complete the convolution.
* The impulse signal is segmented into the three bands, and the input is multiplied by the value of the corresponding impulse signal at a given frequency. 
```
for i = 1:blockSize+zeroPad
    
    if i <= i_ml
        
        YOut(:,i) = YOut(:,i) + l_YImpulse(:,i) .* YIn(:,i);
        continue
   
    end

    if i < i_mh

        YOut(:,i) = YOut(:,i) + m_YImpulse(:,i) .* YIn(:,i);
        continue

    end

    if i < i_h

        YOut(:,i) = YOut(:,i) + h_YImpulse(:,i) .* YIn(:,i);
        continue

    end

    YOut(:,i) = YOut(:,i) + h_YImpulse(:,i) .* YIn(:,i);

end
```
* Because the signal is zero padded, the resulting signal must be interpolated in the frequency domain to have the proper amount of samples for the given window size.
* The DFT is taken of the interpolated output signal.
```
outWet(1, iStart+1:iStart+blockSize) = outWet(1, iStart+1:iStart+blockSize) + interp1(iStart+1:iStart+blockSize+zeroPad, real(ifft(YOut(1,:))), iStart+1:iStart+blockSize);
outWet(2, iStart+1:iStart+blockSize) = outWet(2, iStart+1:iStart+blockSize) + interp1(iStart+1:iStart+blockSize+zeroPad, real(ifft(YOut(2,:))), iStart+1:iStart+blockSize);
```
