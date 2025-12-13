# mumt307-finalproject-convolution-reverb
Author: Andrew Behling\
Date Created: 2025-12-12

A convolution reverb in MATLAB for the MUMT 307 final project.

## Overview
  This repo contains the code that I wrote for the final project in MUMT 307. It contains a folder called tests with the applicable files both written and taken from sources cited in the references section. The `main.m` file contains a main loop that calls pertinent functions for determining the impulse response of a sample on frequency bands, stretches the impulse response to a desired length, and produces a convolved signal that can be mixed according to user specifications. Each distinct section apart from initial audio file loading and user specification uses a MATLAB function contained in a separate file. References are listed in each file where material is taken from external sources as well as at the end of this markdown. This markdown will serve to offer insights to the workings of each function as well as a number of issues faced and improvements. The audio files used in testing are all in a directory called audio. The AIR_1_4 directory contains a number of open source impulse responses available at the [AIR Database](https://www.iks.rwth-aachen.de/en/research/tools-downloads/databases/aachen-impulse-response-database/). The rest of the audio was either recorded or synthesized by me.
  
## Code breakdown
The `main.m` program begins with user specifications as follows:
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
The adjustable parameters are the audio block size (`blockSize`), desired T60 response time for reverberation (`T60_3b`), the desired frequency bands (`freqBands`), and an arbitrary start delay to remove the excitation of the impulse response.

Next, the desired audio files are imported and processed based on a few edge cases: mono impulse response, and sample versus impulse size. Please see the `main.m` for details on how this was completed. Next the `getImpulseT60()` function is called with user defined parameters: impulse file, dB threshold, frequency bands, and sample rate. The dB threshold can be modified as recordings made with less capable microphones or in noisier environments can cause a 60 dB to be impossible. The function works as follows:
* The impulse signal is blocked and windowed with a flat top window for amplitude accuracy.
* For each block, an FFT is taken with zero padding.
```
 YT60 = 20*log10(abs(fft([ impulse(:, iStart+1:iStart+blockSizeT60), zeros(2,zeroPadT60)], [], 2)));
```
* In each frequency band, the peaks are found using `findpeaks.m`, and the maximum value is stored in a previously preallocated matrix.
* In the case that a peak was not found in that frequency band on the specific block, the peaks for that block are set to the initial values.
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
* Because the signal is zero padded, the resulting signal must be interpolated in the frequency domain to have the proper number of samples for the given window size.
* The DFT is taken of the interpolated output signal.
```
outWet(1, iStart+1:iStart+blockSize) = outWet(1, iStart+1:iStart+blockSize) + interp1(iStart+1:iStart+blockSize+zeroPad, real(ifft(YOut(1,:))), iStart+1:iStart+blockSize);
outWet(2, iStart+1:iStart+blockSize) = outWet(2, iStart+1:iStart+blockSize) + interp1(iStart+1:iStart+blockSize+zeroPad, real(ifft(YOut(2,:))), iStart+1:iStart+blockSize);
```
* The MATLAB functions `tic` and `toc` are used to display the time elapsed during the processing of the convolution for an entire sample.
The resulting output signal is finally mixed with a dry signal according to a user specified mix value.
```
mix = .3;
mixWet = mix;
mixDry = 1 - mixWet;

NOut = size(outWet,2);
out = zeros(2,NOut);

for i = 1:NOut

    out(:,i) = mixDry * in(:,i) + mixWet * outWet(:,i);

end

soundsc(out, fs); 
```
## Challenges encountered
  I encountered a number of difficulties during this project. Namely, the `getImpulseT60()` was particularly difficult to implement, the convolved signal does not sound great, and overall code malpractice made development slow.

#### `getImpulseT60()`
  The majority of time spent developing this code was spent trying to get the T60 response for each frequency band. I have a suspicion a great deal of the difficulty came from a suboptimal approach to finding the T60 time. As opposed to finding a T60 value in the time domain, the introduction of frequency bands necessitated analysis in the frequency domain. The `findpeaks()` was used to track the evolution of spectral peaks in the frequency domain on a given frequency band. The use of `findpeaks()` was inspired by homework 6, but instead of an unwindowed approach, a flattop window was chosen for the best amplitude accuracy. 
  
  Right away, I ran into issues with its implantation as it was rather difficult to accurately track peaks in the lower frequency band. Initially, without zero padding the FFT, there was not nearly high enough spectral resolution to produce a peak value in the first frequency band (0-100 Hz). The frequency resolution was roughly 86 Hz, and this was not nearly low enough to track peaks in the lower band. There were two viable solutions to problem: increase the value of the low-mid frequency band or zero pad the FFT. I decided that 100 Hz was a viable low-mid band separator; thus, the windows were zero padded with an arbitrary value of 2048 samples. With a block size of 512 samples, the total frequency resolution was ~17 Hz. This allowed for enough data points in the lower frequency band to see distinct peaks. 
  
  The low frequency band was quite trivial in terms of a peak tracking because there was typically only one peak value per window. The `findpeaks()` method provided trouble for the mid band, however, as the function was not selecting the maximum peak in the band. As was done in homework 6, I initially tried to find a limited number of peaks in the frequency band for a given audio block, and to find the true maximum peak, at the expense of computation time, I opted to find all peaks in a frequency band, and the maximum value was taken as the band's maximum. I unfortunately do not know a less computationally expensive method to find the T60 time, but I think the method shown above is capable of finding the band T60 responses for a properly recorded impulse response.

#### Final convolved signal
  A number of result files are presented in the directory `results`. The fully wet output signal sounds almost bitcrushed, so I have a few suspicions as to why the final result does not sound particularly good. First, the frequency bands were segmented without filtering and mixing, so it is likely high frequency content is created for a given window where the value of the convolution at each band does not match. It would have been a more prudent approach to use the frequency band as the centre frequency for a highpass and lowpass filter at each band interface. The rolloff could be computed to amply scale the values of the output signal at the interfaces. One problem with this approach is the impulse response would have to be calculated 6 times with most of, if not, the entire frequency spectrum, band filtered, and then summed. This would increase the computational cost, but the quality would increase significantly most likely. 
  
  Another problem with the convolution approach taken was the interpolation of the FFT before taking the DFT. As mentioned in the code overview, with a zero padded FFT, the DFT returns a time domain window of size `blockSize + zeroPad`. Thus, it is impossible to fit the DFT into the original size window. To combat this, the values of each channel of the output convolved were interpolated at frequency intervals that correspond with the block size. This method likely caused the creation of further noise as linear interpolation in the frequency domain is likely to create erroneous peak values. Finally, the block size could have been chosen smarter as the FFT was typically completed with a block size of 64 samples. This allows for use as a real time audio plugin, but it is not particularly conducive to good processing for non real time signals. 
  
  Finally, a question that has made the development of this code a bit tough comes in the actual application of the impulse response in the frequency domain: Should the impulse be triggered every input audio block. In a true reverberating space, the perceived reverb is a combination of many excitations from the entire source signal. For example, for a sustained or dynamic signal, the listener hears an impulse response for the entire signal and not just the initial transient sound. When trying to emulate a reverberating space, an initial approach can be taken to discretize the entire input signal. Each individual sample with a non-zero value therefore excites the room and necessitates an impulse response. However, this is not computationally viable nor necessary given the limits of human hearing. Continuing with this logic and decreasing the discretization to nominal audio block sizes, a new impulse signal should be triggered every audio block. This presents some challenges with the code architecture, however.
  
  As the impulse signal needs to be retriggered every audio block, there are multiple impulse signals that would need to be processed in parallel. For an impulse signal that is N samples long, there would be `ceil(N/blockSize)` impulse signals being convolved during ``steady state" operation of the reverb. For lower block sizes and zero padding, this seems incredibly taxing computationally. Given my code already produces so many unwanted high frequencies, adding a steady state response aspect would have likely been a noise machine. Therefore, in its current state, this convolution reverb is best used for plucked and decaying sounds.

#### Code issues
  Generally speaking, there were a number of delays caused by inexperience with MATLAB and poor programming. I initially began this project in Python as it is the language I have dealt with the most in my studies and personal endeavors; however, Python's real time audio and processing packages were not well documented enough to make quick enough progress on this project. Thus, I switched to MATLAB, a language that I am familiar with but inexperienced. The most egregious code sins I committed were in the form of array structuring and indexing. To represent the impulse response of the 6 band and channel combinations, I opted to use a 3D matrix. The shape of this matrix was incredibly difficult to work with as the audio files I was manipulating had time data represented on the N axis of a 2 x N array. The impulse signals were calculated and saved in a N x 3 x 2 array. This caused many errors to throw as I was trying to build the convolution code.

  Second, inexperience with MATLAB caused me to spend quite some time in the documentation and stack overflow looking for syntax and functions to do what I needed. The most pertinent example of this came again from the three-band impulse matrix. Trying to manipulate both channels of a single band caused MATLAB to throw an error pertaining to an attempt to superimpose a 3D matrix with a 0-length axis into a 2D array shape. Thus, I had to look for the squeeze function to remove the 0-length dimension, and this took longer than expected to find. Overall, the code is likely longer and more complex than it needs to be as a result of my inexperience with coding and MATLAB.

## Future improvements
  In no particular order, some of the changes I would like to implement in the future (as a passion project) is the development of a real time audio VST plugin for use in a DAW program. To be able to compute the FFT and convolve the signals quicker, it will require a more intelligent application of the impulse signal processing. Finally, I would like to add a more intuitive and tangible way to change the reverberation time to the user's desires. Currently, it is not particularly clear what changing the reverberation time actually corresponds to, but I think I am hindered by my preconceived notions as to what a reverb is supposed to sound like.

## References
[1] [MUMT 307 Course Notes, Prof. Gary Scavone, McGill University, 2025](https://caml.music.mcgill.ca/~gary/307/index.html)\
[2] [Phase Vocoder, Prof. Dan Ellis, Columbia University, 2003](https://www.ee.columbia.edu/~dpwe/LabROSA/matlab/pvoc/)\
[3] [MATLAB Documentation, MathWorks, 2025](https://www.mathworks.com/help/matlab/index.html)\
<img width="468" height="645" alt="image" src="https://github.com/user-attachments/assets/879a69b0-9f6f-48a3-888a-5d6daf8e9e5a" />
