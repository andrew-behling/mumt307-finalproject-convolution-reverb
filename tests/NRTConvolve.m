function outWet = NRTConvolve(in, impulse_3b, freqBands, blockSize, ...
                              zeroPad, fs)

    tic

    N = size(in, 2);

    f_ml = freqBands(1);
    f_mh = freqBands(2);
    f_h = freqBands(3);
    
    % preallocating space
    outWet = zeros(2,N);
    
    % let's use a hanning window in this example
    % portions taken from [1].
    
    olap = 50;
    
    hop = floor( blockSize * (100-olap) / 100 );
    
    % noWindows = 2 * floor(NBuffer / blockSize) - 1; % for a 50% overlap
    noWindows = floor(N / blockSize);
    w = window(@hanning, blockSize);
    
    yIn = zeros(2,N);
    
    freq = fs/(blockSize+zeroPad) * (0:blockSize+zeroPad-1);
    
    [~, i_ml] = min(abs(freq - f_ml));
    [~, i_mh] = min(abs(freq - f_mh));
    [~, i_h] = min(abs(freq - f_h));
    
    for n = 1:noWindows
        
        iStart = (n-1) * hop;
    
        YOut = zeros(2,blockSize+zeroPad);
        
        l_yImpulse = squeeze(impulse_3b(iStart+1:iStart+blockSize,1,:));
        m_yImpulse = squeeze(impulse_3b(iStart+1:iStart+blockSize,2,:));
        h_yImpulse = squeeze(impulse_3b(iStart+1:iStart+blockSize,3,:));
    
        l_yImpulse = transpose(l_yImpulse);
        m_yImpulse = transpose(m_yImpulse);
        h_yImpulse = transpose(h_yImpulse);
    
        yIn = in(:,iStart+1:iStart+blockSize);
        
        % taking the STFT
        l_YImpulse = fft([ l_yImpulse zeros(2,zeroPad)], [], 2);
        m_YImpulse = fft([ m_yImpulse zeros(2,zeroPad)], [], 2);
        h_YImpulse = fft([ h_yImpulse zeros(2,zeroPad)], [], 2);
        YIn = fft([yIn zeros(2,zeroPad)], [], 2);
        % point-wise multiplication
    
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
    
        outWet(1, iStart+1:iStart+blockSize) = ...
            outWet(1, iStart+1:iStart+blockSize) + ...
            interp1(iStart+1:iStart+blockSize+zeroPad, ...
            real(ifft(YOut(1,:))), ...
            iStart+1:iStart+blockSize);
    
        outWet(2, iStart+1:iStart+blockSize) = ...
            outWet(2, iStart+1:iStart+blockSize) + ...
            interp1(iStart+1:iStart+blockSize+zeroPad, ...
            real(ifft(YOut(2,:))), ... 
            iStart+1:iStart+blockSize);
    
    end

    toc

end
