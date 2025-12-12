function T60 = getImpulseT60(impulse, threshold, freqBands, fs) 

    f_ml = freqBands(1);
    f_mh = freqBands(2);
    f_h = freqBands(3);

    blockSizeT60 = 512;
    N_impulse = size(impulse, 2);

    noBlocksT60 = 2 * floor( N_impulse / blockSizeT60 ) - 1;
    
    % preallocating some space
    Pks_3b = zeros(noBlocksT60,3,2);
    T60 = zeros(2,3);
    
    % let's use a flattop window for amplitude accuracy!
    % portions taken from [1].
    
    olap = 50;
    zeroPadT60 = 2048;
    
    hop = floor( blockSizeT60 * (100-olap) / 100 );
    
    freqT60 = fs/(blockSizeT60+zeroPadT60) * (0:blockSizeT60+zeroPadT60-1);
    
    [~, i_ml] = min(abs(freqT60 - f_ml));
    [~, i_mh] = min(abs(freqT60 - f_mh));
    [~, i_h] = min(abs(freqT60 - f_h));
    
    for n = 1:noBlocksT60
    
        iStart = (n-1) * hop;
    
        YT60 = 20*log10(abs(fft([ impulse(:, iStart+1:iStart+blockSizeT60) ...
                     zeros(2,zeroPadT60)], [], 2)));
        
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
    
    end
    
    % recordings aren't good enoough for T60
    
    % Pks = Pks(block,band,channel)
    % left low band T50 response
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

end