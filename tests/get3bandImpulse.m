function impulse_3b = get3bandImpulse(impulse, T60, T60_3b, buffer, fs)

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

end