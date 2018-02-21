function pos = multinewt(satt, times, steps)

    pos = [0 0 6370 0]';

    for n = 1:steps
        
        [Df, F] = jac(pos, satt, times);
        
        change = Df \ -F;
        pos = pos + change;
        
    end

end
