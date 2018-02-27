function pos = gaussnewt(satt, times, steps)

    pos = [0 0 6370 0]';

    for n = 1:steps
        
        [Df, F] = jac(pos, satt, times);
        
        change = (Df' * Df) \ -(Df' * F);
        pos = pos + change;
        
    end
    condition = cond(Df)

end
