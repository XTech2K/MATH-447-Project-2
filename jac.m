function [Df, F] = jac(pos, satt, times)

    c = 299792.458;
    
    len = length(times);
    F = zeros(4, 1);
    Df = zeros(4, 4);

    for i = 1:len
        
        temp = @(j) (pos(j) - satt(i, j)) ^ 2;
        sphere = sqrt(temp(1) + temp(2) + temp(3));
        
        F(i) = sphere - c * (times(i) - pos(4));
        
        for j = 1:3
            
            Df(i, j) = (2 * pos(j) - 2 * satt(i, j)) ...
                / (2 * sphere);
            
        end
        
        Df(i, 4) = c;
        
    end

end
