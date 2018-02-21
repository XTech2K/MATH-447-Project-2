function [x,y,z,d] = setpos(pos)

    c = num2cell(pos);
    [x,y,z,d] = c{:};

end
