function [pos1, pos2] = quadsolve(satt, times)

    c = 299792.458;
    u1 = [-2*satt(1,1) + 2*satt(2,1); -2*satt(1,1) + 2*satt(3,1); -2*satt(1,1) + 2*satt(4,1)];
    u2 = [-2*satt(1,2) + 2*satt(2,2); -2*satt(1,2) + 2*satt(3,2); -2*satt(1,2) + 2*satt(4,2)];
    u3 = [-2*satt(1,3) + 2*satt(2,3); -2*satt(1,3) + 2*satt(3,3); -2*satt(1,3) + 2*satt(4,3)];
    u4 = [2*(c^2)*times(1) - 2*(c^2)*times(2); 2*(c^2)*times(1) - 2*(c^2)*times(3); 2*(c^2)*times(1) - 2*(c^2)*times(4)];
    cr =[-(satt(1,1)^2)+(satt(2,1)^2)-(satt(1,2)^2)+(satt(2,2)^2)-(satt(1,3)^2)+(satt(2,3)^2)+(c^2)*(times(1)^2)-(c^2)*(times(2)^2);-(satt(1,1)^2)+(satt(3,1)^2)-(satt(1,2)^2)+(satt(3,2)^2)-(satt(1,3)^2)+(satt(3,3)^2)+(c^2)*(times(1)^2)-(c^2)*(times(3)^2);-(satt(1,1)^2)+(satt(4,1)^2)-(satt(1,2)^2)+(satt(4,2)^2)-(satt(1,3)^2)+(satt(4,3)^2)+(c^2)*(times(1)^2)-(c^2)*(times(4)^2)];
    a = [u1,u2,u3,u4,cr];
    R = rref(a);
    r14 = R(1,4);
    r15 = R(1,5);
    r24 = R(2,4);
    r25 = R(2,5);
    r34 = R(3,4);
    r35 = R(3,5);
    
    coeff1 = r14^2 +r24^2 +r34^2 -c^2;
    coeff2 = 2*satt(1,1)*r14 -2*r14*r15 +2*satt(1,2)*r24 -2*r24*r25 +2*satt(1,3)*r34 -2*r34*r35 +2*c^2*times(1); 
    coeff3 = satt(1,1)^2 +r15^2 -2*satt(1,1)*r15 +satt(1,2)^2 +r25^2 -2*satt(1,2)*r25 +satt(1,3)^2 +r35^2 -2*satt(1,3)*r35 -c^2*times(1)^2;
    
    d1 = (-coeff2 + sqrt(coeff2^2 - 4*coeff1*coeff3)) / (2*coeff1);
    d2 = (-coeff2 - sqrt(coeff2^2 - 4*coeff1*coeff3)) / (2*coeff1);
    
    pos1 = [-r14 * d1 + r15; -r24 * d1 + r25; -r34 * d1 + r35; d1];
    pos2 = [-r14 * d2 + r15; -r24 * d2 + r25; -r34 * d2 + r35; d2];