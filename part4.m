p = 26570;

rho = [1.0076909362467;
       0.5254435044599;
       1.3717504919003;
       1.3441600704864];

theta = [1.9925858282077;
         2.1951440493367;
         4.1171174116978;
         4.7144605351959];

c = 299792.458;
[x, y, z, d] = setpos([0 0 6370 0.0001]);

f = @(fr,ft,i) p * fr(rho(i)) * ft(theta(i));

satts = zeros(4, 3);
range = zeros(4, 1);

for i=1:length(rho)
    satts(i,1:3) = [f(@cos,@cos,i) f(@cos,@sin,i) f(@sin,@(x) 1, i)];
    range(i) = sqrt(sum((satts(i, :) - [x y z]) .^ 2));
end

times = arrayfun(@(r) d + r / c, range);

tfactor = 10 ^ -8;

dt = [-tfactor -tfactor tfactor -tfactor]';

[x1, y1, z1, d1] = setpos(multinewt(satts, times + dt, 100));

change = abs([x y z] - [x1 y1 z1]);

EMF = norm(change, Inf) / (c * norm(dt, Inf));