p = 26570;
num_satellites = 12;
p_theta_satellites = 2 + (0:num_satellites-1)*2;
theta = transpose(2 * pi ./ p_theta_satellites);
rho = transpose(pi/2 ./ p_theta_satellites);
%rho = [0.80755085009036;
 %      1.0627655885102;
 %      0.12720873714295;
 %      0.26831526105679;
 %      1.0367076939655;
 %      0.94580206143474;
 %      0.0034406210684407;
 %      0.32567632959023;
 %      0.31914958408994;
 %      0.5027507132258;
 %      1.0485209698316;
 %      0.74499610670609];

%theta = [2.3950417117565;
%         0.061044328278417;
%         2.8534379686105;
 %        2.6757255145701;
 %        0.36788490020106;
 %        1.7933648084073;
 %        2.7851000410435;
 %        6.2051299541654;
 %        4.2100723722438;
 %        2.3679097645394;
 %        0.95042038155274;
 %        1.7846681747887];

c = 299792.458;
[x, y, z, d] = setpos([0 0 6370 0.0001]);

f = @(fr,ft,i) p * fr(rho(i)) * ft(theta(i));

satts = zeros(4, 3);
range = zeros(4, 1);

for i=1:length(rho)
    satts(i,1:3) = [f(@cos,@cos,i) f(@cos,@sin,i) f(@sin,@(x) 1, i)];
    range(i) = sqrt(sum((satts(i, :) - [x y z]) .^ 2));
end

time = [arrayfun(@(r) d + r / c, range);
times = zeros(length(rho), 4);

for i = 1:4
	a = zeros(length(rho))
	a(i) = a(i) + 10^-8
	times(i) = time + a
end

pos1 = gaussnewt(satts, times(1), 10)
pos2 = gaussnewt(satts, times(2), 10)
pos3 = gaussnewt(satts, times(3), 10)
pos4 = gaussnewt(satts, times(4), 10)