

m1 = 0;
m2s = [0:.1:5];

s1 = 1;
s2 = 1;

N = length(m2s);
ov_min = zeros(1,N);
ov_tim = zeros(1,N);
ov_tim2 = zeros(1,N);

f1 = @(x) gaussian(x, m1, s1);
for mi = 1:length(m2s)
    m2 = m2s(mi);
    ov_min(mi) = gaussiansOverlap(m1, s1, m2, s2);

    f2 = @(x) gaussian(x, m2, s2);
    q_calc = quad(@(x) f1(x).*f2(x), -10, 15);
    ov_tim(mi) = q_calc;
    
    A = 1/(2*pi*s1*s2);
    a = 1/(2*s1^2) + 1/(2*s2^2);
    b = m1/(s1^2) + m2/(s2^2);
    c = -m1^2/(2*s1^2) - m2^2/(2*s2^2);
    q_an = A*sqrt(pi/a)*exp(b^2/(4*a)+c);
    ov_tim2(mi) = q_an;
end

figure(1); clf;
plot(m2s, ov_min, 'bo-', m2s, ov_tim, 'go-', m2s, ov_tim2, 'r.:');

figure(2)
plot(ov_min, ov_tim, 'r.');
