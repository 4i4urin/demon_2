% 13 пункт какого-то непонятного говна
syms s;
syms x;
syms w;

B = containers.Map('KeyType', 'double', 'ValueType', 'double');
B(5) = 6.5; B(10) = 6.7; B(20) = 6.9; B(25) = 8.8; B(30) = 11.3; 
B(35) = 14.1; B(40) = 16.9;

k = cell2mat(keys(B));
b = 0;
for i = 1:length(B)
    key = k(i);
    if (key > sigma)
        break;
    end
    b = B(key);
end

% Шиша
% Nzad = 19200;
% sigma = 0.07;
% Emax = 1.2;
% w_mid = 36;
% Kg = 62500;
% Ws = 91428.5714286 / (s^3 + (810*s^2)/7 + (20000)*s/7);
% W_ky = (1.75e-05 * s^3 + 0.002106 * s^2 + 0.081 * s + 1) / (1.161e-05 * s^3 + 0.02897 * s^2 + 18.04 * s + 1);
% W = (0.0005599 * s^3 + 0.06739 * s^2 + 2.592 * s + 32) / (4.063e-09 * s^5 + 1.061e-05 * s^4 + 0.007498 * s^3 + 0.7598 * s^2 + 18.08 * s + 1);

% Журов
% sigma = 7;
% Tmax = 0.038;
% Emax = 1.4;
% Kg = 62500;
% w_mid = 165;
% Nzad = 8000;
% Ws = (56.5*8750) / ((0.026 * s + 1) * (0.0055 * s + 1) * s);

% Ws = input("Введите передаточную функцию разомкнутой системы: ");
% [num, den] = numden(Ws);
% Ws_tf = tf(sym2poly(num), sym2poly(den));

sigma = input("Введите сигму: ");
Tmax = input("Введите Tпмакс: ");
Emax = input("Введите Eмакс: ");

Ws = input("Введите передаточную функцию разомкнутой системы: ");
[num, den] = numden(Ws);
Ws_tf = tf(sym2poly(num), sym2poly(den));

disp("По теореме Котельникова определим период квантования");
w_mid = b / Tmax;
w_high = 6.5 * w_mid;
% disp(w_high);
T = pi / w_high;
disp(T);
disp("Возьмём период квантование не более данного");

disp("Пусть ");
T = round(abs(T - 1/6 * T), 4);
disp("T = ");
% T = 0.001;
disp(T);
N = floor(log2(Nzad * (1 + sigma) / (Emax)) + 2);
disp("Определим число разрядов микропроцессора");
disp("N = ");
disp(N);

disp("Передаточная функция разомкнутой системы с требуемым усилением");
display(Ws_tf);
Wz = c2d(Ws_tf, T);
display(zpk(Wz));

[num, den] = tfdata(Wz);
Wz = poly2sym(cell2mat(num), z) / poly2sym(cell2mat(den), z);
Wz = subs(Wz, z, (1 + w) / (1 - w));
[num, den] = numden(Wz);
Wz = tf(sym2poly(num), sym2poly(den));
display(zpk(Wz));
[num, den] = tfdata(Wz);
num = cell2mat(num); den = cell2mat(den);
k = num(max(size(num))) / den(max(size(den)) - 1);
num = poly2sym(num); den = poly2sym(den);

num = factor(num, 'FactorMode', 'real');
den = factor(den, 'FactorMode', 'real');
% disp(num); disp(den);
vc_arr_num = zeros(max(size(num)) - 1, 0);
vc_arr_den = zeros(max(size(den)) - 1, 0);
for i = 1:(max(size(num)) - 1)
    b = abs(coeffs(num(i + 1)));
    if b(1) < 0.0001 || max(size(b)) < 2 % в случае если скобка x а не x + a1
        b = 0;
    end
    vc_arr_num(i) = b(1);
end
for i = 1:(max(size(den)) - 1)
    b = abs(coeffs(den(i + 1)));
    if b(1) < 0.0001 || max(size(b)) < 2 % в случае если скобка x а не x + a1
        b = 0;
    end
    vc_arr_den(i) = b(1);
end
vc = [vc_arr_den, vc_arr_num];
vc = sort(vc);
vc(vc==0) = []; % удвление 0 из массива
disp("Vc от vc1 до vc5");
disp(vc);
disp(k);
k_high = 20 * log10(k);
disp("20 * log(k) = ");
disp(k_high);
disp("По этим данным надо график посторить");

vc_begin = 1e-6;
% vc1 = 0.1853;
% vc2 = 0.3969;
% vc3 = 1;
% vc4 = 1.4849;
% vc5 = 2.1154;
vc_end = 1e3;

% k_y = 54.4856;

vc_begin_y = k_high + 20 * abs(log10(1 / vc_begin));
vc1_y = vc_begin_y - 20 * abs(log10(vc_begin / vc(1)));
vc2_y = vc1_y - 40 * abs(log10(vc(1) / vc(2)));
vc3_y = vc2_y - 60 * abs(log10(vc(2) / vc(3)));
vc4_y = vc3_y - 40 * abs(log10(vc(3) / vc(4)));
vc5_y = vc4_y - 20 * abs(log10(vc(4) / vc(5)));
vc_end_y = vc5_y;
% vc1

x = [vc_begin, vc(1), vc(2), vc(3), vc(4), vc(5), vc_end];
y = [vc_begin_y, vc1_y, vc2_y, vc3_y, vc4_y, vc5_y, vc_end_y];

semilogx(x, y);
grid on;
hold on;


disp("Праметры желаемой ЛАЧХ");
B = containers.Map('KeyType', 'double', 'ValueType', 'double');
B(5) = 6.5; B(10) = 6.7; B(20) = 6.9; B(25) = 8.8; B(30) = 11.3; 
B(35) = 14.1; B(40) = 16.9;

k = cell2mat(keys(B));
b = 0;
for i = 1:length(B)
    key = k(i);
    if (key > sigma)
        break;
    end
    b = B(key);
end



disp("wсрж = ");
disp(w_mid);
v_mid = tan(w_mid * T / 2);
v_low = 0.16 * v_mid;
v_high = 6.5 * v_mid;
disp("vсрж = ");
disp(v_mid);
disp("vн = ");
disp(v_low);
disp("vв = ");
disp(v_high);


v_mid_y = 0;
v_low_y = v_mid_y + 20 * abs(log10(v_low / v_mid));
v_high_y = v_mid_y - 20 * abs(log10(v_mid / v_high));
vc_begin_y = v_low_y + 40 * abs(log10(v_low / vc_begin));
vc_end_y = v_high_y;

x_l = [vc_begin, v_low, v_mid, v_high, vc_end];
y_l = [vc_begin_y, v_low_y, v_mid_y, v_high_y, vc_end_y];
semilogx(x_l, y_l);
grid on;
