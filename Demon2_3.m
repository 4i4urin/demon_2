% 13 пункт какого-то непонятного говна
syms s;
syms x;
syms w;

% Шиша
Nzad = 19200;
sigma = 0.07;
Emax = 1.2;
w_mid = 36;
Kg = 62500;
Ws = 91428.5714286 / (s^3 + (810*s^2)/7 + (20000)*s/7);
W_ky = (1.75e-05 * s^3 + 0.002106 * s^2 + 0.081 * s + 1) / (1.161e-05 * s^3 + 0.02897 * s^2 + 18.04 * s + 1);
W = (0.0005599 * s^3 + 0.06739 * s^2 + 2.592 * s + 32) / (4.063e-09 * s^5 + 1.061e-05 * s^4 + 0.007498 * s^3 + 0.7598 * s^2 + 18.08 * s + 1);

% Журов
% sigma = 7;
% Tmax = 0.038;
% Emax = 1.4;
% Kg = 62500;
% w_mid = 165;
% Nzad = 8000;
% Ws = (56.5*8750) / ((0.026 * s + 1) * (0.0055 * s + 1) * s);



disp("По теореме Котельникова определим период квантования");
T = pi / w_mid;
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
disp(vpa(Ws, 5));


disp("Проведём z-преобразование:");
[numArr, denArr] = before_tf(Wnch, 2);
Wnch_tf = tf(numArr, denArr);
display(Wnch_tf);

z_trans = c2d(Wnch_tf, T, 'impulse') * 1/T^2;
display(zpk(z_trans));


[n, d] = tfdata(z_trans);

n = round(cell2mat(n), 9); d = cell2mat(d);
n = circshift(n, 1);
d = divide_by_z_1(d);
Wnch_ztrans = tf(n, d);
display(zpk(Wnch_ztrans));

[num, den] = tfdata(Wnch_ztrans);
num = cell2mat(num); den = cell2mat(den);

replace = vpa(poly2sym(num, w) / poly2sym(den, w), 5);
replace = subs(replace, w, (1 + w) / (1 - w));
[rep_num, rep_den] = numden(replace);
replace_tf = zpk(tf(sym2poly(rep_num), sym2poly(rep_den)));
display(replace_tf);

[n, d] = tfdata(replace_tf);
n = round(cell2mat(n), 5); d = round(cell2mat(d), 5);
show_replace = vpa(poly2sym(n, w)/poly2sym(d, w), 5);
display(show_replace);

[num, den] = numden(show_replace);
num = vpa(factor(num, 'FactorMode', 'real'), 5);
den = vpa(factor(den, 'FactorMode', 'real'), 5);
disp(num); disp(den);

% костыль
% построение неизменяемого графика
vc_arr_den = zeros(2, 0);
vc_arr_num = zeros(3, 0);
for i = 1:max(size(vc_arr_den))
    vc_arr_den(i) = vpa(min(abs(coeffs(den(2 + i)))), 5);
end
for i = 1:max(size(vc_arr_num))
    vc_arr_num(i) = vpa(max(abs(coeffs(num(1 + i)))), 5);
end

vc_arr_den = sort(vc_arr_den);
vc_arr_num = sort(vc_arr_num);
%костыль
vc1 = vc_arr_den(1);
vc2 = vc_arr_den(2);
vc3 = vc_arr_num(1);
vc4 = vc_arr_num(2);
vc5 = vc_arr_num(3);

disp("vc1 = " );
disp(vc1);
disp("vc2 = " );
disp(vc2);
disp("vc3 = " );
disp(vc3);
disp("vc4 = " );
disp(vc4);
disp("vc5 = " );
disp(vc5);
k = abs(num(1) * vc4 * vc5 / (den(1) * vc1 * vc2));
disp(k);
disp("20 * log(k) = ");
k = 20* log10(k);
disp(k);

% наверно надо чтобы пользователь вводил нижнюю и верхнюю границу гррафика
x_begin = 0.000001;
x1 = vc1;
x2 = vc2;
x3 = vc3; % = 1;
x4 = vc4;
x5 = vc5;
x_end = 100;

y3 = k;
y2 = y3 + 60 * abs(log10(x3 / x2));
y1 = y2 + 40 * abs(log10(x2 / x1));
y_begin = y1 + 20 * abs(log10(x1 / x_begin));

y4 = y3 - 40 * abs(log10(x4 / x3));
y5 = y4 - 20 * abs(log10(x5 / x4));
y_end = y5;

semilogx([x_begin, x1, x2, x3, x4, x5, x_end], [y_begin, y1, y2, y3, y4, y5, y_end]);
grid on;
hold on;

% построение желаемого графика
disp("Рассчитаем параметры желаемой ЛАЧХ");
disp(w_mid);
v_mid = tan(w_mid * T / 2);
disp(v_mid);
v_low = 0.16 * v_mid;
v_high = 6.5 * v_mid;

y_mid = 0;
y_low = y_mid + 20 * abs(log10(v_low) - log10(v_mid));
y_begin = y_low + 40 * abs(log10(x_begin) - log10(v_low));
y_high = y_mid - 20 * abs(log10(v_mid) - log10(v_high));
y_end = y_high;

semilogx([x_begin v_low, v_mid, v_high, x_end], [y_begin, y_low, y_mid, y_high, y_end], 'g');
hold on;

disp("частоты перегиба");
disp("w1 = ");
disp(v_low);
disp("w2 = ");
disp(vc1);
disp("w3 = ");
disp(vc2);
disp("w4 = ");
disp(v_high);
disp("w5 = ");
disp(vc5);
disp("w6 = ");
disp(vc3);
disp("w7 = ");
disp(vc4);
disp("w8 = ");
disp("Возможно образование дополнительных точек перегиба смотри график");


% флаг указывает на наличие константы в знаменателе 
% если флаг 1 то знаменатель a1 * s^n + a2 * s^(n-1) + ... + an * s^flag
% если флаг 0 то знаменатель a1 * s^n + a2 * s^(n-1) + ... + an * s + const
function [numArr, denArr] = before_tf(W, flag)
    [n, d] = numden(W);
    num = coeffs(n);
    den = coeffs(d);
    numArr = zeros(1, max(size(num)));
    if (flag > 0)
        denArr = zeros(1, max(size(den)) + flag);
    else
        denArr = zeros(1, max(size(den)));
    end
    koef = den(1);
    for i = 1:max(size(num))
        numArr(i) = num(max(size(num)) - i + 1) / koef;
    end
    for i = 1:(max(size(den)))
        denArr(i) = den(max(size(den)) - i + 1) / koef;
    end
    % подставляем в последний член 0 чтобы избавится от константы в
    % знаменателе
    %denArr(max(size(den)) + 1) = 0;
end

function den_new = divide_by_z_1(den_old)
    den_new = zeros(max(size(den_old)) - 1, 0);
    den_new(1) = den_old(1);
    for i = 2:max(size(den_old) - 1)
        den_new(i) = den_new(i - 1) + den_old(i);       
    end
end
