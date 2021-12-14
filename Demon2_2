syms s;
% Второй пункт бля



w = zeros(1, 6);
disp(['Введите точки абсциссы точек перегиба полученного графика слева' ...
    ' направо (их должно быть 6):']);
for i = 1:6
    w(i) = 1 / input("Точка " + i + ": ");
end

W_ky_1 = makeTf(w(2), w(1)); W_ky_2 = makeTf(w(3), w(5)); 
W_ky_3 = makeTf(w(4), w(6));

% Шиша
% Wky = (0.035 * s + 1)(0.02777 * s + 1)(0.01695 * s + 1) / +...
% (31.361 * s + 1)(0.0008615 * s + 1)(0.0004498 * s + 1);
Ws = 32 / ((0.00035 * s^2 + 0.0405 * s + 1));

% Журов полученный нами
% Wky = (0.0476*s + 1)(0.0256*s + 1)(0.0092*s + 1)/+...
%     (33.3*s + 1)(0.001176*s + 1)(0.000133*s + 1)

% Журов из его файла
% Ws = (562.5 * 8750) / ((0.026 * s + 1) * (0.0055 * s + 1) * s);
% Wky = (0.033 * s + 1)*(0.026 * s + 1)*(0.0055 * s + 1) /+...
%     ((21.8 * s + 1)*(0.00042 * s + 1)*(0.00042 * s + 1));
% W_ky_1 = (0.033 * s + 1) / (21.8 * s + 1);
% W_ky_2 = (0.026 * s + 1) / (0.00042 * s + 1);
% W_ky_3 = (0.0055 * s + 1) / (0.00042 * s + 1);

% Серба
% Wky = (0.04*s + 1)(0.038*s + 1)(0.0095*s + 1)/+...
%     (16.67*s + 1)(0.00092*s + 1)(0.00092*s + 1)

[num, den] = numden_crutch(W_ky_1);
C1 = capacitor_from_user(1);

% расчёт фильтра 1
R2 = vpa(num / C1, 4);
R1 = vpa(den / C1 - R2, 4);
disp("С1 = ");
disp(C1);
disp("R1 = ");
R1 = parseE24(R1);
disp(R1)
disp("R2 = ");
R2 = parseE24(R2);
disp(R2);

[w1, w2] = frequency_remake(R2 * C1, (R1 + R2) * C1, 1);

W_ky_remake_1 = vpa((1 / w1 * s + 1) / (1 / w2 * s + 1), 5);
disp("В итоге передаточная функция фильтра 1 будет равна");
display(W_ky_remake_1);

% расчёт фильтра 2
[num, den] = numden_crutch(W_ky_2);
C2 = capacitor_from_user(2);

k2 = vpa(den / num, 4);
R3 = vpa(num / C2, 4);
R4 = vpa(k2 * R3 / (1 - k2), 4);
disp("С2 = ");
disp(C2);
disp("R3 = ");
R3 = parseE24(R3);
disp(R3);
disp("R4 = ");
R4 = parseE24(R4);
disp(R4);

[w3, w4] = frequency_remake(R3 * C2, k2 * R3 * C2, 3);
W_ky_remake_2 = vpa((1 / w3 * s + 1) / (1 / w4 * s + 1), 5);
disp("В итоге передаточная функция фильтра 2 будет равна");
display(W_ky_remake_2);

% расчет фильтра 3
[num, den] = numden_crutch(W_ky_3);
C3 = capacitor_from_user(3);

k3 = vpa(den / num, 4);
R5 = vpa(num / C3, 4);
R6 = vpa(k3 * R5 / (1 - k3), 4);
disp("С3 = ");
disp(C3);
disp("R5 = ");
R5 = parseE24(R5);
disp(R6);
disp("R6 = ");
R6 = parseE24(R6);
disp(R6);

[w5, w6] = frequency_remake(R5 * C3, k3 * R5 * C3, 5);
W_ky_remake_3 = vpa((1 / w5 * s + 1) / (1 / w6 * s + 1), 5);
disp("В итоге передаточная функция фильтра 3 будет равна");
display(W_ky_remake_3);

% Пункт 12 Метро Измайловская Работаем
W_ky = vpa(W_ky_remake_1 * W_ky_remake_2 * W_ky_remake_3, 5);
W = vpa(Ws * W_ky, 5);

[Ws_numArr, Ws_denArr] = before_tf(Ws, 0);
Ws = tf(Ws_numArr, Ws_denArr);
[W_ky_numArr, W_ky_denArr] = before_tf(W_ky, 0);
W_ky = tf(W_ky_numArr, W_ky_denArr);
[numArr, denArr] = before_tf(W, 0);
W = tf(numArr, denArr);

disp("Итоговя передаточная функция коректирующего устройсва ");
display(W_ky);

disp("Передаточная функция системы с коректирующем устройсвом ");
display(W);

W_transient = W / (1 + W);
Ws_transient = Ws / (1 + Ws);
step(Ws_transient);
grid on;
hold on;
step(W_transient);
title('Переходный процесс');
legend('Без коректирущего устройсва', 'С коректирующем устройсвом');


% флаг указывает на наличие константы в знаменателе 
% если флаг 1 то знаменатель a1 * s^n + a2 * s^(n-1) + ... + an * s
% если флаг 0 то знаменатель a1 * s^n + a2 * s^(n-1) + ... + an * s + const
function [numArr, denArr] = before_tf(W, flag)
    [n, d] = numden(W);
    num = coeffs(n);
    den = coeffs(d);
    numArr = zeros(1, max(size(num)));
    if (flag == 1)
        denArr = zeros(1, max(size(den)) + 1);
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


function W = makeTf(wUpper, wLower)
    syms s
    W = (wUpper * s + 1) / (wLower * s + 1);
end

function Rres = parseE24(R)
    E24 = [1, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2, 2.2, 2.4, 2.7, 3, 3.3, + ...
        3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1];
    tt = 10 ^ floor(log10(R));
    R = R / tt;
    res = 0;
    for i = 1:max(size(E24))
        if i == max(size(E24))
            res = E24(i);
            break;
        end
        if E24(i) > R
            diff = abs(R - E24(i - 1)) - abs(R - E24(i));
            if (diff >= 0)
                res = E24(i);
            else
                res = E24(i - 1);
            end
            break;
        end
    end
    Rres = res * tt;
end

function [num, den] = numden_crutch(W_ky)
    [n, d] = numden(W_ky);
    n = coeffs(n);
    d = coeffs(d);
    num = n(2) / n(1);
    den = d(2) / d(1);
end

function C = capacitor_from_user(i)
    disp("Введите значения ёмкости (любое стандартное) в ФАРАДАХ");
    C = input("C" + i + ": ");
end

function [w1Res, w2Res] = frequency_remake(t, T, i)
    w1 = vpa(1 / t, 4);
    w2 = vpa(1 / T, 4);
    disp("Уточнённые значения частот сопряженрия");
    disp("w" + i + " = ");
    disp(w1);
    disp("w" + (i + 1) + " = ");
    disp(w2);
    w1Res = w1;
    w2Res = w2;
end
