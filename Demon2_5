% пункт 15 
% Корневой годограф
% Кто бля?
% Годограф

syms s;
syms x;

% Журов
Ws = 562.5 * 8750 / ((0.026 * s + 1) * (0.0055 * s + 1) * s);
Wky_s = (s * 66.667 + 1) * (s * 47.62 + 1) * (s * 10 + 1) * (s * 2 + 1) / +...
    ((s * 3703703.703 + 1) * (s + 1) * (s * 0.597 + 1) * (s * 0.5571 + 1));

% Чел
% Ws = 528.125 * 9375 / (s * (0.011 * s + 1) * (0.022 * s + 1));

disp("Передаточная функция разомкнутой системы W(s) при t = 0:");
disp(Ws);
[num_coppy, trash] = numden(Ws);
trash = sym2poly(trash);
num_coppy = num_coppy / trash(3);

P = round(solve(1 / Ws), 3);
P = sort(P);
disp("P1 = "); disp(P(1));
disp("P2 = "); disp(P(2));
disp("P3 = "); disp(P(3));

disp("3. Точка пересечения асимптот:")
sigma_asimp = (sum(P) / 3);
disp(sigma_asimp);

disp("5. Подходящая точка пересечения с действительной осью (находящаяся в пределах КГ):");
Cx = 1 / (s - P(1)) + 1 / (s - P(2)) + 1 / (s - P(3));
Cx = round(solve(Cx), 2);
Cx = max(Cx);
display(Cx);

disp("6. Определим границу устойчивости по коэффициенту усиления (точки пересечения с мнимой осью):")
syms lambda;
syms k;
D = subs(Ws, s, lambda);
D = vpa(D * k + 1, 3);
display(D);
[num, den] = numden(D);
num_lambda = subs(num, k, 1);
num_arr = sym2poly(num_lambda);
num = vpa(num / num_arr(3), 5);
% display(num);

% Таблица Рауса
% Костыль
DsCoeffs = vpa(coeffs(num), 4);

MSize = max(size(DsCoeffs));

c11 = DsCoeffs(MSize - 1);
c21 = DsCoeffs(MSize - 3);
c31 = 0;
c12 = DsCoeffs(MSize - 2);
c22 = DsCoeffs(4) * k;
c32 = 0;
c13 = c21 - round((c11 / c12), 3) * c22;
c23 = c31 - (c11 / c12) * c32;
c33 = 0;
c14  = c22 - (c12 / c13) * c23;
c24 = 0;
c34 = 0;

fprintf("i\t|1\t\t\t\t\t|2\t\t\t\t\t\n");
fprintf("1\t|%f\t\t\t|%f\t\t\t\n", c11, c21);
fprintf("2\t|%f\t\t\t|%s\t\t\t\n", c12, c22);
fprintf("3\t|%s\t|%f\t\t\t\n", c13, c23);
fprintf("4\t|%s\t\t|%f\t\t\t\n", c14, c24);

k = vpa(subs(c12 / (c11 * c22), k, 1), 5);
display(k);

disp("Коэффициент усиления системы на границе устойчивости:");
Kkr = vpa(k * num_coppy, 6);
display(Kkr);

syms jw;
equation = [c11, c12 1, Kkr];
equation = vpa(poly2sym(equation, jw), 5);
display(equation);

disp("Выражаем из мнимой части уравнения wkr^2");
wkr = vpa(sqrt(1 / c11), 5);
disp(wkr);

Ws = Ws * k;
[num, den] = numden(Ws);
Ws_tf = tf(sym2poly(num), sym2poly(den));
display(Ws_tf);
display(Ws);

W = Ws_tf;
rlocus(W);
disp("Нажмите любую кнопку чтобы продолжить");
pause;
disp("Возьмём коэффициент усиления в 90% относительно коэффициента усиления на графнице устойчивочти");
disp("K = ");
kek = vpa(0.9 * Kkr, 5);
disp(kek);

Ws = Ws * 0.9;
[num, den] = numden(Ws);
Ws_tf = tf(sym2poly(num), sym2poly(den));
W = Ws_tf / (1 + Ws_tf);

step(W);
