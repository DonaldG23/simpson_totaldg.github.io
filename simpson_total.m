function simpson_total()
    % Definición de la función
    f = @(x) 1 + 2*x - 3*x.^2 + 4*x.^3 - 5*x.^4 + 6*x.^5;
    a = 0;
    b = 1;

    % Número de subintervalos
    n1 = 6; % múltiplo de 2 para Simpson 1/3
    n2 = 6; % múltiplo de 3 para Simpson 3/8

    % Simpson 1/3
    h1 = (b - a) / n1;
    x1 = a:h1:b;
    y1 = f(x1);
    simpson13 = y1(1) + y1(end);
    for i = 2:n1
        if mod(i-1,2) == 0
            simpson13 = simpson13 + 2 * y1(i);
        else
            simpson13 = simpson13 + 4 * y1(i);
        end
    end
    simpson13 = (h1 / 3) * simpson13;

    % Simpson 3/8
    h2 = (b - a) / n2;
    x2 = a:h2:b;
    y2 = f(x2);
    simpson38 = y2(1) + y2(end);
    for i = 2:n2
        if mod(i-1,3) == 0
            simpson38 = simpson38 + 2 * y2(i);
        else
            simpson38 = simpson38 + 3 * y2(i);
        end
    end
    simpson38 = (3*h2 / 8) * simpson38;

    % Integral exacta
    syms x
    fx = 1 + 2*x - 3*x^2 + 4*x^3 - 5*x^4 + 6*x^5;
    exacta = double(int(fx, x, a, b));

    % Cuarta derivada media
    d4 = matlabFunction(diff(fx, x, 4));
    x_eval = linspace(a, b, 1000);
    prom_d4 = mean(abs(d4(x_eval)));

    % Error de truncamiento estimado (Simpson 1/3)
    Et = -((b - a)^5 / (180 * n1^4)) * prom_d4;

    % Error relativo porcentual (Simpson total vs exacta)
    total_aprox = (simpson13 + simpson38) / 2;
    error_rel = abs((exacta - total_aprox) / exacta) * 100;

    % Resultados
    fprintf('Integral aproximada (Simpson 1/3): %.8f\n', simpson13);
    fprintf('Integral aproximada (Simpson 3/8): %.8f\n', simpson38);
    fprintf('Integral aproximada total: %.8f\n', total_aprox);
    fprintf('Valor medio de la cuarta derivada: %.8f\n', prom_d4);
    fprintf('Error de truncamiento estimado: %.8f\n', Et);
    fprintf('Error relativo porcentual: %.8f %%\n', error_rel);
end
