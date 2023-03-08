function CurvaLvL_J(yn, xn, Weigths, M, L)
%         yn : Señal observada
%         xn : referencia a estimar(señal deseada)
%         Weigths : trayectoria de los pesos
%         M  : largo de ventana del filtro


% % Estimamos Ryx e Ry en funcion de la observaciones y los largos
    Ry_ = zeros(M,M); Ryx_ = zeros(M,1);
    for i = M : L; Ry_ = Ry_ + yn(i:-1:i-M+1).*yn(i:-1:i-M+1)'/L; end
    for i = M : L; Ryx_ = Ryx_ + yn(i:-1:i-M+1).*xn(i)'/L; end
    %1. Estimar R utilizando ventanas de largo M > K

    % Calculo de J(w1,w2)
    % Calculamos  la curva de aprendizaje
    nW = 20;
    [W1, W2] = meshgrid(linspace(-nW, nW), linspace(-nW, nW));
    
    J = zeros(size(W1));
    %Recordar: J(W) = sigma2x + W'Ryx - Ryx'W + W'RyW
    for idx = 1:length(W1)^2
        W = [W1(idx); W2(idx)];
        %J(idx) = var(xn) + W'*Ryx_ - Ryx_'*W + W'*Ry_*W;
        J(idx) = var(xn) + Ryx_'*W - W'*Ryx_ + W'*Ry_*W;
    end
    
    
    
    figure()
    contour(W1, W2, J, 50);
    hold on
    plot(Weigths(1,:), Weigths(2,:),'-*');
 
  	%legend("", "$w_1$","", "$w_2$",'Interpreter','latex')
    xlabel('$w_1$','interpreter','latex','FontSize',17);
    ylabel('$w_2$','interpreter','latex','FontSize',17);
    %title("Trayectoria de los pesos para distintos pesos inciales")

end