function [xn, un, vn, varR] = fun(sn, SNR_dB)
    L = length(sn);
    varS = var(sn);
    varR = varS / 10^(SNR_dB/10);  
    vn = sqrt(varR)*randn(L, 1);

    xn = sn + vn;

    %MIC-2: v(n): u(n) = 0.8 v(n)+ 0.2 v(n−1)−0.1 v(n−2).
    h = [0.8 0.2 -0.1];
    un = filter(h, 1, vn);

end