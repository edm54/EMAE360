function [gamma] =  calc_gamma(temp)
    N2_high = [.2896e1, .15155e-2, -.57235e-6, .99807e-10, -.652e-14];
    N2_low = [.3675e1, -.12082e-2, .23240e-5, -.6322e-9, -.2258e-12];
    O2_high = [.362e1, .7362e-3, -.1965e-6, .362e-10, -.2895e-14];
    O2_low = [.36256e1, -.18782e-2, .70555e-5, -.6764e-8, .21556e-11];

    if (temp>1000)
        N = N2_high;
        O = O2_high;
    else
        N = N2_low;
        O = O2_low;
    end    
    
    cp_n = temp_to_gamma(N, temp);
    cp_o = temp_to_gamma(O,temp);
    
    cp = .2095 * cp_o + .7905 * cp_n;
    R = .287;
    cp = cp*R;
    cv = cp - R;
    function [cp] = temp_to_gamma(a, temp)
        cp =  a(1) + a(2)*temp + a(3)*temp^2 + a(4)*temp^3 + a(5) * temp^4;
    end
    gamma = cp/cv;

end