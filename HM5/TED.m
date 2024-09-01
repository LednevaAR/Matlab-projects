function result = TED(Signal, Ts)
    Tm = 2 * Ts;
    Sig1 = real(Signal);
    Sig2 = imag(Signal);
    Xi = 1;
    BnTs = 0.5;
    Kd = 2*pi;
    K0 = 1;
    p = 0;
    for m = 2 : length(Signal) / 2
        ted_output(m) = Sig1(2 * m - 2) * (sign(Sig1(2 * (m - 1) - 1)) - sign(Sig1(2 * m - 1))) + Sig2(2 * m - 2) * (sign(Sig2(2 * (m - 1) - 1)) - sign(Sig2(2 * m - 1)));
        Kp = 2 * Xi * (BnTs / (Xi + 1 / (4 * Xi))) / (Kd * K0);
        Ki = (BnTs / (Xi + 1 / (4 * Xi))) ^ 2 / (Kd * K0);
        ted_filtered(m) = Kp * ted_output(m) + Ki * ted_output(m) + p;
        p = ted_filtered(m) - Kp * ted_output(m);
    end
    result = ted_filtered;
end