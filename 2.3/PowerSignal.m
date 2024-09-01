function Power = PowerSignal(Signal)
    Power = mean(abs(Signal) .* abs(Signal));
end

