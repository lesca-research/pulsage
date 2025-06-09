function ACDelta = calculate_cerebral_autoregulation(CBF_current, CBF_basal, AC_current, AC_basal, ...
                                                     AC_autoreg_gain, AC_time_constant, ...
                                                     AC_sigmoid_bound_high, AC_sigmoid_bound_low)
    % Calculate the change in arterial compliance (ACDelta) as an indicator of cerebral autoregulation.
    %
    % Parameters:
    % - CBF_current: Current cerebral blood flow.
    % - CBF_basal: Basal (default) cerebral blood flow.
    % - AC_current: Current arterial compliance.
    % - AC_basal: Basal (default) arterial compliance.
    % - AC_autoreg_gain: Gain value for AC autoregulation.
    % - AC_time_constant: Time constant for AC autoregulation.
    % - AC_sigmoid_bound_high: High bound of the AC autoregulation sigmoid curve.
    % - AC_sigmoid_bound_low: Low bound of the AC autoregulation sigmoid curve.
    %
    % Returns:
    % - ACDelta: Change in arterial compliance as an indicator of cerebral autoregulation.

    % Normalize CBF relative to basal CBF value
    normalizedCBF = (CBF_current - CBF_basal) / CBF_basal;
    
    % Determine the amplitude of the sigmoid and sigmoid slope constant based on normalized CBF
    if normalizedCBF < 0
        centralSigmoidValue = AC_sigmoid_bound_high;
    else
        centralSigmoidValue = AC_sigmoid_bound_low;
    end
    
    sigmoidSlopeConstant = centralSigmoidValue / 4.0;
    
    eToPower = exp((AC_autoreg_gain * normalizedCBF) / sigmoidSlopeConstant);
    sigma = ((AC_basal + (centralSigmoidValue / 2.0)) + ...
             ((AC_basal - (centralSigmoidValue / 2.0)) * eToPower)) / ...
            (1 + eToPower);
    
    % Calculate change in arterial compliance per unit time
    ACDelta = (sigma - AC_current) / AC_time_constant;
end
