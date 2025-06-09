close all
clear all
clc

%% 
CBF_current = 15; % Current CBF
CBF_basal = 12.5; % Basal CBF
AC_current = 0.15; % Current AC
AC_basal = 0.15; % Basal AC
AC_autoreg_gain = 1.5; % Autoregulation gain
AC_time_constant = 20; % Time constant
AC_sigmoid_bound_high = 0.75; % High sigmoid bound
AC_sigmoid_bound_low = 0.075; % Low sigmoid bound

ACDelta = calculate_cerebral_autoregulation(CBF_current, CBF_basal, AC_current, AC_basal, ...
                                            AC_autoreg_gain, AC_time_constant, ...
                                            AC_sigmoid_bound_high, AC_sigmoid_bound_low);

disp(['Change in Arterial Compliance (ACDelta): ', num2str(ACDelta)]);
