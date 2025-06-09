function w=hanning_car(M);
% w=(1-cos(2*pi*([0:M-1]'+0.5)/M))/2;
w=(1-cos(2*pi*([0:M-1]')/M))/2;