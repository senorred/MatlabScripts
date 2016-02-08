x=[-2:.001:2];
y=(sqrt(cos(x)).*cos(200*x)+sqrt(abs(x))-0.7).*(4-x.*x).^0.01;
plot(x,y,'r')
title('Love you!')