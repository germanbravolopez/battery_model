function z = Intensidad_inicial(t)
global i_salida A
z1=[i_salida/A*ones(1,10) zeros(1,40) -i_salida*0.75/A*ones(1,10)]; % 60 s
z = [i_salida*2/A*ones(1,10) -i_salida*2/A*ones(1,10) zeros(1,80) ...
    i_salida/A*ones(1,10) zeros(1,80) z1 z1 z1 z1 z1 z1 z1 z1 z1 ...
    zeros(1,70)];
% z = [i_salida/A*ones(1,20) -i_salida*2/A*ones(1,10) zeros(1,80) ...
%     i_salida/A*ones(1,10) zeros(1,80) z1 zeros(1,280) z1 ...
%     zeros(1,200)];
figure(7)
plot(z*A)
title('Intensidad de salida')
xlabel('Tiempo, t [s]')
ylabel('Corriente, [A]')

end