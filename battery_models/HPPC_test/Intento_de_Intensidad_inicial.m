function z = Intensidad_inicial(t) % 19200s son 2560
global i_salida A
z1=[i_salida/A*ones(1,10) zeros(1,40) -i_salida*0.75/A*ones(1,10)]; % 60 segundos
z = [i_salida/A*ones(1,460) -i_salida/A*ones(1,460) ...
    zeros(1,480) i_salida/A*ones(1,14) ...
    zeros(1,480) z1 i_salida/A*ones(1,14) ...         
    zeros(1,480) z1 i_salida/A*ones(1,14) ...   
    zeros(1,38)];                     % 2560ptos          
%z = [i_salida/A*ones(1,23) zeros(1,3600) z1 z1 z1 z1 z1 i_salida/A*ones(1,23)];
figure(7)
plot(z)
title('Intensidad de salida')
xlabel('Tiempo, t [s]')
ylabel('Corriente, [A]')

end