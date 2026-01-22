load variables_menos_x_y_ebx_eby.mat

EBx=50
EBy=50

Eminima=0.01
n=250

load escenario_centro.mat
for k=1:15,

    x =X(k,:);
    y=Y(k,:);
%x=rand(1,n)*x_max;
%y=rand(1,n)*y_max;


[FND(k) HND(k) LND(k) noCH(k) noCHFND(k) noCHHND(k)] = saez_best_v3_variables_de_entrada_simples(x,y,x_max,y_max,n,stop,p,EBx,EBy,Eo,ETX,ERX,Efs,Eamp,EDA,packetSize,controlPacketSize,Eminima)

end