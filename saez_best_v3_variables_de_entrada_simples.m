
%Licencia CC
%Citar

%Yuste-Delgado, A. J., Cuevas-Martinez, J. C., & Trivino-Cabrera, A. (2022). 
%Statistical normalization for a guided clustering type-2 fuzzy system for wsn. 
%IEEE Sensors Journal, 22(6), 6187-6195.
%https://ieeexplore.ieee.org/abstract/document/9707772
%DOI: 10.1109/JSEN.2022.3150066


function [FND HND LND noCH noCHFND noCHHND] = saez_best_v3_variables_de_entrada_simples(x,y,x_max,y_max,n,stop,p,EBx,EBy,Eo,ETX,ERX,Efs,Eamp,EDA,packetSize,controlPacketSize,Eminima)

%x: vector con las abcisas en las que se encuentran los sensores
%y: vector con las ordenadas en las que se encuentran los sensores

%x_max: dimensión máxima x
%y_max: dimensión máxima y

%n:  Numer of nodes in the experiment
%stop: número de muertes para llegar al final en % sobre n
%p: variable p del método de LEACH

%EBx: coordenada x de la estación base
%EBy: coordenada y de la estación base

%Eminima: energía mínima antes de que un nodo se considere muerto

disp([ '-- BS x ', num2str(EBx),' BS y ', num2str(EBy)]);

% maximum number of rounds
    rmax=99999; 


%The seed for random function is stablish
seed=2;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));

%Energy Model (all values in Joules)

% Eo Initial Energy

% ETX  default value 50 nJ/bit 


% ERX default value 50 nJ/bit [1]

% Efs Transmit Amplifier types,
% Free space power loss coeficient (d^2): default value 10 pJ/bit/m2 [1]
% Some other authors use 100 pJ/bit/m2

% Eamp Multipath power loss coeficient (d^4): default value 0.0013e-12 pJ/bit/m4 [1]

% EDA Data Aggregation Energy, default value: 5e-9;%nJ/bit/signal [1]

% packetSize Size of the data packets sent in the experiment, default 4000 bits [1]

%controlPacketSize Size of the data packets sent in the experiment, default 4000 bits [1]





%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%% Auxiliar variables

%Computation of do
do=sqrt(Efs/Eamp)
MAXCAM=sqrt(x_max^2+y_max^2);

d1=sqrt( (0-EBx )^2 + (0-EBy )^2 );
d2=sqrt( (x_max- EBx )^2 + (0-EBy)^2 );
d3=sqrt( (0-EBx)^2 + (y_max-EBy)^2 );
d4=sqrt( (x_max-EBx )^2 + (y_max-EBy )^2 );

nodos_vivos=1:n; %(alive nodes)

%Cálculos de las distancias iniciales a la EB

maxdistancia=0;
for i=1:length(nodos_vivos)
	%aprovechamos el bucle para calcular distancias a la estacion base, indica que se normaliza por la distancia mayor 
	D(nodos_vivos(i)) =sqrt( (x(nodos_vivos(i))-EBx )^2 + (y(nodos_vivos(i))-EBy )^2 );
	if (D(nodos_vivos(i))>maxdistancia)
		maxdistancia=D(nodos_vivos(i));
	end;
end	

maxdistancia = max([maxdistancia,d1,d2,d3,d4]);

%Todos los nodos mandan un mensaje a la EB para conocer su localización, además todos los demás nodos lo leen 
%Para que lo vean todos se debe enviar a la mayor de las distancias o MAXCAM o maxdistancia

distancia_maxima = max(MAXCAM,maxdistancia);

E =Eo.*ones(1,n);
for i=1:length(nodos_vivos)
    E(nodos_vivos(i))= E(nodos_vivos(i)) - disminuir_energia_envio_mensaje(0,0,distancia_maxima,0,controlPacketSize,do,ETX,Eamp,Efs);
end



%además cada nodo deber recibir n-1 mensajes
for i=1:length(nodos_vivos)
	E(nodos_vivos(i))= E(nodos_vivos(i)) - (n-1)*disminuir_energia_recibir_mensaje(controlPacketSize,ERX);
end


dead=0;
nclusteres=fix(n*p);

cambio_skip=0;

%cargamos la matriz con las salidas del sistema difuso

load matriz_saez_v3


skip= fix(n * p); 
disp('********');
disp(skip);

FND=-1;
HND=-1;
LND=-1;

skip_count=0;

rCH_trheshold=  n; %nclusteres;

%creamos vectores con los valores de las variablesdel sistema difuso
rCH=zeros(1,n);
vecesCH=zeros(1,n);


orden=1:nclusteres;

nodos_asociados=ones(n,1);

factor=sqrt((x_max*y_max)/(n*p)/(pi));

distancias_nodos_asociados=factor * ones(n,1);
noCH=0;

energia_consumida=zeros(n,1);
energia_previa= Eo.*ones(1,n); 

    
p_1_energia_ch = 0.1 * Eo;
p_2_energia_ch = 0.9 * Eo;

            p_1_rch_ch = 0;
			p_2_rch_ch = n;

   			p_1_vecesCH_ch = 0;
			p_2_vecesCH_ch = 1;


total_CH=1;
chances(n)=0;
energia_ch(n)=0; %variables 1 y 3
        vecesCH_ch(n)=0;
        cl(n)=0;
        rch_ch(n)=0; %variable 4
        rnoCH_ch(n)=0; %variable 4
         

for r=0:rmax,
    if mod(r,1000)==0
		disp([r length(nodos_vivos)  dead FND HND LND])
	end

	%primero eliminamos los nodos sin energía
    i=1;
    while (i<=length(nodos_vivos))
    %for i=1:length(nodos_vivos) Nodos vivos los voy eliminados en el bucle
    %por lo que no puede ser un for
        	
		
		if (E(nodos_vivos(i))<=Eminima)
			%MAP.S(nodos_vivos(i)).E
			%eliminar el nodo vivo de la lista
           
			if (i==1)
				nodos_vivos=nodos_vivos(2:length(nodos_vivos));
				i=i-1;
			elseif i==length(nodos_vivos)
				nodos_vivos=nodos_vivos(1:length(nodos_vivos)-1);
				i=i-1;
			else
				nodos_vivos=[nodos_vivos(1:i-1) nodos_vivos(i+1:length(nodos_vivos))];
				i=i-1;
			end
			
			
			dead=dead+1;
			if (dead==1)
					FND=r;
                     noCHFND=noCH;
					 cambio_skip=1;
                     skip_count=skip;
                     skip=fix(n*20 /100); %LO NUEVO
			end
			if (dead==fix(n*0.5))
				HND=r;
                 noCHHND =noCH;							
            end
           	if (dead==fix(n*0.75))
				 cambio_skip=1;
                 skip=2;				
            end
        end
		
        i=i+1;
    end
    
	%comprobar que no se ha acabado
	if(n-dead<=(n*stop))
		LND=r;
		break;
	end
	
	%Ahora buscamos si los clusters tienen energía para otro skip 
    %no se haría en la iteración inicial ni en la final
	
	if (r>0) & (skip_count< skip)    
	%si ha habido alguna muerte de un cluster se ha cambiado nodos vivos...
	  for i=1:nclusteres,
	    if E((orden(i)))< Eminima
		 % disp([r,skip_count,skip]);
            skip_count = skip; %fin de skip         
			break;
		end
	  end
	end
	
    %ahora se realiza el nuevo skip
 	
	if (skip_count==skip) || (r==0) || (cambio_skip==1)
	
		skip_count=1;
		
		if (cambio_skip)==1
		%lectura del nuevo skip enviado por la EB
			for i=1:length(nodos_vivos)
				E(nodos_vivos(i))= E(nodos_vivos(i)) - disminuir_energia_envio_mensaje(0,0,distancia_maxima,0,controlPacketSize,do,ETX,Eamp,Efs);
			end	
			cambio_skip=0;
		end
	
		%%%se eligen los nuevos clusteres
        
		chances=zeros(n,1);
        
		for i=1:length(nodos_vivos)
       
				if(rnoCH_ch(nodos_vivos(i))> p_2_rch_ch)
                    rnoCH=1;
                elseif (rnoCH_ch(nodos_vivos(i))< p_1_rch_ch)
                    rnoCH=0;
                else
                     rnoCH=(rnoCH_ch(nodos_vivos(i))-p_1_rch_ch)/(p_2_rch_ch-p_1_rch_ch);
 
                end
      
          
         %se normaliza entre p1 y p2
         
         		if E(nodos_vivos(i))> p_2_energia_ch
					energia = 1;
                elseif E(nodos_vivos(i)) < p_1_energia_ch
                    energia=0;
				else
					energia= (E(nodos_vivos(i)) -p_1_energia_ch)/(p_2_energia_ch-p_1_energia_ch);
                end					
         
           
         	if vecesCH(nodos_vivos(i)) > p_2_vecesCH_ch
					vch = 1;
                elseif vecesCH(nodos_vivos(i)) < p_1_vecesCH_ch
                    vch=0;
				else
					vch= (vecesCH(nodos_vivos(i)) -p_1_vecesCH_ch)/(p_2_vecesCH_ch-p_1_vecesCH_ch);
                end					
         
                
           
			X = [ E(nodos_vivos(i))/Eo vch energia rnoCH] ;

				errores = X<=0;
				X(errores)=0;
				%ol
                ol=z_ol(fix(X(1)*9)+1,fix(X(2)*9)+1,fix(X(3)*9)+1,fix(X(4)*9)+1);
                %or:
                or=z_or(fix(X(1)*9)+1,fix(X(2)*9)+1,fix(X(3)*9)+1,fix(X(4)*9)+1);
                
                %If the node has been a CH in the previous round it chooses the
                %lower probability
                if (cl(nodos_vivos(i))==1) || (cl(nodos_vivos(i))==2) % || rnoCH <5
                    chance = ol;
                else
                    chance = or;
                end
				
				% aleatoriedad para ser CH
		  temp_rand = rand;
				if (temp_rand > chance)
					chances(nodos_vivos(i))=0;
                    rnoCH_ch(nodos_vivos(i))=rnoCH_ch(nodos_vivos(i))+skip; 

					cl(nodos_vivos(i))=0;
                    rCH(nodos_vivos(i))=rCH(nodos_vivos(i))+skip;
                    
                else
                    
                    vecesCH(nodos_vivos(i))=vecesCH(nodos_vivos(i)) +1;
					total_CH=total_CH+1;
                    vecesCH_ch(nodos_vivos(i))=vecesCH(nodos_vivos(i));
        
                                                                                   
					chances(nodos_vivos(i))=chance + 0.000001 * rand *E(nodos_vivos(i));  % + 0.00001 * rand;
						%descontamos el envio del CH
					energia_ch(nodos_vivos(i)) = E(nodos_vivos(i));
    
                    rch_ch(nodos_vivos(i))=rCH(nodos_vivos(i));
                     rCH(nodos_vivos(i))=0;
    
                    
                     
					E(nodos_vivos(i))= E(nodos_vivos(i)) - disminuir_energia_envio_mensaje(0,0,MAXCAM,0,controlPacketSize,do,ETX,Eamp,Efs);
					
					%luego todos los nodos tienen que eliminar la escucha de todos los clusteres
        
                end
				
				
		end

		%Aquí habría que mirar cuantos clusteres hay
		clusteres_potenciales = sum(chances(nodos_vivos)>0);
        nclusteres = min(fix(n*p), clusteres_potenciales);
		[indices orden] =sort(-chances(nodos_vivos)); 
      
		orden2= nodos_vivos(orden); %para resetear nodos que ahora no serán noCH
        orden = nodos_vivos(orden(1:nclusteres));
        
		if (nclusteres ==0)
            noCH=noCH+1;
            disp(r);
            skip_count=skip;
		else
			p_1_energia_ch = prctile(energia_ch(orden2(1:clusteres_potenciales)),10);
			p_2_energia_ch = prctile(energia_ch(orden2(1:clusteres_potenciales)),90);
            
           p_1_rch_ch = prctile(rnoCH_ch(orden2(1:clusteres_potenciales)),10);
			p_2_rch_ch = prctile(rnoCH_ch(orden2(1:clusteres_potenciales)),90);

            
            if (p_2_rch_ch ==0)
            
             p_2_rch_ch=1;
            end
            
            
            if (p_1_rch_ch == p_2_rch_ch)
                disp([r,p_1_rch_ch]);
                p_1_rch_ch=0;
            end
            
            
            
   			p_1_vecesCH_ch = prctile(vecesCH_ch(orden2(1:clusteres_potenciales)),10);
			p_2_vecesCH_ch = prctile(vecesCH_ch(orden2(1:clusteres_potenciales)),90);
    
            if (p_1_vecesCH_ch == p_2_vecesCH_ch)
                disp([r,p_1_vecesCH_ch]);
                p_1_vecesCH_ch=0;
            end
            
            
            if clusteres_potenciales ==1
                disp('cero');
                p_1_energia_ch = 0;
                p_1_vecesCH_ch=0;
                p_1_rch_ch=0;
                
            end
            
		end;
        
        
		%añadimos solo los clusteres finales 
		for i=1:nclusteres
			cl((orden(i)))=1;
			rCH((orden(i)))=0;
            rnoCH_ch(orden(i))=0; %variable 4
    
            %de esto no estoy seguro
            %actualizo las distancias y los nodos asociados
            nodos_asociados(orden(i))=1;
            distancias_nodos_asociados(orden(i))=0;
            
        end
		
        %eliminamos los cl por si en la contienda anterior estuvieran a 1.
      %El cl=2 significa que ha sido elegido cluster pero desechado, así
      %entra como nodo normal y se adhiere a otro cluster
      
   	 for i=nclusteres+1:clusteres_potenciales
			cl((orden2(i)))=2;
            rCH((orden2(i)))=-skip; %luego se añade uno abajo
            nodos_asociados(orden2(i))=1;
            distancias_nodos_asociados(orden2(i))=0;
     end

       
     
	
	else
	%se mantienen los clusteres
		skip_count = skip_count +1;
	
	
	end

	%gasto energetico de enviar los avisos de CH, 
	%gasto energetico de enviar datos al cluster asociado
	%gasto energetico de enviar el CH datos a la EB, se incluye la agregación
	
	%hay que ver el número de nodos asociados a cada CH
    %al inicio

	for i=1:length(nodos_vivos)
		%quitar gasto energético del aviso de los CH, sólo la primera vez
          if (skip_count ==1)    
		      E(nodos_vivos(i))= E(nodos_vivos(i)) - disminuir_energia_recibir_mensaje(clusteres_potenciales * controlPacketSize,ERX);
		  
          end
		if (cl(nodos_vivos(i))~=1)
            %El rCH se debería actualizar sólo la primera vez que se entra
            %...
			%MAP.S(nodos_vivos(i)).rCH=MAP.S(nodos_vivos(i)).rCH+1;
			min_dista=10e6; %%OJOJOJOJ
			
            cluster_asignado=0;
			for k=1:nclusteres,
			   if (cl((orden(k)))==1)
					distancia= sqrt( (x(nodos_vivos(i))-x((orden(k)))) ^2 + (y(nodos_vivos(i))-y((orden(k))))^2);
					if (min_dista>distancia)
					  cluster_asignado = orden(k);					  
	%				  k_asignado=k;
					  min_dista=distancia;
					end
				end
            end
		

            if (cluster_asignado~=0)
                
                %se actualiza Ctr la primera vez
                if (skip_count==1)
                    nodos_asociados(cluster_asignado)=nodos_asociados(cluster_asignado)+1;
                    distancias_nodos_asociados(cluster_asignado)=distancias_nodos_asociados(cluster_asignado)+min_dista;
                    
                %se envían control para unirse al cluster, solo primera vez.
                E(nodos_vivos(i))= E(nodos_vivos(i)) - disminuir_energia_envio_mensaje(x(nodos_vivos(i)),y(nodos_vivos(i)),x(cluster_asignado),y(cluster_asignado),controlPacketSize,do,ETX,Eamp,Efs);

                end
                
                
                %se envían datos. 
                E(nodos_vivos(i))= E(nodos_vivos(i)) - disminuir_energia_envio_mensaje(x(nodos_vivos(i)),y(nodos_vivos(i)),x(cluster_asignado),y(cluster_asignado),packetSize,do,ETX,Eamp,Efs);
                %el cluster asignado recibe los datos
                E(cluster_asignado)= E(cluster_asignado) - disminuir_energia_recibir_mensaje(packetSize,ERX);
                %el cluster tiene que agregar los datos
                E(cluster_asignado)= E(cluster_asignado)- disminuir_energia_agregacion(packetSize,EDA);
                
            else
                %envio directo a la EB solo datos
                E(nodos_vivos(i))= E(nodos_vivos(i)) - disminuir_energia_envio_mensaje(x(nodos_vivos(i)),y(nodos_vivos(i)),EBx,EBy,packetSize,do,ETX,Eamp,Efs);
               
            end
        end 
	  end

	
      
	  %envío a la EB
	  for i=1:nclusteres,
       E((orden(i)))=  E((orden(i))) - disminuir_energia_envio_mensaje(x((orden(i))),y((orden(i))),EBx,EBy,packetSize,do,ETX,Eamp,Efs);
      end
	
	
      %actualizamos las variables  para un solo skip
      if (skip_count==1)
        for i=1:length(nodos_vivos),          
            energia_consumida (nodos_vivos(i))= E(nodos_vivos(i))- energia_previa(nodos_vivos(i));
            energia_previa(nodos_vivos(i))=E(nodos_vivos(i));
         end
      end
      
%     ssctr = distancias_nodos_asociados ./ nodos_asociados;  
        

end %r
	
cprintf('_red',['- FND=', num2str(FND),' HND=', num2str(HND),' LND=', num2str(LND),' noCH=', num2str(noCH),' EB (',num2str(EBx),',',num2str(EBy),')\n']);
cprintf('blue',[' V3 - SKIP ',num2str(skip), ' Control ', num2str(controlPacketSize),' \n']);
disp('');

end %function




function perdidas= disminuir_energia_envio_mensaje(xa,ya,xb,yb,longitud,do,ETX,Eamp,Efs)

 

 distancia = sqrt((xa-xb)^2+(ya-yb)^2);
 
 if (distancia>do)
    perdidas=( ETX*(longitud) + Eamp*longitud*(distancia^4));
 else
    perdidas=( ETX*(longitud) + Efs*longitud*(distancia^2)); 
 end

 end


function perdidas= disminuir_energia_recibir_mensaje(longitud,ERX) 

 perdidas= ERX*longitud;
end


function perdidas= disminuir_energia_agregacion(longitud,EDA)

perdidas= longitud * EDA;

end

