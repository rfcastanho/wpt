% Algoritmo para determinar correntes e tensões em um transformador
% monofásico.

% Informações Iniciais (recuperadas de outras rotinas)
clear
clc
n1 = 80;
n2 = 80;
R1 = 0.2;
R2 = 0.2;
L1 = 186.65e-6; %120e-6;
L2 = 186.65e-6; %120e-6;
C1 = 601.95e-9;
C2 = 601.95e-9;
v1 = 9.20;          % Tensão eficaz da comp. fundamental na saída do inversor. 
w = 94342.12173;
LL= 0;            % Impedância da carga.
CL = inf;
RL = 10;
ZL = RL + 1i*w*LL + 1/(1i*w*CL);

ip_vazio = 74.5;	% Corrente do primário, a vazio.

% open FEMM and initialize the problem
openfemm;
try
    vv=ver;
    opendocument([cd,'/circuito_magnetico.fem']);
catch
    opendocument('circuito_magnetico.fem');
end
mi_saveas('tempmag.fem');

% Aqui calcula impedância do transformador com base em ip_vazio

mi_setcurrent('primario',ip_vazio); %Força ip_vazio no circuito primário do modelo do FEMM.
mi_setcurrent('secundario',0);      %Força corrente nula no circuito secundário do modelo do FEMM.
mi_analyze;
mi_loadsolution;
% p = mo_getcircuitproperties('primario');
% fluxo_mutuo = abs(p(3));
% i1 = abs(p(1));
% M = n2*fluxo_mutuo/i1;
% Zmag=[0,0;0,0];
% Zmag(1,1)= R1 + 1i*w*L1 + 1/(1i*w*C1);
% Zmag(1,2)=-1i*w*M;
% Zmag(2,1)=Zmag(1,2);
% Zmag(2,2)= R2 + 1i*w*L2 + 1/(1i*w*C2);
Zt=[0,0;0,0];
r=mo_getcircuitproperties('primario');
Zmag(1,1)=r(2)/ip_vazio;
r=mo_getcircuitproperties('secundario');
Zmag(1,2)=r(2)/ip_vazio;
Zmag(2,1)=Zmag(1,2);
Zmag(2,2)=R2 + 1i*w*L2 + 1/(1i*w*C2) ; %Zmag(1,1)*(n2/n1)^2;

% Create load impedance matrix and voltage vector
Zcarga = [0,0;0,ZL];          % load impedance matrix
v = sqrt(2)*[v1;0];           % applied voltage amplitude

% compute initial current estimate based on
% approximate transformer and load impedance
ic=(Zmag+Zcarga)\v;

for k=1:30
    mi_setcurrent('primario',ic(1));
    mi_setcurrent('secundario',ic(2));
    mi_analyze;
    mi_loadsolution;
	vt=[0;0];
	r=mo_getcircuitproperties('primario');
	vt(1)=r(2);
	r=mo_getcircuitproperties('secundario');
	vt(2)=r(2);
	u=(v-vt-Zcarga*ic);
	ic=ic + (Zmag+Zcarga)\u;
	disp(sprintf('%i  %f',k,abs(sqrt(u'*u))));
	if (sqrt(u'*u)<0.1)
		break;
	end
end

disp('load currents are:');
disp(ic);

closefemm
