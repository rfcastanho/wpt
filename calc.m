

% some preliminary information about the transformer and load
n1 = 7;             % turns in primary
n2 = 7;             % turns in secondary
v1 = 20;            % RMS voltage applied to primary
Z2 = 3;             % load impedance
ix = 1;			% guesstimated no-load primary current

% open FEMM and initialize the problem
openfemm;
try
    vv=ver;
    opendocument([cd,'/calc.fem']);
catch
    opendocument('calc.fem');
end
mi_saveas('temp.fem');

% build approximate transformer impedance
% based on the nominal no-load primary current
mi_setcurrent('emitter',ix);
mi_setcurrent('receiver',0);
mi_analyze;
mi_loadsolution;
Zt=[0,0;0,0];
r=mo_getcircuitproperties('emitter');
Zt(1,1)=r(2)/ix;
r=mo_getcircuitproperties('receiver');
Zt(1,2)=r(2)/ix;
Zt(2,1)=Zt(1,2);
Zt(2,2)=Zt(1,1)*(n2/n1)^2;

% Create load impedance matrix and voltage vector
Zl=[0,0;0,Z2];          % load impedance matrix
v =sqrt(2)*[v1;0];      % applied voltage amplitude

% compute initial current estimate based on
% approximate transformer and load impedance
ic=(Zt+Zl)\v;

for k=1:100
    mi_setcurrent('emitter',ic(1));
    mi_setcurrent('receiver',ic(2));
    mi_analyze;
    mi_loadsolution;
	vt=[0;0];
	r=mo_getcircuitproperties('emitter');
	vt(1)=r(2);
	r=mo_getcircuitproperties('receiver');
	vt(2)=r(2);
	u=(v-vt-Zl*ic);
	ic=ic + (Zt+Zl)\u;
	disp(sprintf('%i  %f',k,abs(sqrt(u'*u))));
	if (sqrt(u'*u)<0.1)
		break;
	end
end

disp('load currents are:');
disp(ic);

closefemm
