% clear alld
% close all
% 
% k = 2*pi*1000/343;
% a = 0.1;
% r = 10^(-5);
% par = abs( sin(1/2*k*r*(sqrt(1+(a/r)^2)-1)) );
% p = 2*1.21*343*par;
% dB = 20*log10(p/sqrt(2)/(20*10^(-6)));
% 
% p_1 = 10^(148.58/20)*sqrt(2)*(20*10^(-6));
% par_1 = p_1/(2*1.21*343);

% inverse 
% pressure = sqrt(2)*20*10^(-6)*10^(dB/20)

% %% Spherical Radiator
% rho = 1.21;
% c = 343;
% u = 1;
% a = 1;
% r = 1; % 어디서부터인지 정확히 확인할 것
% freq = 100/1.8318;
% ka = 2*pi*freq/c*a;
% p = rho*c*u*a/r*abs(1j*ka/(1+1j*ka));

%% 

% rho = 1.21;
% c = 343;
% freq = 10:10:5000;
% omega = 2*pi*freq;
% x = 0.1;
% 
% imp = rho.*c.*(1j*omega*x./(c+1j*omega*x));
% 
% plot(freq,abs(imp))

% c = 343;
% D1 = 50*10^(-3);
% D2 = 30*10^(-3);
% kl = 3.22211336044481;
% N = D1/D2;
% l_p = 500*10^(-3);
% f = c/(2*pi)*(kl/l_p);


%%
% x = 0:.1:10;
% y = x.*x + randn(size(x));
% w = [100; ones(length(x)-2,1); 100];
% x = x(:);
% y = y(:);
% w = w(:);
% %plot data
% plot(x,y,'.');
% %fit
% ft = fittype('poly2');
% cf = fit(x,y,ft,'Weight',w);
% % Plot fit
% hold on
% plot(cf,'fit',0.95);


% c = 343;
% m = 1;
% n = 0:10;
% L = 1;
% f = c/(2*pi).*sqrt(m^2/4 + (n*pi./L).^2);

% A = 720;
% B = [675, 707, 701];
% A = 670;
% B = 680;
% E = abs(B-A)./A*100;
% E_avg = mean(E);

% F2 = gradient(FDCurveCase2);
% F4 = gradient(FDCurveCase4);
% 
% p_1 = 0.5; p_2 = 2.4; r_1 = 0.7; r_2 = 1.5;
% fun = @(r) ((p_2-p_1)./(r_2-r_1).*(r-r_1)+p_1).*r;
% fun = @(r) ((p_2-p_1)./(r_2-r_1).*(r-r_1)+p_1);
% q_1 = integral(fun,r_1,r_2);
% q_2 = (r_2-r_1)*( (1/3*r_2+1/6*r_1)*p_2 + (1/6*r_2+1/3*r_1)*p_1);
% q_2 = 1/2*(r_2-r_1)*(p_2+p_1);

F = 1e3; % [N]
d = 6e-3; % [m] 
A = pi*d^2/4; %[m^2]
sigma = F/A; %normal stress
E = 193e9;  % Young Modulus [Pa]
epsilon = sigma/E; %strain
L = 60;% length of thread [mm]
delta = epsilon*L;% [mm]

%% Relative error
A = FEM_KaRatio(FEM_locs(3));
B = Lee_KaRatio(LPM_locs(2));
error = abs((A-B)/A)*100;

% %% Interpolation
% start_num
% step
% end_num
% Int_range = (-start_num:step:end_num);
% Int_Load = interp1(Displacement,Load,Int_range);



    
