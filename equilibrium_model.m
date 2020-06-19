%% Evaluating Validity of Equilibrium Reaction Expression

%% Equilibrium Reaction Model Regression

% Extracted Raw Data from Figure S1 (Wu and Chien)
    xl = 4;
    xr = 166;
    yb = 358;
    yt = 287;
    x_pix = xr-xl;
    y_pix = yb-yt;
    point = [134 236 288 346 401;306 247 136 115 107];
    temp = 1./(0.0024+(point(1,:)-xl)./x_pix*0.0001); % Experimental Data
    logkeq = (-6+(yb-point(2,:))./y_pix*0.4);

%Plotting Data Points from Figure S1
    figure(1); hold on;
    plot(1./temp,logkeq,'ko','linewidth',1.2);

% Plotting Wu and Chien Model
    T = linspace(min(temp),max(temp),1e5);
    oldmodel = -201.9 + 43088./T +0.22*T;
    plot(1./T,oldmodel,'r-.','linewidth',1.2);

% MATLAB lsqcurvefit Model
    modelfun = @(b,t) b(1) + b(2)./t +b(3).*t; 
    p = [-200;40000;0.2]; % initial guesses 
    x=zeros(1,3);
    [x(1,:),resnorm,resid,~,~,~,J] = lsqcurvefit(modelfun,p,temp,logkeq);
    rsquare = 1- resnorm; % R-square Value
    fit = x(1)+x(2)./T+x(3).*T;
    plot(1./T,fit,'b-.','linewidth',1.2)

% Miscellaneous Labels    
xlabel('$\bf\frac{1}{T}[K^{-1}]$','interpreter','latex')
ylabel('$\bf ln(K_{eq}) $','interpreter','latex')
legend('\bf Experimental Data','\bf Wu and Chien Model',...
    '\bf lsqcurvefit Model, $\bf R^2 = 0.9156$','interpreter','latex','location','best')
xlim([1/max(temp)-5e-6 1/min(temp)+5e-6])
set(gca,'fontsize',11)


%% Raw Data from Figure 3a (Bian et al.)

Temps = [352.855337387354;362.998937651110;372.975859037292;377.976225364496;382.954945084222;387.976958018902;392.977324346106;403.014689720859];
X = [0.0504780954239438;0.0589667284613133;0.101702605132208;0.101702605132208;0.100771250410240;0.0943848180310274;0.0724314567274856;0.0615110667457238];

% Converting Raw Data to Keq
Keq = X.^2.*(3-X)./(4*(1-X).^3)./12; % Partial Pressure Basis (12 bar)
Keq1 = X.^2.*(3-X)./(4*(1-X).^3)./15; % Partial Pressure Basis (15 bar)
Keq2 = X.^2.*(3-X)./(4*(1-X).^3); % Mole Fraction Basis 

T = Temps(1):0.001:Temps(8);
Tmat = [ones(size(T)) T 1./T];

% Plotting
figure(2); hold on
ylim([-9,-3.5])
plot(1./temp,logkeq,'ro','linewidth',1.2);
plot(1./Temps,log(Keq2),'ko','linewidth',1.2);
plot(1./Temps,log(Keq),'bo','linewidth',1.2);
xlabel('$\bf\frac{1}{T}$, $\bf K^{-1}$','interpreter','latex')
ylabel('$\bf ln(K_{eq}) $','interpreter','latex')
legend('\bf Wu and Chien Figure S1 Data','\bf Bian et al. Figure 3a Data, Mole Fraction Basis','\bf Bian et al. Figure 3a Data, Partial Pressure Basis','interpreter','latex')
set(gca,'fontsize',11)

figure(3)
plot(Temps,X.*100,'ko','linewidth',1.2)
ylim([0 30])
xlabel('\bf Temperature, K','interpreter','latex')
ylabel('\bf Methanol Conversion, \%','interpreter','latex')
set(gca,'fontsize',11)




