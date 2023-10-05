%% All

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420; %420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;

        p.d_V = 5.45;
        p.V0 = 0.000069;
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

patient = log10(curves);
p13 = readtable('p13.csv');
T13 = table2array(p13(:,4));
V13 = table2array(p13(:,10));
R13 = table2array(p13(:,11));
A13 = table2array(p13(:,12));
figure(1)
plot(t, patient,'LineWidth',0.5,'color',[37,37,37]/255);
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
scatter(T13,R13,6,'o','filled','MarkerEdgeColor','#CD2626','MarkerFaceColor','#CD2626','LineWidth',0.75,'Color','#CD2626');
hold off
ylim([1.5 5])
xlim([0 31])
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='	#CD2626';
box off
figfile2 = fullfile(pathname,'P13_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.5]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.5]);
saveas(gcf,figfile2,'pdf')
%
figure(2)
plot(t, patient,'LineWidth',0.5,'color',[37,37,37]/255);
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
scatter(T13,A13,6,'o','filled','MarkerEdgeColor','#CD950C','MarkerFaceColor','#CD950C','LineWidth',0.75,'Color','#CD950C');
hold off
ylim([1 5])
xlim([0 31])
xlabel('Time from symptom onset (days)')
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#CD950C';
box off
figfile2 = fullfile(pathname,'P13_ang');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.5]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.5]);
saveas(gcf,figfile2,'pdf')

%% P146

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420; %420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;

   p.d_V = 3.65;
        p.V0 =0.00000031;
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

patient = log10(curves);
p146= readtable('p146.csv');
T146= table2array(p146(:,4));
V146= table2array(p146(:,10));
R146= table2array(p146(:,11));
A146= table2array(p146(:,12));
figure(5)
plot(t, patient,'LineWidth',0.5,'color',[37,37,37]/255);
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
scatter(T146,R146,6,'o','filled','MarkerEdgeColor','#CD2626','MarkerFaceColor','#CD2626','LineWidth',0.75,'Color','#CD2626');
hold off
ylim([1.5 4.5])
xlim([0 31])
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='	#CD2626';
box off
figfile2 = fullfile(pathname,'P146_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.5]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.5]);
saveas(gcf,figfile2,'pdf')
%
figure(6)
plot(t, patient,'LineWidth',0.5,'color',[37,37,37]/255);
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
scatter(T146,A146,6,'o','filled','MarkerEdgeColor','#CD950C','MarkerFaceColor','#CD950C','LineWidth',0.75,'Color','#CD950C');
hold off
ylim([1.5 4.5])
xlim([0 31])
xlabel('Time from symptom onset (days)')
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#CD950C';
box off
figfile2 = fullfile(pathname,'P146_ang');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.5]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.5]);
saveas(gcf,figfile2,'pdf')
%% P90
pathname=fileparts('/Users/dengxiaoyan/Desktop/ViralLoadsFittingCode/Paper Code_reversion/Matlab code/Plot/');

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420; %420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;

   p.d_V = 0.52;
        p.V0 =0.00000071;
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

patient = log10(curves);
p90= readtable('p90.csv');
T90= table2array(p90(:,4));
V90= table2array(p90(:,10));
R90= table2array(p90(:,11));
A90= table2array(p90(:,12));
figure(3)
plot(t, patient,'LineWidth',0.5,'color',[37,37,37]/255);
hold on
ylim([0 3.5])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
scatter(T90,R90,6,'o','filled','MarkerEdgeColor','#CD2626','MarkerFaceColor','#CD2626','LineWidth',0.75,'Color','#CD2626');
hold off
ylim([0.5 4.5])
xlim([0 31])
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='	#CD2626';
box off
figfile2 = fullfile(pathname,'P90_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.5]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.5]);
saveas(gcf,figfile2,'pdf')
%
figure(4)
plot(t, patient,'LineWidth',0.5,'color',[37,37,37]/255);
hold on
ylim([0 3.5])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
scatter(T90,A90,6,'o','filled','MarkerEdgeColor','#CD950C','MarkerFaceColor','#CD950C','LineWidth',0.75,'Color','#CD950C');
hold off
ylim([0 4])
xlim([0 31])
xlabel('Time from symptom onset (days)')
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#CD950C';
box off
figfile2 = fullfile(pathname,'P90_ang');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.5]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.5]);
saveas(gcf,figfile2,'pdf')
