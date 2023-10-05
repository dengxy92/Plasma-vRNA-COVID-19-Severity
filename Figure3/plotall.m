clear all
input = readtable('all_ci.xlsx');
men= readtable('men_ci.xlsx');
women = readtable('women_ci.xlsx');
severe = readtable('severe_ci.xlsx');
moderate = readtable('moderate_ci.xlsx');
critical = readtable('critical_ci.xlsx');
alive= readtable('alive_ci.xlsx');
deceased = readtable('deceased_ci.xlsx');

dv = [2.01; table2array(input(:,3))];
V0 = [1.47e-6; table2array(input(:,4))];

dv_moderate = [4.5;table2array(moderate(:,3))];
V0_moderate =[2.43e-6;table2array(moderate(:,4))];

dv_severe = [2.96;table2array(severe(:,3))];
V0_severe =[4.18e-9;table2array(severe(:,4))]; %	0.00000000091

dv_critical = [1.12;table2array(critical(:,3))];
V0_critical =[1.82e-7;table2array(critical(:,4))];

dv_alive= [2.55;table2array(alive(:,3))];
V0_alive=[2.88e-6;table2array(alive(:,4))];

dv_deceased = [0.74;table2array(deceased(:,3))];
V0_deceased =[9.24e-10;table2array(deceased(:,4))];

dv_women = [2.48;table2array(women(:,3))];
V0_women =[2.43e-6;table2array(women(:,4))];


dv_men= [1.74;table2array(men(:,3))];
V0_men=[1.14e-6;table2array(men(:,4))];

all_data = readtable('all_mean.csv');
all_mean = table2array(all_data(:,3));
all_time = table2array(all_data(:,2));
all_std = table2array(all_data(:,4));
women_data = readtable('women_mean.csv');
women_mean = table2array(women_data(:,3));
women_time = table2array(women_data(:,2));
women_std = table2array(women_data(:,4));
men_data = readtable('men_mean.csv');
men_mean = table2array(men_data(:,3));
men_time = table2array(men_data(:,2));
men_std = table2array(men_data(:,4));
moderate_data = readtable('moderate_mean.csv');
moderate_mean = table2array(moderate_data(:,3));
moderate_time = table2array(moderate_data(:,2));
moderate_std = table2array(moderate_data(:,4));
severe_data = readtable('severe_mean.csv');
severe_mean = table2array(severe_data(:,3));
severe_time = table2array(severe_data(:,2));
severe_std = table2array(severe_data(:,4));
critical_data = readtable('critical_mean.csv');
critical_mean = table2array(critical_data(:,3));
critical_time = table2array(critical_data(:,2));
critical_std = table2array(critical_data(:,4));
alive_data = readtable('alive_mean.csv');
alive_mean = table2array(alive_data(:,3));
alive_time = table2array(alive_data(:,2));
alive_std = table2array(alive_data(:,4));
deceased_data = readtable('deceased_mean.csv');
deceased_mean = table2array(deceased_data(:,3));
deceased_time = table2array(deceased_data(:,2));
deceased_std = table2array(deceased_data(:,4));

%% All

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420; %420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;
   
patient = zeros(length(dv),1000);
for i = 1:length(dv)
        p.d_V = dv(i);
        p.V0 = V0(i);
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

        patient(i,:) = real(log10(curves));

end

patient_CI = zeros(1000,2);

for i = 1:1000
   patient_CI(i,:) =prctile(patient(2:101,i),[2.5,97.5]); % Calculates the 95% CB at every simulated point
end

%% women

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;
   
patient_women = zeros(length(dv_women)+1,1000);
for i = 1:length(dv_women)
        p.d_V = dv_women(i);
        p.V0 = V0_women(i);
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

        patient_women(i,:) = curves;

end

%patient(patient(:,1)>1,:) = [];
women_ci = zeros(1000,2);
for i = 1:1000
   women_ci(i,:) =prctile(real(log10(patient_women(2:size(patient_women,1),i))),[2.5,97.5]); % Calculates the 95% CB at every simulated point
end


%% Men

    p.T0 = 1.27;
    p.p = 420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;
   
patient_men = zeros(length(dv_men)+1,1000);
for i = 1:length(dv_men)
        p.d_V = dv_men(i);
        p.V0 = V0_men(i);
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

        patient_men(i,:) = curves;

end


men_ci = zeros(1000,2);
for i = 1:1000
   men_ci(i,:) =prctile(real(log10( patient_men(2:101,i))),[2.5,97.5]); % Calculates the 95% CB at every simulated point
end
%% Moderate

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;
   
patient_moderate = zeros(length(dv_moderate)+1,1000);
for i = 1:length(dv_moderate)
        p.d_V = dv_moderate(i);
        p.V0 = V0_moderate(i);
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

        patient_moderate(i,:) = curves;

end

%patient(patient(:,1)>1,:) = [];
moderate_ci = zeros(1000,2);
for i = 1:1000
   moderate_ci(i,:) =prctile(real(log10(patient_moderate(2:size(patient_moderate,1),i))),[2.5,97.5]); % Calculates the 95% CB at every simulated point
end


%% Severe
    p.T0 = 1.27;
    p.p = 420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;
   
patient_severe = zeros(length(dv_severe)+1,1000);
for i = 1:length(dv_severe)
        p.d_V = dv_severe(i);
        p.V0 = V0_severe(i);
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

        patient_severe(i,:) = curves;

end

severe_ci = zeros(1000,2);
for i = 1:1000
   severe_ci(i,:) =prctile(real(log10(patient_severe(2:length(dv_severe ),i))),[2.5,97.5]); % Calculates the 95% CB at every simulated point
end



%% critical

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;
   
patient_critical = zeros(length(dv_critical)+1,1000);
for i = 1:length(dv_critical)
        p.d_V = dv_critical(i);
        p.V0 = V0_critical(i);
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

        patient_critical(i,:) = curves;

end

%patient(patient(:,1)>1,:) = [];
critical_ci = zeros(1000,2);
for i = 1:1000
   critical_ci(i,:) =prctile(real(log10(patient_critical(2:size(patient_critical,1),i))),[2.5,97.5]); % Calculates the 95% CB at every simulated point
end



%% deceased

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;
   
patient_deceased = zeros(length(dv_deceased)+1,1000);
for i = 1:length(dv_deceased)
        p.d_V = dv_deceased(i);
        p.V0 = V0_deceased(i);
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

        patient_deceased(i,:) = curves;

end

%patient(patient(:,1)>1,:) = [];
deceased_ci = zeros(1000,2);
for i = 1:1000
   deceased_ci(i,:) =prctile(real(log10(patient_deceased(2:size(patient_deceased,1),i))),[2.5,97.5]); % Calculates the 95% CB at every simulated point
end


%% alive

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;
   
patient_alive = zeros(length(dv_alive)+1,1000);
for i = 1:length(dv_alive)
        p.d_V = dv_alive(i);
        p.V0 = V0_alive(i);
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

        patient_alive(i,:) = curves;

end


alive_ci = zeros(1000,2);
for i = 1:1000
   alive_ci(i,:) =prctile(real(log10( patient_alive(2:101,i))),[2.5,97.5]); % Calculates the 95% CB at every simulated point
end


%%
figure(1)
h1 = patch([t, fliplr(t)],[patient_CI(:,1)', fliplr(patient_CI(:,2)')],1,'facecolor','#E5E7EC','edgecolor','none'); %CI
hold on
%h2 = scatter(all_time,all_mean,6,'o','MarkerEdgeColor','#001253','LineWidth',0.75,'Color','#001253');
hold on
h3 = plot(t,patient(1,:),'LineWidth',0.5,'color',[37,37,37]/255);
hold on
errorbar(all_time,all_mean,all_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#001253');
hold on
yline(log10(68),'color','#878787','LineWidth',0.25)
hold on
yline(log10(13),'color','#878787','LineWidth',0.25,'LineStyle', '--')
%hold on
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 5])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.7]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.7]);
saveas(gcf,'all_mean','pdf')
%%
figure(2)
patch([t, fliplr(t)],[women_ci(:,1)', fliplr(women_ci(:,2)')],1,'facecolor','#E5E7EC','edgecolor','none'); %CI
hold on
errorbar(women_time,women_mean,women_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#001253');
hold on
plot(t,log10(patient_women(1,:)),'LineWidth',0.5,'color',[37,37,37]/255)
hold on
yline(log10(68),'color','#878787','LineWidth',0.25)
hold on
yline(log10(13),'color','#878787','LineWidth',0.25,'LineStyle', '--')
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 5])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.7]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.7]);
saveas(gcf,'women_mean','pdf')
%%
figure(3)
patch([t, fliplr(t)],[men_ci(:,1)', fliplr(men_ci(:,2)')],1,'facecolor','#E5E7EC','edgecolor','none'); %CI
hold on
errorbar(men_time,men_mean,men_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#001253');
hold on
plot(t,log10(patient_men(1,:)),'LineWidth',0.5,'color',[37,37,37]/255)
hold on
yline(log10(68),'color','#878787','LineWidth',0.25)
hold on
yline(log10(13),'color','#878787','LineWidth',0.25,'LineStyle', '--')
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 5])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.7]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.7]);
saveas(gcf,'men_mean','pdf')
%%
figure(4)
patch([t, fliplr(t)],[moderate_ci(:,1)', fliplr(moderate_ci(:,2)')],1,'facecolor','#E5E7EC','edgecolor','none'); %CI
hold on
errorbar(moderate_time,moderate_mean,moderate_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#001253');
hold on
plot(t,log10(patient_moderate(1,:)),'LineWidth',0.5,'color',[37,37,37]/255)
hold on
yline(log10(68),'color','#878787','LineWidth',0.25)
hold on
yline(log10(13),'color','#878787','LineWidth',0.25,'LineStyle', '--')
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 5])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 5.2 3.7]);
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.7]);
saveas(gcf,'moderate_mean','pdf')
%%
figure(5)
patch([t, fliplr(t)],[severe_ci(:,1)', fliplr(severe_ci(:,2)')],1,'facecolor','#E5E7EC','edgecolor','none'); %CI
hold on
errorbar(severe_time,severe_mean,severe_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#001253');
hold on
plot(t,log10(patient_severe(1,:)),'LineWidth',0.5,'color',[37,37,37]/255)
hold on
yline(log10(68),'color','#878787','LineWidth',0.25)
hold on
yline(log10(13),'color','#878787','LineWidth',0.25,'LineStyle', '--')
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 5])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.7]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.7]);
saveas(gcf,'severe_mean','pdf')
%%
figure(6)
patch([t, fliplr(t)],[critical_ci(:,1)', fliplr(critical_ci(:,2)')],1,'facecolor','#E5E7EC','edgecolor','none'); %CI
hold on
errorbar(critical_time,critical_mean,critical_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#001253');
hold on
plot(t,log10(patient_critical(1,:)),'LineWidth',0.5,'color',[37,37,37]/255)
hold on
yline(log10(68),'color','#878787','LineWidth',0.25)
hold on
yline(log10(13),'color','#878787','LineWidth',0.25,'LineStyle', '--')
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 5])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.7]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.7]);
saveas(gcf,'critical_mean','pdf')
%%
figure(7)
patch([t, fliplr(t)],[alive_ci(:,1)', fliplr(alive_ci(:,2)')],1,'facecolor','#E5E7EC','edgecolor','none'); %CI
hold on
errorbar(alive_time,alive_mean,alive_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#001253');
hold on
plot(t,log10(patient_alive(1,:)),'LineWidth',0.5,'color',[37,37,37]/255)
hold on
yline(log10(68),'color','#878787','LineWidth',0.25)
hold on
yline(log10(13),'color','#878787','LineWidth',0.25,'LineStyle', '--')
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 5])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters');
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.7]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.7]);
saveas(gcf,'alive_mean','pdf')
%%
figure(8)
patch([t, fliplr(t)],[deceased_ci(:,1)', fliplr(deceased_ci(:,2)')],1,'facecolor','#E5E7EC','edgecolor','none'); %CI
hold on
errorbar(deceased_time,deceased_mean,deceased_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#001253');
hold on
plot(t,log10(patient_deceased(1,:)),'LineWidth',0.5,'color',[37,37,37]/255)
hold on
yline(log10(68),'color','#878787','LineWidth',0.25)
hold on
yline(log10(13),'color','#878787','LineWidth',0.25,'LineStyle', '--')
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 5])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters');
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 5.2 3.7]); 
set(gca,'FontSize',5)
set(gcf, 'PaperSize', [5.2 3.7]);
saveas(gcf,'deceased','pdf')


