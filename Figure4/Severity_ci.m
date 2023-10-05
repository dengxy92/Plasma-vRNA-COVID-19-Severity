clear all

severe = readtable('severe_ci.xlsx');
moderate = readtable('moderate_ci.xlsx');
critical = readtable('critical_ci.xlsx');



dv_moderate = [4.5;table2array(moderate(:,3))];
V0_moderate =[0.0000214;table2array(moderate(:,4))];

dv_severe = [2.96;table2array(severe(:,3))];
V0_severe =[4.18e-9;table2array(severe(:,4))]; %	0.00000000091

dv_critical = [1.12;table2array(critical(:,3))];
V0_critical =[1.82e-7;table2array(critical(:,4))];



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


%dv_severe = [3.43;table2array(severe(:,3))];
%V0_severe =[0.0013;table2array(severe(:,4))];



%%

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

%
moderate_data = readtable('moderate_mean.csv');
moderate_time = table2array(moderate_data(:,2));
moderate_mean = table2array(moderate_data(:,3));
moderate_std = table2array(moderate_data(:,4));

figure
patch([t, fliplr(t)],[moderate_ci(:,1)', fliplr(moderate_ci(:,2)')],1,'facecolor','#001253','edgecolor','none','facealpha', 0.1); %CI
hold on
errorbar(moderate_time,moderate_mean,moderate_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',8,'LineWidth',1.0,'Color','#001253');
hold on
plot(t,log10(patient_moderate(1,:)),'LineWidth',1,'color',[37,37,37]/255)
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 4])
xlim([0 30])
set(gcf, 'PaperPosition', [0 0 7 4.5]); 
set(gcf, 'PaperSize', [7 4.5]);


%%

% Initial condition and parameter guesses --------------------------------

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
severe_data = readtable('severe_mean.csv');
severe_mean = table2array(severe_data(:,3));
severe_time = table2array(severe_data(:,2));
severe_std = table2array(severe_data(:,4));
%
figure
patch([t, fliplr(t)],[severe_ci(:,1)', fliplr(severe_ci(:,2)')],1,'facecolor','#001253','edgecolor','none','facealpha', 0.1); %CI
hold on
errorbar(severe_time,severe_mean,severe_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',8,'LineWidth',1.0,'Color','#001253');
hold on
plot(t,log10(patient_severe(1,:)),'LineWidth',1,'color',[37,37,37]/255)
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 4])
xlim([0 30])
set(gcf, 'PaperPosition', [0 0 7 4.5]); 
set(gcf, 'PaperSize', [7 4.5]);


%%

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

%
critical_data = readtable('critical_mean.csv');
critical_time = table2array(critical_data(:,2));
critical_mean = table2array(critical_data(:,3));
critical_std = table2array(critical_data(:,4));

figure
patch([t, fliplr(t)],[critical_ci(:,1)', fliplr(critical_ci(:,2)')],1,'facecolor','#001253','edgecolor','none','facealpha', 0.1); %CI
hold on
errorbar(critical_time,critical_mean,critical_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',8,'LineWidth',1.0,'Color','#001253');
hold on
plot(t,log10(patient_critical(1,:)),'LineWidth',1,'color',[37,37,37]/255)
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 4])
xlim([0 30])
set(gcf, 'PaperPosition', [0 0 7 4.5]); 
set(gcf, 'PaperSize', [7 4.5]);

%%
figure(4)
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
h3 = patch([t, fliplr(t)],[moderate_ci(:,1)', fliplr(moderate_ci(:,2)')],1,'facecolor','#1874CD','edgecolor','none','facealpha', 0.1); %CI
hold on
h2 = plot(t,log10(patient_moderate(1,:)),'LineWidth',0.5,'color',[16,78,139]/255);
hold on
h6 = patch([t, fliplr(t)],[severe_ci(:,1)', fliplr(severe_ci(:,2)')],1,'facecolor','#CD9B1D','edgecolor','none','facealpha', 0.1); %CI
hold on
h5 = plot(t,log10(patient_severe(1,:)),'LineWidth',0.5,'color',[139,117,0]/255);
hold on
h9 = patch([t, fliplr(t)],[critical_ci(:,1)', fliplr(critical_ci(:,2)')],1,'facecolor','#B22222','edgecolor','none','facealpha', 0.1); %CI
hold on
h8 = plot(t,log10(patient_critical(1,:)),'LineWidth',0.5,'color',[178,34,34]/255);
hold on
h4 = scatter(severe_time,severe_mean,6,'o','MarkerEdgeColor',[139,117,0]/255,'LineWidth',0.75,'Color',[139,117,0]/255);
hold on
h1 = scatter(moderate_time,moderate_mean,6,'o','MarkerEdgeColor','#104E8B','LineWidth',0.75,'Color','#104E8B');
hold on
h7 = scatter(critical_time,critical_mean,6,'o','MarkerEdgeColor','#B22222','LineWidth',0.75,'Color','#B22222');
hold off
xlim([0 30])
ylim([0 3])
%legend([h1,h4,h7,h2,h5,h8],'FontSize',14,'NumColumns',2)
%legend boxoff
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 7 4.5]); 
set(gca,'FontSize',6)
set(gcf, 'PaperSize', [7 4.5]);
saveas(gcf,'Severity_mean','pdf')
