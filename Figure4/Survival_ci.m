clear all

alive= readtable('alive_ci.xlsx');
deceased = readtable('deceased_ci.xlsx');
dv_alive= [2.55;table2array(alive(:,3))];
V0_alive=[2.88e-6;table2array(alive(:,4))];

dv_deceased = [0.74;table2array(deceased(:,3))];
V0_deceased =[9.24e-10;table2array(deceased(:,4))];


alive_data = readtable('alive_mean.csv');
alive_mean = table2array(alive_data(:,3));
alive_time = table2array(alive_data(:,2));
alive_std = table2array(alive_data(:,4));

deceased_data = readtable('deceased_mean.csv');
deceased_mean = table2array(deceased_data(:,3));
deceased_time = table2array(deceased_data(:,2));
deceased_std = table2array(deceased_data(:,4));

%%

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

%
%deceased = readtable('deceased_censor_new.csv');
%deceased_mean = table2array(deceased(:,4));
%deceased_time = table2array(deceased(:,3));
%deceased_std = table2array(deceased(:,5));
%deceased_se= table2array(deceased(:,6));
figure
patch([t, fliplr(t)],[deceased_ci(:,1)', fliplr(deceased_ci(:,2)')],1,'facecolor','#001253','edgecolor','none','facealpha', 0.1); %CI
hold on
%errorbar(deceased_time,deceased_mean,deceased_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',8,'LineWidth',1.0,'Color','#001253');
hold on
plot(t,log10(patient_deceased(1,:)),'LineWidth',1,'color',[37,37,37]/255)
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 4])
xlim([0 30])
set(gcf, 'PaperPosition', [0 0 7 4]); 
set(gcf, 'PaperSize', [7 4]);


%%

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

%
figure
patch([t, fliplr(t)],[alive_ci(:,1)', fliplr(alive_ci(:,2)')],1,'facecolor','#001253','edgecolor','none','facealpha', 0.1); %CI
hold on
plot(t,log10(patient_alive(1,:)),'LineWidth',1,'color',[37,37,37]/255)
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 4])
xlim([0 30])
set(gcf, 'PaperPosition', [0 0 7 4]); 
set(gcf, 'PaperSize', [7 4]);


%%
figure(3)
alive_data = readtable('alive_mean.csv');
alive_mean = table2array(alive_data(:,3));
alive_time = table2array(alive_data(:,2));
alive_std = table2array(alive_data(:,4));
deceased_data = readtable('deceased_mean.csv');
deceased_mean = table2array(deceased_data(:,3));
deceased_time = table2array(deceased_data(:,2));
deceased_std = table2array(deceased_data(:,4));

h1 = patch([t, fliplr(t)],[deceased_ci(:,1)', fliplr(deceased_ci(:,2)')],1,'facecolor','#B22222','edgecolor','none','facealpha', 0.1); %CI
hold on
h2 = plot(t,log10(patient_deceased(1,:)),'LineWidth',0.5,'color',[178,34,34]/255);
hold on
h3 = scatter(deceased_time,deceased_mean,6,'o','MarkerEdgeColor','#B22222','LineWidth',0.75,'Color','#B22222');
hold on
h4 = patch([t, fliplr(t)],[alive_ci(:,1)', fliplr(alive_ci(:,2)')],1,'facecolor','#1874CD','edgecolor','none','facealpha', 0.1); %CI
hold on
h5 = plot(t,log10(patient_alive(1,:)),'LineWidth',0.5,'color',[16,78,139]/255);
hold on
h6 = scatter(alive_time,alive_mean,6,'o','MarkerEdgeColor','#104E8B','LineWidth',0.75,'Color','#104E8B');
hold on
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 3])
%legend( [h6, h3, h5,h2],'FontSize',14,'NumColumns',2)
%legend boxoff
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gca,'FontSize',6)
set(gcf, 'PaperPosition', [0 0 7 4.5]); 
set(gcf, 'PaperSize', [7 4.5]);
saveas(gcf,'Survival','pdf')
