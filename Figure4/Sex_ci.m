clear all

men= readtable('men_ci.xlsx');
women = readtable('women_ci.xlsx');

dv_women = [2.48;table2array(women(:,3))];
V0_women =[2.43e-6;table2array(women(:,4))];

dv_men= [1.74;table2array(men(:,3))];
V0_men=[1.14e-6;table2array(men(:,4))];

women_data = readtable('women_mean.csv');
women_mean = table2array(women_data(:,3));
women_time = table2array(women_data(:,2));
women_std = table2array(women_data(:,4));

men_data = readtable('men_mean.csv');
men_mean = table2array(men_data(:,3));
men_time = table2array(men_data(:,2));
men_std = table2array(men_data(:,4));

%%

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


figure
patch([t, fliplr(t)],[women_ci(:,1)', fliplr(women_ci(:,2)')],1,'facecolor','#001253','edgecolor','none','facealpha', 0.1); %CI
hold on
%errorbar(women_time,women_mean,women_std,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',8,'LineWidth',1.0,'Color','#001253');
hold on
plot(t,log10(patient_women(1,:)),'LineWidth',0.5,'color',[37,37,37]/255)
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 3])
xlim([0 30])
set(gcf, 'PaperPosition', [0 0 8 5]); 
set(gcf, 'PaperSize', [8 5]);


%%

% Initial condition and parameter guesses --------------------------------

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

%
figure
patch([t, fliplr(t)],[men_ci(:,1)', fliplr(men_ci(:,2)')],1,'facecolor','#001253','edgecolor','none','facealpha', 0.1); %CI
hold on
plot(t,log10(patient_men(1,:)),'LineWidth',0.5,'color',[37,37,37]/255)
%scatter(DSO, Observation,'filled','MarkerFaceColor',[99,99,99]/255)
hold off
xlim([0 30])
ylim([0 3])
xlim([0 30])
set(gcf, 'PaperPosition', [0 0 8 5]); 
set(gcf, 'PaperSize', [8 5]);

%%

men_data = readtable('men_mean.csv');
men_mean = table2array(men_data(:,3));
men_time = table2array(men_data(:,2));
men_std = table2array(men_data(:,4));
women_data = readtable('women_mean.csv');
women_mean = table2array(women_data(:,3));
women_time = table2array(women_data(:,2));
women_std = table2array(women_data(:,4));
%%
figure(3)
h1 = patch([t, fliplr(t)],[women_ci(:,1)', fliplr(women_ci(:,2)')],1,'facecolor','#B22222','edgecolor','none','facealpha', 0.1); %CI
hold on
h2 = plot(t,log10(patient_women(1,:)),'LineWidth',0.5,'color',[178,34,34]/255);
hold on
h3 = scatter(women_time,women_mean,6,'o','MarkerEdgeColor','#B22222','LineWidth',0.75,'Color','#B22222');
hold on
h4 = patch([t, fliplr(t)],[men_ci(:,1)', fliplr(men_ci(:,2)')],1,'facecolor','#1874CD','edgecolor','none','facealpha', 0.1); %CI
hold on
h5 = plot(t,log10(patient_men(1,:)),'LineWidth',0.5,'color',[16,78,139]/255);
hold on
h6 = scatter(men_time,men_mean,6,'o','MarkerEdgeColor','#104E8B','LineWidth',0.75,'Color','#104E8B');
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
set(gcf, 'PaperPosition', [0 0 7 4.5]); 
set(gca,'FontSize',6)
set(gcf, 'PaperSize', [7 4.5]);
saveas(gcf,'Sex','pdf')
