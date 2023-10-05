clear all
par1 = readtable("elimination_new.xlsx");

all_senior = readtable('all_senior_log10.csv');
moderate_senior =readtable('moderate_senior_log10.csv');
severe_senior =  readtable('severe_senior_log10.csv');
critical_senior =  readtable('critical_senior_log10.csv');
women_senior=readtable('women_senior_log10.csv');
men_senior= readtable('men_senior_log10.csv');
alive_senior = readtable('alive_senior_log10.csv');
deceased_senior = readtable('deceased_senior_log10.csv');

all_young = readtable('all_young_log10.csv');
moderate_young =readtable('moderate_young_log10.csv');
severe_young =  readtable('severe_young_log10.csv');
critical_young =  readtable('critical_young_log10.csv');
women_young=readtable('women_young_log10.csv');
men_young= readtable('men_young_log10.csv');
alive_young = readtable('alive_young_log10.csv');
deceased_young = readtable('deceased_young_log10.csv');

% Simulate viral loads and calculate mean and SE of dynamics

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;
    %p.V0 = 0.00000023;

%all
%all
dv_patient = table2array(par1(1:16,2));
V0_patient = table2array(par1(1:16,3));
%p.V0 = 0.00603484006202506;
solutions_patient = zeros(16,1000);
for i = 1:16
        p.d_V = dv_patient(i);
        p.V0 = V0_patient(i);
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

        solutions_patient(i,:) = real(log10(curves));

end
%%

figure(1)

subplot(4,2,1)
scatter(table2array(all_young(:,2)),table2array(all_young(:,3)),8,'filled','MarkerFaceColor','#1874CD');
hold on
scatter(table2array(all_senior(:,2)),table2array(all_senior(:,3)),8,'filled','MarkerFaceColor','#B22222');
hold on
plot(t,solutions_patient(1,:),'color','#1874CD');
hold on
plot(t,solutions_patient(9,:),'color','#B22222');
hold off
xlim([0 30])
ylim([0 3])
%legend('Young data',' Senior data','Pred young','Pred senior')
%title('All patients')
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gca,'FontSize',6);

subplot(4,2,2)
h1 = scatter(table2array(women_young(:,2)),table2array(women_young(:,3)),8,'filled','MarkerFaceColor','#1874CD');
hold on
h2 = scatter(table2array(women_senior(:,2)),table2array(women_senior(:,3)),8,'filled','MarkerFaceColor','#B22222');
hold on
h3 = plot(t,solutions_patient(5,:),'color','#1874CD');
hold on
h4 = plot(t,solutions_patient(13,:),'color','#B22222');
hold off
xlim([0 30])
ylim([0 3])
%legend('Young data',' Senior data','Pred young','Pred senior')
%title('Women')
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gca,'FontSize',6);

subplot(4,2,3)
h1 = scatter(table2array(men_young(:,2)),table2array(men_young(:,3)),8,'filled','MarkerFaceColor','#1874CD');
hold on
h2 = scatter(table2array(men_senior(:,2)),table2array(men_senior(:,3)),8,'filled','MarkerFaceColor','#B22222');
hold on
h3 = plot(t,solutions_patient(6,:),'color','#1874CD');
hold on
h4 = plot(t,solutions_patient(14,:),'color','#B22222');
hold off
xlim([0 30])
ylim([0 3])
%legend('Young data',' Senior data','Pred young','Pred senior')
%title('Men')
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gca,'FontSize',6);

subplot(4,2,4)
h1 = scatter(table2array(moderate_young(:,2)),table2array(moderate_young(:,3)),8,'filled','MarkerFaceColor','#1874CD');
hold on
h2 = scatter(table2array(moderate_senior(:,2)),table2array(moderate_senior(:,3)),8,'filled','MarkerFaceColor','#B22222');
hold on
h3 = plot(t,solutions_patient(2,:),'color','#1874CD');
hold on
h4 = plot(t,solutions_patient(10,:),'color','#B22222');
hold off
xlim([0 30])
ylim([0 3])
%legend('Young data',' Senior data','Pred young','Pred senior')
%title('Moderate')
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gca,'FontSize',6);

subplot(4,2,5)
h1 = scatter(table2array(severe_young(:,2)),table2array(severe_young(:,3)),8,'filled','MarkerFaceColor','#1874CD');
hold on
h2 = scatter(table2array(severe_senior(:,2)),table2array(severe_senior(:,3)),8,'filled','MarkerFaceColor','#B22222');
hold on
h3 = plot(t,solutions_patient(3,:),'color','#1874CD');
hold on
h4 = plot(t,solutions_patient(11,:),'color','#B22222');
hold off
xlim([0 30])
ylim([0 3])
%legend('Young data',' Senior data','Pred young','Pred senior')
%title('Severe')
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gca,'FontSize',6);

subplot(4,2,6)
h1 = scatter(table2array(critical_young(:,2)),table2array(critical_young(:,3)),8,'filled','MarkerFaceColor','#1874CD');
hold on
h2 = scatter(table2array(critical_senior(:,2)),table2array(critical_senior(:,3)),8,'filled','MarkerFaceColor','#B22222');
hold on
h3 = plot(t,solutions_patient(4,:),'color','#1874CD');
hold on
h4 = plot(t,solutions_patient(12,:),'color','#B22222');
hold off
xlim([0 30])
ylim([0 3])
%legend('Young data',' Senior data','Pred young','Pred senior')
%title('Critical')
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gca,'FontSize',6);


subplot(4,2,7)
h1 = scatter(table2array(alive_young(:,2)),table2array(alive_young(:,3)),8,'filled','MarkerFaceColor','#1874CD');
hold on
h2 = scatter(table2array(alive_senior(:,2)),table2array(alive_senior(:,3)),8,'filled','MarkerFaceColor','#B22222');
hold on
h3 = plot(t,solutions_patient(7,:),'color','#1874CD');
hold on
h4 = plot(t,solutions_patient(15,:),'color','#B22222');
hold off
xlim([0 30])
ylim([0 3])
%legend('Young data',' Senior data','Pred young','Pred senior')
%title('Alive')
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gca,'FontSize',6);

subplot(4,2,8)
h1 = scatter(table2array(deceased_young(:,2)),table2array(deceased_young(:,3)),8,'filled','MarkerFaceColor','#1874CD');
hold on
h2 = scatter(table2array(deceased_senior(:,2)),table2array(deceased_senior(:,3)),8,'filled','MarkerFaceColor','#B22222');
hold on
h3 = plot(t,solutions_patient(8,:),'color','#1874CD');
hold on
h4 = plot(t,solutions_patient(16,:),'color','#B22222');
hold off
xlim([0 30])
ylim([0 3])
%legend('Young data',' Senior data','Pred young','Pred senior', 'location','southeast')
%title('Deceased')
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gca,'FontSize',6);

figfile2 = fullfile(pathname,'age');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',6);
set(gcf, 'PaperPosition', [0 0 18 20]); 
set(gcf, 'PaperSize', [18 20]);
set(gca, 'LooseInset', get(gca,'TightInset'))
saveas(gcf,figfile2,'pdf')
