clear all
data = table2array(readtable('average.csv'));

% Simulate viral loads and calculate mean and SD of dynamics

% Initial condition and parameter guesses --------------------------------

    p.T0 = 0.16;
    p.p = 741.2;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.144*0.5; 
    p.t_inf = 0;
    p.bet = 0.289;
    p.V0 =4.5; %4.5;%12
    p.d_V = 14.5;
       
    
        p.IC = [p.T0,p.I0,p.V0];
    
       
        
        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,30]);
        
        t = linspace(0,25,1000);
        curves = deval(sol,t,3);

%figure(1)
%plot(data(:,1),data(:,2),'.','LineWidth', 2)
%hold on 
%plot(t,curves,'LineWidth',2)
%hold off
%ylabel('(copies/mL)')
%title('Average')
%set(gcf, 'PaperPositionMode', 'manual'); 
%set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPosition',[0 0 10 8]); 

%saveas(gcf,'average','epsc')
%saveas(gcf,'average','fig')

%
P1= table2array(readtable('P1.csv'));
P2= table2array(readtable('P2.csv'));
P3= table2array(readtable('P3.csv'));
P4= table2array(readtable('P4.csv'));
P5= table2array(readtable('P5.csv'));
P6= table2array(readtable('P6.csv'));
P7= table2array(readtable('P7.csv'));
P8= table2array(readtable('P8.csv'));
color = [140,81,10;
191,129,45;
223,194,125;
246,232,195;
199,234,229;
128,205,193;
53,151,143;
1,102,94]./255;

figure
plot(t,curves,'LineWidth',2,'Color','k')
hold on
h1 = plot(P1(:,1),P1(:,2),'o','Color',color(1,:),'LineWidth', 2);
hold  on
h2 = plot(P2(:,1),P2(:,2),'o','Color',color(2,:),'LineWidth', 2);
hold  on
h3 = plot(P3(:,1),P3(:,2),'o','Color',color(3,:),'LineWidth', 2);
hold  on
h4 = plot(P4(:,1),P4(:,2),'o','Color',color(4,:),'LineWidth', 2);
hold  on
h5 = plot(P5(:,1),P5(:,2),'o','Color',color(5,:),'LineWidth', 2);
hold  on
h6 = plot(P6(:,1),P6(:,2),'o','Color',color(6,:),'LineWidth', 2);
hold  on
h7 = plot(P7(:,1),P7(:,2),'o','Color',color(7,:),'LineWidth', 2);
hold  on
h8 = plot(P8(:,1),P8(:,2),'o','Color',color(8,:),'LineWidth', 2);
hold off
legend([h1, h2, h3 ,h4, h5,h6, h7, h8],{'P1','P2','P3','P4','P5','P6','P7','P8'},'NumColumns',2)
xlim([0 25])
ylim([0 10])
ylabel('Nasopharyngeal viral load (log_{10}(copies/mL))')
xlabel('Time from symptom onset (days)')
set(gca,'FontSize',15)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition',[0 0 20 14]); 
set(gcf, 'PaperSize', [20 14]);
saveas(gcf,'Goyal','pdf')
saveas(gcf,'Goyal','fig')

%% Simulate viral loads and calculate mean and SD of dynamics

% Initial condition and parameter guesses --------------------------------

    p.T0 = 0.16;
    p.p = 741.2;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.144; 
    p.t_inf = 0;
    p.bet = 0.289;
    p.V0 =4.5; %4.5;%12
    p.d_V = 10;
       
    
        p.IC = [p.T0,p.I0,p.V0];
    
       
        
        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,30]);
        
        t = linspace(0,30,1000);
        curves = deval(sol,t,3);

        solutions_patient = curves;

figure
plot(P1(:,1),P1(:,2),'.','LineWidth', 2)
hold on
plot(t,curves,'LineWidth',2)
hold off
title('Patient 1')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition',[0 0 10 8]); 

saveas(gcf,'P1','epsc')
saveas(gcf,'P1','fig')

figure
plot(P1(:,1),log10(P1(:,2)),'.','LineWidth', 2)
hold on
plot(t,log10(curves),'LineWidth',2)
hold off
ylabel('log_{10} (copies/mL)')
title('Patient 1')
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition',[0 0 10 8]); 

saveas(gcf,'P1_log10','epsc')
saveas(gcf,'P1_log10','fig')


   
