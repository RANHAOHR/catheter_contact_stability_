function [velocity_samples] = blood_flow
    subplot(2,2,1);
    
    s = what('cwru_catheter_kinematics');
    s_data = strcat(s.path,'/data/pulmonary.mat');
    load(s_data);
    [~,idx] = unique(data(:,1));
    data = data(idx,:);
    x = data(:,1);
    pul_flow_curve = data(:,2); % cm^3/s
    
    diam = 0.027; % 0.027m
    cross_area = pi * diam^2 / 4; % m^2

    flow_curve = pul_flow_curve * 10^(-6) / cross_area;
    plot(x,flow_curve,'r','LineWidth',2);
    xlabel('t /s');
    ylabel('Pulmonary blood flow velocity m/s');    

    %interpolation
    sample_size = 10000;
    dmin = min(x);
    dmax = max(x);
    delta = (dmax - dmin) / (sample_size - 1);
    query = [dmin:delta:dmax];
    flow_curve = interp1(x, flow_curve, query);

    %turn this distribution into a pdf
    positive_curve = abs(flow_curve);
% 

    area = sum( delta * positive_curve );
    pdf_curve = positive_curve / area;
%         figure(2);
%     plot(query,pdf_curve,'r','LineWidth',2);

    %inverse transform sampling
    cumulative_curve = zeros(1, sample_size);
    for i = 1: sample_size
        cumulative_curve(i) = delta * sum(pdf_curve(1:i));
    end
    % plot(x,cumulative_curve,'r','LineWidth',2); %check this

    %Inverse transform sampling
    rng('default'); 
    velocity_samples = zeros(1,sample_size);
    for i = 1: sample_size
        u = rand;
        velocity_samples(i) = flow_curve(find(cumulative_curve < u, 1, 'last')); %use the original curve for negtive velocity
    end
    subplot(2,2,2);
%     for i = 1: sample_size
%         line([velocity_samples(i),velocity_samples(i)],[0,0.5]);
%         hold on;
%     end
%     dmin = min(velocity_samples) -0.1;
%     dmax = max(velocity_samples) + 0.1;
%     set(gca,'ytick',[]) %need not show the y axis.
%     axis([dmin dmax 0 1]);
    histogram(velocity_samples);

    xlabel('Distribution of pulmonary blood velocity samples $||v||$ m/s','Interpreter','latex');
    Setfontsize = 15;
    set(get(gca,'Xlabel'),'FontSize',Setfontsize);
    set(get(gca,'Ylabel'),'FontSize',Setfontsize);
end

 