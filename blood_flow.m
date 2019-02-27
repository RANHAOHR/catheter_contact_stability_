function [velocity_samples] = blood_flow(selection)
    figure(1);
    selection = 2;
    if selection ==1
        load  /home/ranhao/Downloads/catheter_project/possible_data/AORTA_ROOT_physio.mat
        Npatient = [1,4,10,14,16]; %average of mutiple patients
        Nsample = 956;
        flow_curve = zeros(1, Nsample);

        diam = 0.0303; % m
        cross_area = pi * diam^2 / 4; % m^2

        for i = 1: size(Npatient,2)
            sample = AORTA_ROOT_PHYSIO(Npatient(i)); % AORTA_ROOT_PHYSIO
            sample = sample{1};
            one_cycle = sample.ONE_CYCLE;
            flow_curve = flow_curve + one_cycle(:,3)';  % [one_cycle(:,1);one_cycle(:,2);one_cycle(:,3);one_cycle(:,4)];
        end

        flow_curve = flow_curve / size(Npatient,2); % m^3/s
        flow_curve = flow_curve / cross_area; % velocity m/s

        x = [1:Nsample];
        x  = x / Nsample;
        plot(x,flow_curve,'r','LineWidth',2);
        xlim([0 max(x)])
        xlabel('t /s');
        ylabel('Aortic blood flow velocity m/s');
    elseif selection == 2
        load  /home/ranhao/Downloads/catheter_project/possible_data/pulmonary.mat
        [~,idx] = unique(data(:,1));
        data = data(idx,:);
        x = data(:,1);
        pul_flow_curve = data(:,2);
        Nsample = size(x,2);

        diam = 0.027; % m
        cross_area = pi * diam^2 / 4; % m^2

        flow_curve = pul_flow_curve * 10^(-6) / cross_area;
        plot(x,flow_curve,'r','LineWidth',2);
        xlabel('t /s');
        ylabel('Pulmonary blood flow velocity m/s');    
    end

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
%     figure(3);
%     subplot(2,1,1);
% %     for i = 1: sample_size
% %         line([velocity_samples(i),velocity_samples(i)],[0,0.5]);
% %         hold on;
% %     end
% histogram(velocity_samples);
%     dmin = min(velocity_samples) -0.1;
%     dmax = max(velocity_samples) + 0.1;
% %     set(gca,'ytick',[]) %need not show the y axis.
% %     axis([dmin dmax 0 1]);
%     xlabel('Distribution of pulmonary blood velocity samples $||v||$ m/s','Interpreter','latex');
% Setfontsize = 15;
% set(get(gca,'Xlabel'),'FontSize',Setfontsize);
% set(get(gca,'Ylabel'),'FontSize',Setfontsize);
end

 