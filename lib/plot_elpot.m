function [VoutArm, Vout, h1] = plot_elpot(Imuscles, nnode, electrode_index, LR, selected_el, Vreal_arm, ifplot, title_string,savepath, name, ifsave)
    
    nelectrode = nnode - length(Imuscles);
    
    % reordering the current according to the graph order 
    Iest_fullnode = zeros(nnode,1);
    indm = 1 ;      % index on the current vector
    for nodeID = 1 : nnode
        if electrode_index(nodeID) == 0 
            Iest_fullnode(nodeID,1) = Imuscles(indm,1);
            indm = indm +1;
        end
    end
    
    % estimated potential on the electrode space
    Vout = inv(LR)*(Iest_fullnode.*1e3); % [mV]
%     Vout = inv(LR)*(Iest_fullnode); % [V]
       
    % reorder: graph --> sequential
    VoutArm = zeros(nelectrode,1);
    for row_ind = 1 : length(selected_el)
        el_ind = find(electrode_index == selected_el(row_ind));
        VoutArm(row_ind,1) = Vout(el_ind,1);
    end

    if ifplot == 1
        % Electrodes evenly distributed around 360 degrees
        theta = linspace(0, 2*pi, nelectrode+1);

        % Data estimated potential
        rho1 = VoutArm;
        rho1 = [rho1 ; rho1(1)];

         % Data measured potential
        rho2 = Vreal_arm;
        rho2 = [rho2 ; rho2(1)];


        h1 =  figure();
        polarplot(theta,rho1, 'b-o','Linewidth',1.5);
        hold on 
        polarplot(theta,rho2, 'r--','Linewidth',1.5);
        title(title_string)
        legend('Estimated electrode voltage', 'Measured EMG Voltage')
    else
        h1 = [];
    end
    
    if ifsave
        filename = strcat('\estimatedvoltages_', name);
        saveas(h1,[strcat(savepath, filename),'.fig'])
        saveas(h1,[strcat(savepath, filename),'.pdf'])
        saveas(h1,[strcat(savepath, filename),'.jpg'])
    end
    
end
