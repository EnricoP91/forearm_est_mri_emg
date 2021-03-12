function [Vout, Iest, Vfullgraph] = estimateCurrents(LFM, rms_graph,rms_armcirc,electrode_index, selected_el, LR, estID, muscleareas_norm, ifplot)

% *************************************************************************
% This function calculate the estimated currents of the RMS of the input
% signal once the LFM of the system is known, using different estimation 
% methods inclusing L2 regularized minimum norm, L2 regularizion with a 
% penalyzing term proportional to the inverse of the area, L1 regularization
% @param
%       LFM
%       rms_graph
%       rms_armcirc
%       electrode index
%       selected_el
%       LR
%       estimationID
%       areavec
%       ifplot
%       
% @ out
%       Vout
%       Iest
% 
% *************************************************************************

[nelectrode, nmuscle] = size(LFM);
nnode = nelectrode + nmuscle;

switch estID 
    case 'L2reg'
        % 2] L2- REGULARIZED MINIMUM-NORM SOLUTION  ===============================
        % % 2.1] squared residuals vector
        title_L2R = 'Estimated Potentials (L2-REGULARIZED MINIMUM NORM)';
        
        RSS_vec = [];
        indRSS = 1 ;
        lambdavec = [];
        for lambda = 10 : 10 : 200
            lambdavec = [lambdavec, lambda];
   
            % current estimation 
            Iest_muscle_L2MN = inv(LFM'*LFM + lambda*eye(nmuscle))*(LFM'*rms_graph);
           
            % reprojection of the currents on the electrode space
            [Vout_tmp, ~, ~] = plot_elpot(Iest_muscle_L2MN, nnode, electrode_index, LR, selected_el, rms_armcirc, 0, title_L2R, '','LRReg', 0);
            
            % Residual calculation
            RSS_vec(indRSS) = sqrt(sum((rms_armcirc - Vout_tmp).^2)); 
            indRSS = indRSS+1;
        
        end
        
        % find the lambda that guarantee the minimum RSS
        ind = (RSS_vec == min(RSS_vec));
        lambda_min = lambdavec(ind);
        
        % Re-estimate the minimum current
        Iest_muscle = inv(LFM'*LFM + lambda_min*eye(nmuscle))*(LFM'*VsEMG_graph);
        
        % reprojection on the electrode space
        [Vout, Vfullgraph, ~] = plot_elpot(Iest_muscle, nnode, electrode_index, LR, selected_el, Vreal_arm, ifplot, title_L2R,'','L2Reg',0);
        
    case 'L2reg_area'
        if nargin > 8
            % 2.2]  penalty function based on the inverse of the area
            title_L2RAreaPen = 'Estimated Potentials (L2-REGULARIZED MINIMUM NORM SOLUTION with AREA PENALTY)';

            % selection of muscle area for penalty term 
            lambda = (1./muscleareas_norm);
            
            % estimation of the currents
            Iest = inv(LFM'*LFM + diag(lambda))*(LFM'*rms_graph);
            
            % estimation of the 
            [Vout, Vfullgraph, ~] = plot_elpot(Iest, nnode, electrode_index, LR, selected_el, rms_armcirc, ifplot, title_L2RAreaPen, '','L2Reg_Area',0);
        else
            disp('wrong number of input parameter')
            return
        end
        
    case 'L1reg'
        % 4] L1- REGULARIZED MINIMUM NORM ========================================
%         title_L1 = 'L1-REGULARIZED LEAST SQUARE';
%         
%         % i) fixed LAMBDA
%         lambda_L1 = 0.01;
%         [Iest, ~] = l1_ls(LFM, rms_graph, lambda_L1);
%         [Vout, Vfullgraph, ~] = plot_elpot(Iest, nnode, electrode_index, LR, selected_el, rms_armcirc, ifplot, title_L1, '','L1Reg',0);

    case 'L1reg_2'
        % 1) variable LAMBDA
        %     
        % for lambda_L1 = 0.01:0.01:1
        %     
        %         [Iest_muscle_L1reg(:,testID), status] = l1_ls(LFM, VsEMG_graph(:,testID), lambda_L1);
        % 
        %         ifplot = 0;
        %         [Vout_arm_L1Reg,h5] = plot_elpot(Iest_muscle_L1reg(:,testID), nnode, electrode_index, LR, selected_el, Vreal_arm(:,testID), ifplot, title_L1, savepath,'L1Reg',savePotPlot);
        %         Vest_L1Reg(:,testID) = Vout_arm_L1Reg;
        %         
        %     end
        
    case 'L1_3'
%         % 4.2) with non negative constraint
%         title('L1 regularization with nonnegativity constraint')
%         lambda_L1 = 0.01;
%         [Iest, ~] = l1_ls_nonneg(LFM, rms_graph, lambda_L1);
%        
%         [Vout, Vfullgraph, ~] = plot_elpot(Iest, nnode, electrode_index, LR, selected_el, rms_armcirc, 0, title_L1, '','L1Reg',0);
        
        
    otherwise 
        error("The estimation %s is not existing.", estID)
end