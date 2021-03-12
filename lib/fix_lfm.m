function [LFM,LFM_full] = fix_lfm(LR, electrode_index) 
    
    MIN_Z_EIGENVALUE = 1e-3;
    MIN_CONDITION_NUMBER = 1e3;
        
    % Fix the LFM matrix if the conditioning number is too high
    nnode = size(LR,1);
    nelectrode = nnz(electrode_index);
    nmuscle = nnode-nelectrode;
    
    % Recalculate the unit current - electrode voltage
    Iinput_mat = zeros(nnode, nmuscle);
    musind = 1;
    for nodeID = 1 : nnode
        if electrode_index(nodeID) == 0
            Iinput_mat(nodeID,musind) = 1;
            musind = musind + 1;
        end
    end

    % fix ill conditioned LR matrix
    LR_tmp = LR;
    Lcond = cond(LR_tmp);
    fprintf("Conditioning number : %d\n", Lcond)
    while Lcond > MIN_CONDITION_NUMBER
        [Vv, Dd] = eig(LR_tmp);
        Dd = Dd + MIN_Z_EIGENVALUE.*eye(nnode);
        LR_tmp = Vv*Dd*Vv';
        Lcond = cond(LR_tmp);
    end

    % Calculate the voltages at ALL NODES
    Vmat = zeros(nnode, nmuscle);
    for muscleID = 1 : nmuscle
        Im = Iinput_mat(:,muscleID);
        Vtmp = inv(LR_tmp)*Im;
        Vmat(:,muscleID) = Vtmp;
    end
    
    % Definition of the lead field matrix (LFM). In our case it is the V matrix
    % 1st attemp: creation of the LFM ignoring the lines about the muscle
    % voltages
    LFM = Vmat;    % LFM is in Volts/Ampere  
    LFM_full = LFM;
    for nodeID = nnode : -1 : 1
        if electrode_index(nodeID) == 0
            LFM(nodeID,:)=[];
        end
    end