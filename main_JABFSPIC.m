tic
rng shuffle;
%% system parameters
% these parameters can be modified for different simulation setups
SNR = 3;
num_active = 3;
numAPs = 5;% Number of APs
Num_MonteCarlo = 100e4;
numClusters = 3; % Number of clusters of users   
numUser_cluster = 40;%the number of mUEs of each cluster
numSubCarriers = 20; %the number of subcarriers
num_slots = 7; 
load(['','angle.mat'],'angle')
%% other parameters
activity_rate = num_active/numUser_cluster;
numIterations = 100;
numIter = 20;
num_Iter = 2;
MaSL = 8; % maximum sparsity level
% channel_gain = zeros(numAPs,numUser_cluster,numSubCarriers,numClusters); % channel gain vector
% angle = zeros(numUser_cluster,numClusters);
% load([addr,'PL.mat'],'PL')
spreading_matrix = gen_ZadoffChu_Seq(numSubCarriers,numSubCarriers-1,numUser_cluster*numClusters);
% results in each Monte
error_detection = zeros(Num_MonteCarlo,numClusters);
Rel_error_Monte = zeros(Num_MonteCarlo,numClusters);
ser_Monte = zeros(Num_MonteCarlo,numClusters);
Rel_error_Monte_IC = zeros(Num_MonteCarlo,numClusters);
ser_Monte_IC = zeros(Num_MonteCarlo,numClusters);
num_false_alarm = zeros(Num_MonteCarlo,numClusters);
num_false_dismissal = zeros(Num_MonteCarlo,numClusters);
noise_power = 1e-4;
T_power = noise_power*10.^(SNR./10);
% parpool(12) % recommend using HPC for large number of pools
disp('Monte Carlo is started.')

for i_Monte = 1:Num_MonteCarlo
    if rem(i_Monte, 1000) == 0
        disp(i_Monte)
    end
    error_detection1 = zeros(1,numClusters);
    Rel_error_Monte1 = zeros(1,numClusters);
    ser_Monte1 = zeros(1,numClusters);
    Rel_error_Monte_IC1 = zeros(1,numClusters);
    ser_Monte_IC1 = zeros(1,numClusters);
    num_false_alarm1 = zeros(1,numClusters);
    num_false_dismissal1 = zeros(1,numClusters);
%     noise_power = transmit_power(1)/10.^(SNR(i_SNR)/10); % minimum receiver noise
    channel_gain_average = zeros(numAPs,numClusters);
    channel_gain_equ = zeros(numAPs,numUser_cluster,numSubCarriers,numClusters); % equivalent channel gain vector
    channel_gain_cat = zeros(numSubCarriers*numAPs,numUser_cluster,numClusters); % the catenating of equivalent channel gain on subcarriers
    num_err = 0;
    steering_mat = zeros(numAPs,numUser_cluster,numClusters);
    %%%%%%steering vectors
    for i_c = 1:numClusters
        steering_mat(:,:,i_c) = exp(-1j*pi*(0:numAPs-1).'...
            *sind(angle(:,i_c)).');
    end
    channel_gain_corr = zeros(numAPs,numAPs,numClusters);
    active_set = zeros(numUser_cluster,numClusters);
    %%%%%%generate the active set
    for i_cluster = 1:numClusters
        active_temp=sort(User_Active_Set_Generate(numUser_cluster, ...
            num_active, num_active, num_slots));
        active_set(active_temp(:,1) ,i_cluster)=1;
    end
    %%%%%%%%compute the equivalent channel matrix
    fading_random = 1/sqrt(2)*(randn(numUser_cluster,numSubCarriers,numClusters)+...
        1j*randn(numUser_cluster,numSubCarriers,numClusters));%fading
    % fading_random = randn(numUser_cluster,numSubCarriers,numClusters);%fading
    for i_cluster = 1:numClusters
        for i_subcarrier = 1:numSubCarriers
            channel_gain_equ(:,:,i_subcarrier,i_cluster)=steering_mat(:,:,i_cluster)*...
                diag(fading_random(:,i_subcarrier,i_cluster)...
                .*spreading_matrix(i_subcarrier,...
                (i_cluster-1)*numUser_cluster+1:i_cluster*numUser_cluster).');
            channel_gain_corr(:,:,i_cluster) = ...
                channel_gain_corr(:,:,i_cluster)+...
                channel_gain_equ(:,:,i_subcarrier,i_cluster)*...
                channel_gain_equ(:,:,i_subcarrier,i_cluster)';
            channel_gain_cat((i_subcarrier-1)*numAPs+1:i_subcarrier*numAPs,:,i_cluster)=...
                channel_gain_equ(:,:,i_subcarrier,i_cluster);
        end
%             channel_gain_average(:,i_cluster) = squeeze(mean(mean(channel_gain_equ(:,:,:,i_cluster),3),2));
        channel_gain_average(:,i_cluster)= mean(steering_mat(:,:,i_cluster),2);
    end
    transmit_power=T_power*ones(numUser_cluster,numClusters);
    %% generate the transmitted and received signal
    N_sig = 16; % 16QAM
    sig_N = zeros(numUser_cluster,num_slots,numClusters);
    sig = zeros(numUser_cluster,num_slots,numClusters);
    received_sig = zeros(numAPs,num_slots,numSubCarriers);
    received = zeros(numSubCarriers*numAPs,num_slots);
%         auto_corr = zeros(numAPs,numAPs,numClusters);
    received_sig_cluster = zeros(numAPs,num_slots,numClusters,numSubCarriers);
    for i_cluster=1:numClusters
        sig_N(:,:,i_cluster) = randi([0 N_sig-1],numUser_cluster,num_slots);
        sig(:,:,i_cluster) = diag(active_set(:,i_cluster).*...
            sqrt(transmit_power(:,i_cluster)/10))*...
            qammod(sig_N(:,:,i_cluster),N_sig);% signal
    end
    for i_subcarrier = 1:numSubCarriers
        for i_cluster=1:numClusters
           received_sig_cluster(:,:,i_cluster,i_subcarrier) = ...
                channel_gain_equ(:,:,i_subcarrier,i_cluster)*sig(:,:,i_cluster);
        end 
        received_sig(:,:,i_subcarrier) = squeeze(sum(received_sig_cluster(:,:,:,i_subcarrier),3))...
            +sqrt(noise_power/2)*(randn(numAPs,num_slots)...
            +1j*randn(numAPs,num_slots));
        received((i_subcarrier-1)*numAPs+1:i_subcarrier*numAPs,:)=received_sig(:,:,i_subcarrier);
    end
    %% statistical beamforming
    opt_factor=1;
    diag_loading_stat = numSubCarriers*5e-2/activity_rate;
    diag_loading_dyna = 1e-5;
    response = eye(numClusters);
    b_SBF = zeros(numAPs,numClusters);
    sum_channel = sum(channel_gain_corr,3);
    for i_cluster = 1:numClusters
        channel_gain_corr_sum = sum_channel-...
            opt_factor*channel_gain_corr(:,:,i_cluster)+diag_loading_stat*eye(numAPs);
        inv_channel_product=inv(channel_gain_corr_sum)*channel_gain_average(:,i_cluster);
        b_SBF(:,i_cluster) = inv_channel_product/...
            (channel_gain_average(:,i_cluster)'*inv_channel_product);
    end
     %% J-ABF-SP
    maxNum_beam_update = 20;
    alpha = 1;% averaging factor
    stopping_factor = 1e-5;
    signal_est = zeros(numUser_cluster,num_slots,numClusters);
    active_set_est = zeros(numUser_cluster,numClusters);
    error_cluster = zeros(numClusters,1);
    for i_cluster = 1:numClusters
        project_matrix_OABF = eye(numAPs)-channel_gain_average(:,i_cluster)*...
            channel_gain_average(:,i_cluster)'/norm(channel_gain_average(:,i_cluster))^2;
        Gamma = zeros(numUser_cluster,MaSL);
        c_s = zeros(numUser_cluster,num_slots,MaSL);
        error_s = zeros(MaSL,1);
        c_b = zeros(numUser_cluster,num_slots,maxNum_beam_update);
        Gamma_temp = zeros(numUser_cluster,1);
        for i_s = 1:MaSL
            Residual = zeros(num_slots*numSubCarriers,numIterations+1);
            Gamma_b = zeros(numUser_cluster,maxNum_beam_update);
            Gamma_amp = zeros(numUser_cluster,numIterations+1);
            Error_amp = zeros(numIterations+1,1);
            Gamma_amp(:,1) = Gamma_temp;
            b_n = b_SBF(:,i_cluster);
%             b_n = b_CBF(:,i_cluster);
%             b_n = b_ZFBF(:,i_cluster);
            b_temp=kron(eye(numSubCarriers),b_n);
            Y_n = (b_temp'*received).';
            B_n = b_temp'*channel_gain_cat(:,:,i_cluster);
            y_n = Y_n(:);
            D_n = kron(B_n,eye(num_slots));
            Residual(:,1) = y_n;
            Error_amp(1) = norm(Residual(:,1));
            c_n = zeros(numUser_cluster,num_slots,numIterations);
            error_b = zeros(maxNum_beam_update,1);
            count_beam_update=1;
            error_b(count_beam_update) =1e20;
            for iter = 1:numIterations
                c_n(:,:,iter) = zeros(numUser_cluster,num_slots);
                w = zeros(num_slots,numUser_cluster);
                energy_users_slots=abs(D_n'*Residual(:,iter)).^2;
                energy_user = zeros(numUser_cluster,1);
                for i_user = 1:numUser_cluster
                    energy_user(i_user) = sum(energy_users_slots((i_user-1)*num_slots+1:i_user*num_slots));
                end
                [~,set_sort_energy] = sort(energy_user,'descend');
                Lambda = union(find(Gamma_amp(:,iter)),set_sort_energy(1:i_s)); % Support estimate
                w(:,Lambda)=Inverse_vectorization...
                    (pinv(kron(B_n(:,Lambda),eye(num_slots)))*y_n,num_slots); %LS estimate
                energy = sum(abs(w).^2);
                [~,set_sort_LS]=sort(energy,'descend');
                Gamma_amp(set_sort_LS(1:i_s),iter+1) = 1; %Support pruning
                c_n(logical(Gamma_amp(:,iter+1)),:,iter)=Inverse_vectorization...
                    (pinv(kron(B_n(:,logical(Gamma_amp(:,iter+1))),eye(num_slots)))*y_n,...
                    num_slots).'; %Signal estimation
                c_n_temp=c_n(:,:,iter).';
                Residual(:,iter+1) = y_n-D_n*c_n_temp(:); %Residual update
                Error_amp(iter+1) = norm(Residual(:,iter+1));
                if Error_amp(iter+1)>Error_amp(iter)...
                        ||Error_amp(iter+1)==Error_amp(iter)
                    Error_amp(iter+1) = Error_amp(iter);
                    c_n(:,:,iter) = c_n(:,:,iter-1);
                    Residual(:,iter+1) = Residual(:,iter);
                    Gamma_amp(:,iter+1) = Gamma_amp(:,iter);
                    if abs(Error_amp(iter)-error_b(count_beam_update))/error_b(count_beam_update)<stopping_factor...
                            ||count_beam_update==maxNum_beam_update
                        break
                    else
                        count_beam_update = count_beam_update+1;
                        c_b(:,:,count_beam_update) = c_n(:,:,iter-1);
                        error_b(count_beam_update) = Error_amp(iter);
                        Gamma_b(:,count_beam_update) = Gamma_amp(:,iter);
                        R_n = zeros(numAPs);
                        for i_subcarrier = 1:numSubCarriers
                            int=received_sig(:,:,i_subcarrier)-...
                                channel_gain_equ(:,:,i_subcarrier,i_cluster)*...
                                c_b(:,:,count_beam_update);
                            R_n = R_n+int*int';
                        end
                        R_n = R_n+diag_loading_dyna*eye(numAPs);
                        inv_R_c=inv(R_n)*channel_gain_average(:,i_cluster);
                        b_n = alpha*inv_R_c/(channel_gain_average(:,i_cluster)'*inv_R_c);% beamforming update
                        Y_n = (kron(eye(numSubCarriers),b_n)'*received).';
                        B_n = kron(eye(numSubCarriers),b_n)'*channel_gain_cat(:,:,i_cluster);
                        y_n = Y_n(:);
                        D_n = kron(B_n,eye(num_slots)); % measurement update
                    end
                end
            end
            c_s(:,:,i_s) = c_b(:,:,count_beam_update);
            error_s(i_s) = error_b(count_beam_update);
            Gamma(:,i_s) = Gamma_b(:,count_beam_update);
            Gamma_temp = Gamma(:,i_s);
        end
        relative_error_s = [0;abs(diff(error_s))]./abs(error_s);
        relative_error_s1 = [0;abs(diff(error_s(end:-1:1)))]./abs(error_s(end:-1:1));
        [~,order_Rel_error] = sort(relative_error_s,'descend');
        num_candidate = 3;
        candidate_set = 1:MaSL;
        Range=zeros(MaSL,1);
        consecutive_average_range = zeros(MaSL,1);
        for i_s = 1:MaSL
            norm_row = sum(abs(c_s(logical(Gamma(:,i_s)),:,i_s)),2);
            Range(i_s) = max(norm_row)/min(norm_row);
            if Range(candidate_set(i_s))>3
                candidate_set(i_s)=0;
            end

        end
        candidate_set = candidate_set(candidate_set~=0);
        s_opt = candidate_set(error_s(candidate_set)==min(error_s(candidate_set)));
        active_set_est(:,i_cluster) = Gamma(:,s_opt(1));
        signal_est(:,:,i_cluster) = c_s(:,:,s_opt(1));
        error_cluster(i_cluster) = error_s(s_opt(1));
    end
    %% re-estimation IC
    s_est = zeros(numUser_cluster,num_slots,numClusters);
    for i_c = 1:numClusters
        s_est(:,:,i_c) = signal_est(:,:,i_c);
    end
    Err = zeros(numIter,numClusters);
    b_n = zeros(numAPs,numClusters);
    for i_c = 1:numClusters
        Err(1,i_c) = error_cluster(i_c);
        R_n = zeros(numAPs);
        for i_subcarrier = 1:numSubCarriers
            int = received_sig(:,:,i_subcarrier)-...
                channel_gain_equ(:,:,i_subcarrier,i_c)*...
                s_est(:,:,i_c);
            R_n = R_n+int*int';
        end
        R_n = R_n+diag_loading_dyna*eye(numAPs);
        inv_R_c=inv(R_n)*channel_gain_average(:,i_c);
        b_n(:,i_c) = inv_R_c/(channel_gain_average(:,i_c)'*inv_R_c);% beamforming initialisation
    end
    signal_temp = zeros(numUser_cluster,num_slots,numClusters,num_Iter+1);
    signal_temp(:,:,:,1) = s_est;
    for iter = 1:num_Iter
        for  i_c = 1:numClusters% interference cancellation
            diff_set=setdiff(1:numClusters,i_c);
            c_n_enhance = zeros(numUser_cluster,num_slots); 
            received_int_subcarrier = zeros(numAPs,num_slots,numSubCarriers);
            received_int = zeros(numSubCarriers*numAPs,num_slots);
            for i_sub=1:numSubCarriers
                for i_diff_clu=diff_set
                   received_int_subcarrier(:,:,i_sub) =...
                       received_int_subcarrier(:,:,i_sub)+...
                        channel_gain_equ(:,:,i_sub,i_diff_clu)*...
                        signal_temp(:,:,i_diff_clu,iter);
                end 
                received_int((i_sub-1)*numAPs+1:i_sub*numAPs,:)=...
                    received_int_subcarrier(:,:,i_sub);
            end 
            received_IC = received-received_int;
            for i_iter = 1:numIter
                B_n = kron(eye(numSubCarriers),b_n(:,i_c))'*channel_gain_cat(:,:,i_c);
                D_n = kron(B_n,eye(num_slots)); 
                Y_n = (kron(eye(numSubCarriers),b_n(:,i_c))'*received_IC).';
                Y_n1 = (kron(eye(numSubCarriers),b_n(:,i_c))'*received).';
                y_n = Y_n(:);% measurement update
                y_n1 = Y_n1(:);% measurement update
                c_n_enhance(logical(active_set_est(:,i_c)),:)=Inverse_vectorization...
                    (pinv(kron(B_n(:,logical(active_set_est(:,i_c))),eye(num_slots)))*y_n,...
                    num_slots).'; %Signal estimation
                c_n_temp=c_n_enhance.';
                Err(i_iter+1,i_c) = norm(y_n-D_n*c_n_temp(:));
                if Err(i_iter+1,i_c)<Err(i_iter,i_c)
                    s_est(:,:,i_c) = c_n_enhance;
                    R_n = zeros(numAPs);
                    for i_sub = 1:numSubCarriers
                        int=received_sig(:,:,i_sub)-...
                            channel_gain_equ(:,:,i_sub,i_c)*...
                            s_est(:,:,i_c);
                        R_n = R_n+int*int';
                    end
                    R_n = R_n+diag_loading_dyna*eye(numAPs);
                    inv_R_c=inv(R_n)*channel_gain_average(:,i_c);
                    b_n(:,i_c) = inv_R_c/(channel_gain_average(:,i_c)'*inv_R_c);% beamforming initialisation
                else
                    Err(1,i_c)=Err(i_iter,i_c);
                    break
                end
            end
            signal_temp(:,:,i_c,iter+1)=s_est(:,:,i_c);
        end 
    end
    for i_cluster = 1:numClusters
        set_active=find(active_set(:,i_cluster));
        set_detected = find(active_set_est(:,i_cluster));
        set_detected_active = intersect(set_detected,set_active);
        set_false_alarm = setdiff(set_detected,set_detected_active);
        set_false_dismissal = setdiff(set_active,set_detected_active);
        num_false_alarm1(i_cluster)=length(set_false_alarm);
        num_false_dismissal1(i_cluster)=length(set_false_dismissal);
        error_detection1(i_cluster) = (num_false_alarm1(i_cluster)+...
            num_false_dismissal1(i_cluster))/numUser_cluster;
        signal_equivanlence_IC=diag(sqrt(10./transmit_power(set_detected_active,i_cluster)))*...
            s_est(set_detected_active,:,i_cluster);
        sig_N_temp_IC = qamdemod(signal_equivanlence_IC,N_sig);
        signal_equivanlence=diag(sqrt(10./transmit_power(set_detected_active,i_cluster)))*...
            signal_est(set_detected_active,:,i_cluster);
        sig_N_temp = qamdemod(signal_equivanlence,N_sig);
        num_SE_active = length(find(sig_N_temp-...
            sig_N(set_detected_active,:,i_cluster)));
        num_SE_active_IC = length(find(sig_N_temp_IC-...
            sig_N(set_detected_active,:,i_cluster)));
        ser_Monte_IC1(i_cluster) = (num_SE_active_IC+...
            (num_false_dismissal1(i_cluster)+...
            num_false_alarm1(i_cluster))*num_slots)...
            /(num_slots*numUser_cluster);
        Rel_error_Monte_IC1(i_cluster) = ...
            sqrt(norm(s_est(set_detected,:,i_cluster)...
            -sig(set_detected,:,i_cluster))^2+...
            norm(sig(set_false_dismissal,:,i_cluster))^2)./...
            sqrt(norm(sig(set_active,:,i_cluster))^2+...
            norm(s_est(set_false_alarm,:,i_cluster))^2);
        ser_Monte1(i_cluster) = (num_SE_active+...
            (num_false_dismissal1(i_cluster)+...
            num_false_alarm1(i_cluster))*num_slots)...
            /(num_slots*numUser_cluster);
        Rel_error_Monte1(i_cluster) = ...
            sqrt(norm(signal_est(set_detected,:,i_cluster)...
            -sig(set_detected,:,i_cluster))^2+...
            norm(sig(set_false_dismissal,:,i_cluster))^2)./...
            sqrt(norm(sig(set_active,:,i_cluster))^2+...
            norm(signal_est(set_false_alarm,:,i_cluster))^2);
    end
    error_detection(i_Monte,:) = error_detection1;
    Rel_error_Monte(i_Monte,:) = Rel_error_Monte1;
    ser_Monte(i_Monte,:) = ser_Monte1;
    Rel_error_Monte_IC(i_Monte,:) = Rel_error_Monte_IC1;
    ser_Monte_IC(i_Monte,:) = ser_Monte_IC1;
    num_false_alarm(i_Monte,:) = num_false_alarm1;
    num_false_dismissal(i_Monte,:) = num_false_dismissal1;
end

disp('Monte Carlo is completed.')

disp('Processing results')

detect_error_rate= mean(error_detection);
ser = mean(squeeze(ser_Monte));
ser_IC = mean(squeeze(ser_Monte_IC));
Rel_error= mean(squeeze(Rel_error_Monte));
Rel_error_IC = mean(squeeze(Rel_error_Monte_IC));
total_num_false_alarm = sum(squeeze(num_false_alarm));
total_num_false_dismissal = sum(squeeze(num_false_dismissal));
%%
disp('Simulation completed.')
toc



 




