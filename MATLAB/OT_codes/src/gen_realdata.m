function [realdata] = gen_realdata(simparams)
    % About: Generates the real state values and the observations
    % delt: interval between two observations
    % Nt: total number of observations
    % Assumption:
    % time starts from T=0.
    % first observation is available at T=1
    %% Final truth data after propagation 
    x_start = [simparams.R0;simparams.rd; simparams.th0;simparams.thd];
    [~,xburn] = ode45(simparams.dyn,[0 simparams.burnperiod],x_start,[]);
    x_start = xburn(end,:)';
    x_true = zeros(simparams.Nt,4);
    for i = 1:simparams.Nt
        [t,x1] = ode45(simparams.dyn,[0 i*simparams.delt],x_start,[]);
        x_true(i,:) = x1(end,:);
        if i == simparams.Nt
           fullorbit.loc = x1;
           fullorbit.time = t;
        end
    end
    realdata.x_true = x_true;
    realdata.x_start = x_start;
    % Observations
    z_nonoise = (simparams.H*x_true')';
    meas_noise = mvnrnd([0;0],simparams.meas_sig,simparams.Nt);
    z_noise = z_nonoise + meas_noise;
    realdata.z_noise = z_noise;
    realdata.fullorbit = fullorbit;
    realdata.delt = simparams.delt;
    realdata.Nt = simparams.Nt ;
end