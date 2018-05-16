function vec = simulateDomsetBees(N,parameters)

    global casu_pos

    n = parameters(1);
    pIntMin = parameters(2);
    pIntMax = parameters(3);
    max_ctrl = parameters(4);
    start_heat = parameters(5);
    stop_heat = parameters(6);
    heat_crossover = parameters(7);
    start_cool = parameters(8);
    stop_cool = parameters(9);
    cool_crossover = parameters(10);
    rho = parameters(11);
    dT_heat = parameters(12);
    dT_cool = parameters(13);
    
    T = 28 * ones(size(N));
    P = zeros(size(N));
    pInt = P;

    nBees = 15;
    posD = randn(sum(sum(N))/2,nBees);
    posA = randn(sum(sum(N))/2,nBees);
    vel = randn(sum(sum(N))/2,nBees);

    ctrl = zeros(length(N),1);
    flag = ctrl;

    heat = ctrl';
    heat_float = ctrl';
    cool = ctrl';
    cool_float = ctrl';
    
    maxP = heat;
    avgP = heat;
%     vec = zeros(n, length(N));
    
    %% (n * k seconds)
    for i = 1 : n
        iPair = 0;
        for iSmaller = 1 : length(N)
            for iLarger = iSmaller + 1 : length(N)
                if N(iSmaller,iLarger) == 1
                    iPair = iPair + 1;
                    Tsubarena = [T(iSmaller),T(iLarger)];
                    k = 10;
                    for iTime = 1 : k
                        %% move bees k times
                        [posD(iPair,:),posA(iPair,:),vel(iPair,:)] = beeSimulation...
                            (posD(iPair,:),...
                            posA(iPair,:),...
                            vel(iPair,:),...
                            Tsubarena, 0);                    
                    end
                    % percentage of bees around a casu
                    P(iSmaller,iLarger) = ...
                        sum((posA(iPair,:).*cos(posD(iPair,:)) + casu_pos).^2 +...
                        (posA(iPair,:).*sin(posD(iPair,:))).^2 < 4) / nBees;
                    P(iLarger,iSmaller) = ...
                        sum((posA(iPair,:).*cos(posD(iPair,:)) - casu_pos).^2 +...
                        (posA(iPair,:).*sin(posD(iPair,:))).^2 < 4) / nBees;
                end
            end
        end
        %% calculate action depending on bees around a casu
        
        pInt = pInt + P;
        ctrl = (ctrl + (sum(pInt,2) > pIntMin) & (sum(pInt,2) < pIntMax)) .* (1-flag);
        flag = ((flag == 0) & (ctrl == max_ctrl)) | flag;

        arenaCnt = sum(N,2);
        for iNode = 1 : length(N)
            maxP(iNode) = max(P(iNode, N(iNode, :)==1));
            avgP(iNode) = mean(P(iNode, N(iNode, :)==1));
        end

        progress_heat = 1 - exp(heat_crossover * (1 - 1 / (1 - i / n)));
        progress_cool = 1 - exp(cool_crossover * (1 - 1 / (1 - i / n)));
        scaling_heat = (1-progress_heat) * start_heat + stop_heat * progress_heat;
        scaling_cool = (1-progress_cool) * start_cool + stop_cool * progress_cool;

        heat_float = heat_float * (1-rho) + (avgP > scaling_heat | ctrl') * rho;
        heat = heat_float > 0.5;

        cool_float = cool_float * (1-rho) + (maxP < scaling_cool .* (heat == 0)) * rho;
        cool = cool_float > 0.5;

        T(:,1) = min(36, T(:,1) + dT_heat * arenaCnt .* (heat') .* (sum(pInt,2) > 1));% - arenaCnt == zeros(size(heat))));
        T(:,1) = max(26, T(:,1) - dT_cool * (cool') .* (sum(pInt,2) > 1));% - arenaCnt == zeros(size(cool))));

        vec(i,:) = T(:,1);

        for c = 1 : length(N)
            T(c,:) = T(c,1);
        end
    end

end