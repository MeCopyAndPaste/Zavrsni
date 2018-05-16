% Script calculates percetage of correct mds, suboptimal ds or incorrect ds
% found by algorithm run in function simulateDomset. 

% Graph is chosen from a set of n graphs in function chooseGraph. 
% Function simulateDomset runs n iterations of decisions on which node 
% becomes dominating set. Bees deciding in binary arenas simulated in 
% beeSimulation function within simulateDomset. Statistics calculated 
% on iter iterations. 

col = {[0 0.4470 0.7410];
[0.8500 0.3250 0.0980];
[0.9290 0.6940 0.1250]; 
[0.4940 0.1840 0.5560]; 
[0.4660 0.6740 0.1880];
[0 0.4470 0.7410]};

global casu_pos;
casu_pos = 4.5;

correct = 0;
incorrect = 0;
subopt = 0;
suboptStats = []; 

figure();
% iGraph = 6;
for iGraph = 7 : 7
    rho = 0.85;
    [N, mds] = chooseGraph(iGraph);
    n = 300;
    iter = 10;
    correct = 0;
    incorrect = 0;
    subopt = 0;
    
    subplot(4,2,iGraph)
    plot(graph(N));
    title(num2str(iGraph))
        
    scaling_heat = [];
    scalong_cool = [];
    stop_heat = 0.7;
    start_heat = 0.1;
    stop_cool = 0.5;
    start_cool = 0.2;
    for i = 1 : n
        progress_smooth_heat = 1 - exp(0.17 * (1 - 1 / (1 - min(1, 3 * i / n ))));
        progress_smooth_cool = 1 - exp(0.85 * (1 - 1 / (1 - min(1, 2 * i / n))));
        
        scaling_heat(i) = (1-progress_smooth_heat) * start_heat + stop_heat * progress_smooth_heat;
        scaling_cool(i) = (1-progress_smooth_cool) * start_cool + stop_cool * progress_smooth_cool;
    end
    
    suboptStats = [];
    
%     for i = 1 : iter
    i = 1;
    while i < iter
%         [vec, pStat] = simulateDomsetBees(N,n,rho);
        [vec, pStat, maxs, avgs] = simulateDomsetBees_crtl_vars(N,n,rho);
        
        [cor, sub, inc] = calculateStats(vec, N, mds);
%         if cor > 0
%             result = vec(end,:) > 35
%         end

        correct = correct + (cor > 0);
        subopt = subopt + (sub > 0);
        incorrect = incorrect + (inc > 0);
        if sub > 0
            suboptStats(subopt) = (sub - 1);
            i = iter;
            
            figure();
            n = length(vec);
            for c = 1 : length(N)
                subplot(length(N),1,c)
                scatter(linspace(1,n,n), vec(:,c), '.', 'markeredgecolor', col{mod(c,6)+1});
                axis([1, n, 26, 36])
            end

            figure();
            n = length(maxs);
            for c = 1 : length(N)
                subplot(length(N),1,c)
                scatter(linspace(1,n,n), maxs(:,c), '.', 'markeredgecolor', col{1});
                hold on
                scatter(linspace(1,n,n), scaling_cool, '.', 'markeredgecolor', col{2});
                scatter(linspace(1,n,n), avgs(:,c), '.', 'markeredgecolor', col{3});
                hold on
                scatter(linspace(1,n,n), scaling_heat, '.', 'markeredgecolor', col{4});
                axis([1, n, 0, 1])
            end
        else 
            i = i + 1;
        end
    end
    disp(strcat('graph: ', num2str(iGraph), ', rho = ', num2str(rho)));
    disp(strcat('cor:',num2str(correct*2),',so:',num2str(subopt*2),',inc:',num2str(incorrect*2)))
    disp(strcat('mean: ',num2str(mean(suboptStats)), ', var: ', num2str(var(suboptStats))))
    disp('--------------------------------')

end
% title(strcat('graph ',num2str(iGraph)));
% f = figure();
% plot(graph(N))

