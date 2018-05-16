% Script calculates percetage of correct mds, suboptimal ds or incorrect ds
% found by algorithm run in function simulateDomset. 

% Graph is chosen from a set of n graphs in function chooseGraph. 
% Function simulateDomset runs n iterations of decisions on which node 
% becomes dominating set. Bees deciding in binary arenas simulated in 
% beeSimulation function within simulateDomset. Statistics calculated 
% on iter iterations. 

function [correct, subopt, incorrect] = MdsStats(parameters, iter)

    global casu_pos;
    casu_pos = 4.5;

    correct = 0;
    incorrect = 0;
    subopt = 0;

    for iGraph = 1 : 7
        
        [N, mds] = chooseGraph(iGraph);

        for i = 1 : iter
        
            vec = simulateDomsetBees(N, parameters);

            [cor, sub, inc] = calculateStats(vec, N, mds);
            correct = correct + (cor > 0);
            subopt = subopt + (sub > 0);
            incorrect = incorrect + (inc > 0);
        end
    end
end