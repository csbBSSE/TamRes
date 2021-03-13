using Random
using Distributions
using Statistics
using StatsBase
using Plots
pyplot()
#%%
const num_replicates = 10
const max_time = 100
const initial_cells = 1000
const carrying_capacity = 10^5
const drug_coeff = 1.0
const drug_const = 0.6

mutable struct cell
    id::Int64
    score1::Float64
    score2::Float64
    t_rate_1_to_2::Float64
    t_rate_2_to_1::Float64
    p_rate::Float64
    d_rate::Float64
    current_score::Int64 ## 1 -> score1, 2 -> score2
end

function sample_gaussian(mean,std)
    return rand(Normal(mean,std))
end

function init_double_gaussian(mean1,std1,mean2,std2,t_rate_12,t_rate_21,initial_cells,prol_rate,death_rate,current_score_of_cell)
    List_of_cells = []
    if current_score_of_cell == "r"
        for i in 1:initial_cells
            push!(List_of_cells,cell(i,sample_gaussian(mean1,std1),sample_gaussian(mean2,std2),t_rate_12,t_rate_21,prol_rate,death_rate,rand(1:2)))
        end
    else
        for i in 1:initial_cells
            push!(List_of_cells,cell(i,sample_gaussian(mean1,std1),sample_gaussian(mean2,std2),t_rate_12,t_rate_21,prol_rate,death_rate,current_score_of_cell))
        end
    end

    return List_of_cells
end

function diversity(List_of_cells)
    if length(List_of_cells) == 0
        return fit(Histogram,[],-6:0.5:6).weights
    end

    a = []
    for i in List_of_cells
        if i.current_score == 1
            append!(a,i.score1)
        elseif i.current_score == 2
            append!(a,i.score2)
        end
    end
    return fit(Histogram,a,-6:0.5:6).weights
end

function run_sim(List_of_cells,initial_cells,num_replicates,max_time,carrying_capacity,prol_rate,death_rate,mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,current_score_of_cell,drug_status)
    summary_stats = []
    n = 4 ## no. of random numbers to be generated
    pop_size = []
    num_death = []
    num_prol = []
    num_switch = []
    diversity_array=[]
    decision = ["p","d","s"]
    for iter in 1:max_time
        list_to_add = []
        list_to_remove =[]
        count = initial_cells
        sus_dead = 0
        res_dead = 0
        sus_prol = 0
        res_prol = 0
        sus_to_res = 0
        res_to_sus = 0
        if length(List_of_cells) == 0
            break
        end
        shuffle!(List_of_cells)
        for (i,j) in enumerate(List_of_cells)
            rand_array = rand(n)
            shuffle!(decision)
            for dec in decision
                if dec == "p"
                    if rand_array[1]< j.p_rate *( 1 - length(List_of_cells)/carrying_capacity ) && (length(List_of_cells)+length(list_to_add) < carrying_capacity)
                        count += 1
                        if j.current_score == 1
                            sus_prol += 1
                        else
                            res_prol += 1
                        end
                        if current_score_of_cell == "r"
                            push!(list_to_add,cell(i,sample_gaussian(mean_dist1,std_dist1),sample_gaussian(mean_dist2,std_dist2),t_rate_12,t_rate_21,prol_rate,death_rate,rand(1:2)))
                        else
                            push!(list_to_add,cell(i,sample_gaussian(mean_dist1,std_dist1),sample_gaussian(mean_dist2,std_dist2),t_rate_12,t_rate_21,prol_rate,death_rate,current_score_of_cell))
                        end
                    end
                elseif dec == "d"
                    if j.current_score == 1
                        if drug_status != 0
                            if (ℯ^(drug_coeff*j.score1))/((ℯ^(drug_coeff*j.score1)) + drug_const)  < rand_array[2]
                                push!(list_to_remove,i)
                                sus_dead += 1
                                break
                            end
                        end
                        if rand_array[3] < j.d_rate
                            push!(list_to_remove,i)
                            sus_dead += 1
                            break
                        end
                    else
                        if drug_status != 0
                            if (ℯ^(drug_coeff*j.score2))/((ℯ^(drug_coeff*j.score2)) + drug_const)  < rand_array[2]
                                push!(list_to_remove,i)
                                res_dead += 1
                                break
                            end
                        end
                        if rand_array[3] < j.d_rate
                            push!(list_to_remove,i)
                            res_dead += 1
                            break
                        end
                    end
                else
                    if j.current_score == 1
                        if rand_array[4] < (j.t_rate_1_to_2)
                            j.current_score = 2
                            sus_to_res += 1
                        end
                    else
                        if rand_array[4] < (j.t_rate_2_to_1)
                            j.current_score = 1
                            res_to_sus += 1
                        end
                    end
                end
            end
        end

        append!(num_switch,[[sus_to_res,res_to_sus]])
        append!(num_death,[[sus_dead,res_dead]])
        append!(num_prol,[[sus_prol,res_prol]])
        deleteat!(List_of_cells,list_to_remove)
        List_of_cells = vcat(List_of_cells,list_to_add)
        append!(pop_size,length(List_of_cells))
        append!(diversity_array,[diversity(List_of_cells)])
    end
    if length(pop_size) != max_time
        for i in 1:(max_time - length(pop_size))
            append!(pop_size, 0)
            append!(num_switch, [[0,0]])
            append!(num_death, [[0,0]])
            append!(num_prol, [[0,0]])
            append!(diversity_array,[zeros(25)])
        end
    end
    append!(summary_stats,(pop_size,num_prol,num_death,num_switch,diversity_array))
    return summary_stats
end

pop_initial = 1 # 1 for sus,  2 for res, 3 for random
if pop_initial == 3
	pop_initial = "r"
end
#%%
drug_status = 1
mean_dist1, std_dist1, mean_dist2, std_dist2, t_rate_12, t_rate_21, prol_rate, death_rate, current_score_of_cell = (-2, 0, 2, 0, 0, 0, 0.9, 0.1, pop_initial)
lock = Threads.SpinLock()
summary_stats_rep = []
Threads.@threads for rep in 1:num_replicates
    List_of_cells = init_double_gaussian(mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,initial_cells,prol_rate,death_rate,current_score_of_cell)
    summary_stats = run_sim(List_of_cells,initial_cells,num_replicates,max_time,carrying_capacity,prol_rate,death_rate,mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,current_score_of_cell,drug_status)
    Threads.lock(lock)
    append!(summary_stats_rep, [summary_stats])
    Threads.unlock(lock)
end

plast_00_summary = summary_stats_rep
plast_00_pop = []
for i in plast_00_summary
    push!(plast_00_pop,i[1])
end
plast_00_pop_mean = mean(plast_00_pop)
plast_00_pop_std = 0.434*std(plast_00_pop)./plast_00_pop_mean
plast_00_pop_mean

#%%
drug_status = 1
mean_dist1, std_dist1, mean_dist2, std_dist2, t_rate_12, t_rate_21, prol_rate, death_rate, current_score_of_cell = (-2, 0, 2, 0, 0, 1, 0.9, 0.1, pop_initial)
lock = Threads.SpinLock()
summary_stats_rep = []
Threads.@threads for rep in 1:num_replicates
    List_of_cells = init_double_gaussian(mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,initial_cells,prol_rate,death_rate,current_score_of_cell)
    summary_stats = run_sim(List_of_cells,initial_cells,num_replicates,max_time,carrying_capacity,prol_rate,death_rate,mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,current_score_of_cell,drug_status)
    Threads.lock(lock)
    append!(summary_stats_rep, [summary_stats])
    Threads.unlock(lock)
end

plast_01_summary = summary_stats_rep
plast_01_pop = []
for i in plast_01_summary
    push!(plast_01_pop,i[1])
end
plast_01_pop_mean = mean(plast_01_pop)
plast_01_pop_std = 0.434*std(plast_01_pop)./plast_01_pop_mean
plast_01_pop_mean

#%%
drug_status = 1
mean_dist1, std_dist1, mean_dist2, std_dist2, t_rate_12, t_rate_21, prol_rate, death_rate, current_score_of_cell = (-2, 0, 2, 0, 1, 0, 0.9, 0.1, pop_initial)
lock = Threads.SpinLock()
summary_stats_rep = []
Threads.@threads for rep in 1:num_replicates
    List_of_cells = init_double_gaussian(mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,initial_cells,prol_rate,death_rate,current_score_of_cell)
    summary_stats = run_sim(List_of_cells,initial_cells,num_replicates,max_time,carrying_capacity,prol_rate,death_rate,mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,current_score_of_cell,drug_status)
    Threads.lock(lock)
    append!(summary_stats_rep, [summary_stats])
    Threads.unlock(lock)
end

plast_10_summary = summary_stats_rep
plast_10_pop = []
for i in plast_10_summary
    push!(plast_10_pop,i[1])
end
plast_10_pop_mean = mean(plast_10_pop)
plast_10_pop_std = 0.434*std(plast_10_pop)./plast_10_pop_mean
plast_10_pop_mean

#%%
drug_status = 1
mean_dist1, std_dist1, mean_dist2, std_dist2, t_rate_12, t_rate_21, prol_rate, death_rate, current_score_of_cell = (-2, 0, 2, 0, 1, 1, 0.9, 0.1, pop_initial)
lock = Threads.SpinLock()
summary_stats_rep = []
Threads.@threads for rep in 1:num_replicates
    List_of_cells = init_double_gaussian(mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,initial_cells,prol_rate,death_rate,current_score_of_cell)
    summary_stats = run_sim(List_of_cells,initial_cells,num_replicates,max_time,carrying_capacity,prol_rate,death_rate,mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,current_score_of_cell,drug_status)
    Threads.lock(lock)
    append!(summary_stats_rep, [summary_stats])
    Threads.unlock(lock)
end

plast_11_summary = summary_stats_rep
plast_11_pop = []
for i in plast_11_summary
    push!(plast_11_pop,i[1])
end
plast_11_pop_mean = mean(plast_11_pop)
plast_11_pop_std = 0.434*std(plast_11_pop)./plast_11_pop_mean

#%%
plot(1:max_time,plast_00_pop_mean.+1,ribbon=plast_00_pop_std,fmt=:png, yaxis=:log,dpi = 800,legend = false, grid = false)
plot!(1:max_time,plast_01_pop_mean.+1,ribbon=plast_01_pop_std,fmt=:png, yaxis=:log,dpi = 800,legend = false, grid = false)
plot!(1:max_time,plast_10_pop_mean.+1,ribbon=plast_10_pop_std,fmt=:png, yaxis=:log,dpi = 800,legend = false, grid = false)
plot!(1:max_time,plast_11_pop_mean.+1,ribbon=plast_11_pop_std,fmt=:png, yaxis=:log,dpi = 800,legend = false, grid = false)
ylims!((0.8,10^5))
xlims!((0,40))
savefig("Figure5C.png")
