# This code is written for doing the population dynamics simulation
# using the Julia programming language. We start wtih all sensitive cells and
# see the population dynamics (population size, diversity, etc.) under the presence
# and absence of drug.

## To Run the code use the following command
#  'julia Tamres_Population_Dynamics_Simulation.jl P_SR P_RS Heterogeneity InitialPopulation DrugInducedPlasticity METInducedPlasticity InitialFracResistant OutputFileName'
#   See the bash script for more details on the arguments to the script
#%%
# Import required packages
using Random
using Distributions
using Statistics
using StatsBase
#%%
# Define simulation constants
const num_replicates = 20
const max_time = 100
const initial_cells = 100
const carrying_capacity = 10^5
const drug_coeff = 1.0
const drug_const = 0.6

# Define the cells as a class with atributes of id, score, transition rates, proliferation and death rates
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

#Function to randomly sample score from a gaussian distribution
function sample_gaussian(mean,std)
    return rand(Normal(mean,std))
end

#Initializing function that returns a list of cells with scores randomly sampled from the specified distribution
function init_double_gaussian(mean1,std1,mean2,std2,t_rate_12,t_rate_21,initial_cells,prol_rate,death_rate,current_score_of_cell,init_fraction_of_res)
    List_of_cells = []
    if current_score_of_cell == "r"
        for i in 1:initial_cells
            push!(List_of_cells,cell(i,sample_gaussian(mean1,std1),sample_gaussian(mean2,std2),t_rate_12,t_rate_21,prol_rate,death_rate,rand(1:2)))
        end
    else
        for i in 1:initial_cells
            if rand() < init_fraction_of_res #initial fraction of resistant cells
                push!(List_of_cells,cell(i,sample_gaussian(mean1,std1),sample_gaussian(mean2,std2),t_rate_12,t_rate_21,prol_rate,death_rate,2))
            else
                push!(List_of_cells,cell(i,sample_gaussian(mean1,std1),sample_gaussian(mean2,std2),t_rate_12,t_rate_21,prol_rate,death_rate,current_score_of_cell))
            end

        end
    end

    return List_of_cells
end

#Function that returns the population diversity (# of cells for a particular score)
function diversity(List_of_cells)
    if length(List_of_cells) == 0
        #println("intial ",fit(Histogram,[],-5:0.5:5).weights)
        return fit(Histogram,[],-6:0.1:6).weights
    end

    a = []
    for i in List_of_cells
        if i.current_score == 1
            append!(a,i.score1)
        elseif i.current_score == 2
            append!(a,i.score2)
        end
    end
    return fit(Histogram,a,-6:0.1:6).weights
end

#Function for doing the simulation. Includes the death, proliferation and swtiching of the cells
function run_sim(List_of_cells,initial_cells,num_replicates,max_time,carrying_capacity,prol_rate,death_rate,mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,current_score_of_cell,drug_plast_prob_sus,drug_plast_prob_res)
    summary_stats = []
    n = 4 ## no. of random numbers to be generated
    pop_size = []
    num_death = []
    num_prol = []
    num_switch = []
    diversity_array=[]
    decision = ["p","d","s"] # p -> proliferation, d-> death, s-> swtiching
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
        if length(List_of_cells) == 0 # all cells have died
            break
        end
        shuffle!(List_of_cells) # randomly order the cells so that there is no inherent bias of selection of cells
        for (i,j) in enumerate(List_of_cells)
            #println(i.second.score)
            rand_array = rand(n)
            shuffle!(decision) # randomly order the proliferation, death and swtiching to reduce the bias of having a particular order
            for dec in decision
                if dec == "p"  ## proliferation
                    if rand_array[1]< j.p_rate *( 1 - length(List_of_cells)/carrying_capacity ) && (length(List_of_cells)+length(list_to_add) < carrying_capacity)
                        count += 1
                        if j.current_score == 1
                            sus_prol += 1
                        else
                            res_prol += 1
                        end
                        if current_score_of_cell == "r" ## mode of proliferation is random (daughter cells have S/R type randomly selected)
                            push!(list_to_add,cell(i,sample_gaussian(mean_dist1,std_dist1),sample_gaussian(mean_dist2,std_dist2),t_rate_12,t_rate_21,prol_rate,death_rate,rand(1:2)))
                        else ## mode of proliferation is parental (S/R of daughter cells depend on S/R state of parental cell)
                            push!(list_to_add,cell(i,sample_gaussian(mean_dist1,std_dist1),sample_gaussian(mean_dist2,std_dist2),t_rate_12,t_rate_21,prol_rate,death_rate,current_score_of_cell))
                        end
                    end
                elseif dec == "d" ## Death
                    if j.current_score == 1
                        if (ℯ^(drug_coeff*j.score1))/((ℯ^(drug_coeff*j.score1)) + drug_const)  < rand_array[2] # death due to drug
                            push!(list_to_remove,i)
                            sus_dead += 1
                            break
                        else ## Drug induced plasticity
                            if rand() < drug_plast_prob_sus
                                j.current_score = 2
                                sus_to_res += 1
                            end
                        end
                        if rand_array[3] < j.d_rate # Basal death rate
                            push!(list_to_remove,i)
                            sus_dead += 1
                            break
                        end
                    else
                        if (ℯ^(drug_coeff*j.score2))/((ℯ^(drug_coeff*j.score2)) + drug_const)  < rand_array[2]
                            push!(list_to_remove,i)
                            res_dead += 1
                            break
                        else ## MET induced plasticity
                            if rand() < drug_plast_prob_res
                                j.current_score = 1
                                res_to_sus += 1
                            end
                        end
                        if rand_array[3] < j.d_rate
                            push!(list_to_remove,i)
                            res_dead += 1
                            break
                        end
                    end
                else #switching
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
            append!(diversity_array,[zeros(24)])
        end
    end
    append!(summary_stats,(pop_size,num_prol,num_death,num_switch,diversity_array))
    return summary_stats
end

#%%

pop_initial = parse(Int64,ARGS[4]) # initial population type susceptible/resistant/randomly selected
if pop_initial == 3
	pop_initial = "r"
end
#%%
#println(pop_initial)
init_fraction_of_res = parse(Float64,ARGS[7]) #initial fraction of resistant cells  
drug_plast_prob_sus = parse(Float64,ARGS[5]) #drug induced plasticity
drug_plast_prob_res = parse(Float64,ARGS[6]) #MET induced plasticity

mean_dist1, std_dist1, mean_dist2, std_dist2, t_rate_12, t_rate_21, prol_rate, death_rate, current_score_of_cell = (-2, parse(Float64,ARGS[3]), 2, parse(Float64,ARGS[3]), parse(Float64,ARGS[1]), parse(Float64,ARGS[2]), 0.91, 0.1, pop_initial)
lock = Threads.SpinLock()
summary_stats_rep = []
Threads.@threads for rep in 1:num_replicates
    List_of_cells = init_double_gaussian(mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,initial_cells,prol_rate,death_rate,current_score_of_cell,init_fraction_of_res)
    summary_stats = run_sim(List_of_cells,initial_cells,num_replicates,max_time,carrying_capacity,prol_rate,death_rate,mean_dist1,std_dist1,mean_dist2,std_dist2,t_rate_12,t_rate_21,current_score_of_cell,drug_plast_prob_sus,drug_plast_prob_res)
    Threads.lock(lock)
    append!(summary_stats_rep, [summary_stats])
    Threads.unlock(lock)
end
#%%
#Saving the output to file

filename = ARGS[8]
f1 = open(filename*"_summary.txt","a")
f2 = open(filename*"_diversity.txt","a")
global count_rep = 0
for i in summary_stats_rep
    global count_rep;
    count_rep += 1
    write(f1,string("Replicate","\t",count_rep,"\n"))
    write(f2,string("Replicate","\t",count_rep,"\n"))
    write(f1,string("Population size","\t"))
    for j in i[1]
        write(f1,string(j,"\t")) # Write Population size time series
    end
    write(f1,"\n")
    write(f1,string("Prol_sus","\t"))
    for j in i[2]
        write(f1,string(j[1],"\t")) # Write Number of proliferation from susceptible population time series
    end
    write(f1,"\n")
    write(f1,string("Prol_res","\t"))
    for j in i[2]
        write(f1,string(j[2],"\t")) # Write Number of proliferation from resistant population time series
    end
    write(f1,"\n")
    write(f1,string("Death_sus","\t"))
    for j in i[3]
        write(f1,string(j[1],"\t")) # Write Number of death from susceptible population time series
    end
    write(f1,"\n")
    write(f1,string("Death_res","\t"))
    for j in i[3]
        write(f1,string(j[2],"\t")) # Write Number of death from resistant population time series
    end
    write(f1,"\n")
    write(f1,string("switch_sus_to_res","\t"))
    for j in i[4]
        write(f1,string(j[1],"\t")) # Write Number of switches from susceptible to resistant population time series
    end
    write(f1,"\n")
    write(f1,string("switch_res_to_sus","\t"))
    for j in i[4]
        write(f1,string(j[2],"\t")) # Write Number of switches from resistant to susceptible population time series
    end
    write(f1,"\n")
    for j in i[5]
        for k in j
            write(f2,string(k,"\t")) # Write diversity array time series to second file
        end
        write(f2,"\n")
    end
    write(f1,"\n")
    write(f1,"###\n")
    write(f2,"###\n")
end
close(f1)
close(f2)
