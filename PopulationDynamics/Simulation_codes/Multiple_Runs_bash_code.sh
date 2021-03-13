export JULIA_NUM_THREADS=10
p1_min=0 #S->R switching prob
p1_max=1
p2_min=0 #R->S switching prob
p2_max=1
p1_num_points=10 # Total 11 points including 0
p2_num_points=10 # Total 11 points including 0
p1=$p1_min
p2=$p2_min
dp1=0.1
dp2=0.1
dhet=0.25
het_num_points=4 #Heterogeneity
het=0
initial_pop=1
init_frac=0 #Initial fraction of resistant cells
for ini1 in $( eval echo {0.1,0.50,1.0}); do
    for ini2 in $( eval echo {0.1,0.50,1.0}); do
	    p1=$p1_min
	    p2=$p2_min
	    het=0
	    initial_pop=1
	    for d in $( eval echo {0..$het_num_points}); do
		    p1=$p1_min
		    p2=$p2_min
		    for i in $( eval echo {0..$p1_num_points}); do
			    p2=$p2_min
			    for j in $( eval echo {0..$p2_num_points}); do
				    julia tamres_pop_dynamics_final_code_v6_drug_met_both_switching.jl $p1 $p2 $het $initial_pop ${ini1} ${ini2} ${init_frac} out_tamres_met_both_${ini1}_${ini2}_${het}_${p1}_${p2}
				    echo $(date) "Drug Induced plasticity prob =" $ini1 " MET Induced plasticity prob =" $ini2 "Hetrogenity = "$het "p1 = "$p1 "p2 = "$p2
				    p2=`echo $p2 + $dp2 | bc`
			    done
			    p1=`echo $p1 + $dp1 | bc`
		    done
		    het=`echo $het + $dhet | bc`
	    done
    done	
done
