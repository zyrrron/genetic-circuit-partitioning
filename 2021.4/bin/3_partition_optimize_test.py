#!/usr/bin/env python

#	Copyright (C) 2021 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.


# Supporting modules
import argparse
import genetic_partition_test as gp
import os
import csv
import time
import copy

# set the max run time for one sample as 1 minute = 60 seconds
max_runtime = 60

def main():

	# graphviz = GraphvizOutput()
	# graphviz.output_file = 'partition_RCA4.png'

	begin_time = time.time()

	# Parse the command line inputs
	parser = argparse.ArgumentParser(description="perform graph partition using metis")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	parser.add_argument("-samples",  dest="samples",  required=True,  help="1,2", metavar="string")
	args = parser.parse_args()
	
	# Run the command
	samples = args.samples.split(',')
	settings = gp.load_settings(args.settings)

	for s in samples:
		begin_time_sample = time.time()
		print ('Processing sample', s)
		# print (settings[s])
		# obtain user-defined params 
		S_bounds        = settings[s]['S_bounds'].split(',')
		target_n        = settings[s]['target_n'].split(',')
		primitive_only  = settings[s]['primitive_only']
		maxNodes        = int(settings[s]['max_nodes_to_shift'])
		high_constraint = settings[s]['high_constraint'].split(',')
		low_constraint  = settings[s]['low_constraint'].split(',')
		loop_free       = settings[s]['loop_free']
		priority        = settings[s]['priority']
		trajectories    = int(settings[s]['trajectories'])
		out_path        = settings[s]['output_path']
		timestep = 10000

		begin_time_current_step = time.time()
		# load graph 
		G   = gp.load_graph_undirected (settings, s)
		DAG = gp.load_graph (settings, s)

		in_nodes, out_nodes, nonprimitives = gp.get_nonprimitive_nodes (DAG)

		if primitive_only == 'TRUE':
			G_primitive = gp.get_G_primitive (DAG, nonprimitives)
		else:
			G_primitive = copy.deepcopy (DAG)


		outdir = out_path + '/nparts/'
		order = gp.rank_connectivity (DAG, primitive_only, outdir)
		# print('median degree of connectivity in subgraphs', order)

		order = gp.rank_constraint_met (DAG, primitive_only, high_constraint, loop_free, outdir)
		# print('percent of cells with unmet constraint under high constraint', order)

		order = gp.rank_constraint_met (DAG, primitive_only, low_constraint, loop_free, outdir)
		# print('percent of cells with unmet constraint under low constraint', order)

		# optimize graph partition results to satisfy constraints
		if target_n == ['']:
			nparts = [int(n) for n in os.listdir(outdir) if n.isdigit()]
			Smin, Smax = int(S_bounds[0]), int(S_bounds[1])
			if len(G_primitive.nodes()) <= Smax:
				print("*****************************************************")
				print("Stop partitioning! You can put all nodes in one cell!")
				print("*****************************************************")

			else:
				for npart in sorted(nparts):
					print('target npart', npart)
					if npart > 8:
						part_sol_path = outdir + str(npart) + '/part_solns.txt'
						cut, partDict = gp.load_metis_part_sol (part_sol_path)
						sol_file = outdir + str(npart) +'/optimized_lc/part_solns.txt'
						# print('optimizing by subnetworks')
						gp.optimize_signal_subnetwork_tmp (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, 3, 2, True, low_constraint, loop_free, priority, timestep, trajectories, outdir + str(npart) +'/optimized_lc/')
						# print('execution time', time.time()-begin_time_current_step)
					else:
						# print('target npart', npart)
						part_sol_path = outdir + str(npart) + '/part_solns.txt'
						cut, partDict = gp.load_metis_part_sol (part_sol_path)
						sol_file = outdir + str(npart) +'/optimized_lc/part_solns.txt'
						# print('optimizing by subnetworks')
						gp.optimize_signal_subnetwork_tmp (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, 3, 2, True, low_constraint, loop_free, priority, timestep, trajectories, outdir + str(npart) +'/optimized_lc/')
						# print('execution time', time.time()-begin_time_current_step)
						gp.optimize_signal_subnetwork_tmp (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, 3, 2, True, high_constraint, loop_free, priority, timestep, trajectories, outdir + str(npart) +'/optimized_hc/')
						# print('execution time', time.time()-begin_time_current_step)


				# check solutions, if no solution, consider splitting a cell
				solDict = gp.load_opt_part_sol (sol_file)
				for iteration in solDict:
					part = solDict[iteration]['part']
					cut = solDict[iteration]['cut']
					part_opt = [gp.get_part(part, n) for n in G_primitive.nodes()]
					matrix, partG = gp.partition_matrix (G_primitive, part_opt)
					cell_unmet_const, cell_met_const = gp.get_cells_unmet_constraint (matrix, partG, low_constraint, loop_free)
					loop_free_i, motif_allowed = gp.check_constraint (matrix, partG, low_constraint)
					print('cell_unmet', cell_unmet_const)
					end_time_sampel = time.time()

					# set timer for the iteration
					if end_time_sampel - begin_time_sample > max_runtime:
						print("Time runs out!")
						break

					if 2 >= len(cell_unmet_const) > 0:
						gp.split_cells (DAG, primitive_only, S_bounds, cut, partDict, iteration, maxNodes, low_constraint, loop_free, priority, trajectories, outdir + str(npart) +'/optimized_lc/')

					if len(cell_unmet_const) == 0 and loop_free_i and motif_allowed:
						print('has solution')
						break

		else:
			for n in target_n:
				print('target npart', n)
				outdir = out_path + '/nparts/'
				part_sol = outdir + n + '/part_solns.txt'
				cut, partDict = gp.load_metis_part_sol (part_sol)
				if len(list(G.nodes()))>=15:
					# print('optimizing by subnetworks')
					# print('optimizing with low constraint')
					gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, 3, True, low_constraint, loop_free, priority, timestep, trajectories, outdir + n + '/optimized_lc/')
				else:
					print('brute-force optimizing')
					gp.optimize_signal_bruteforce (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, low_constraint, outdir + n + '/optimized_lc/')
					gp.optimize_signal_bruteforce (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, low_constraint, outdir + n + '/optimized_hc/')

		gp.determine_best_solution (DAG, primitive_only, high_constraint, low_constraint, outdir)
		runtime = time.time() - begin_time_current_step
		PATH = '/home/cidar-lab/genetic-circuit-partitioning/2021.4'
		outcsvpath = f"{PATH}/runs/results/oriole_time.csv"
		with open(outcsvpath, 'a', newline='') as f:
			fieldnames = ['ID', 'runtime']
			writer = csv.DictWriter(f, fieldnames=fieldnames)
			writer.writerow({'ID': s, 'runtime': runtime})

	runtime_total = time.time() - begin_time
	print('TOTAL execution time', runtime_total)




if __name__ == "__main__":
	main()


	