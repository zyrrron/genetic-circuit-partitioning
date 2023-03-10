#!/usr/bin/env python

#	Copyright (C) 2021 by
#	Jing Zhang <jgzhang@bu.edu>, CIDAR Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Supporting modules
import argparse
import genetic_partition_test as gp 
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
import os
import datetime
import shutil
import copy

def main():

	begin_time = datetime.datetime.now()

	# Parse the command line inputs
	parser = argparse.ArgumentParser(description="perform graph partition using metis")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	parser.add_argument("-samples",  dest="samples",  required=True,  help="1,2", metavar="string")
	args = parser.parse_args()
	
	# Run the command
	samples = args.samples.split(',')
	settings = gp.load_settings(args.settings)
	for s in samples:
		print ('Processing sample', s)
		# obtain user-defined params 
		S_bounds        = settings[s]['S_bounds'].split(',')
		primitive_only  = settings[s]['primitive_only']
		maxNodes        = int(settings[s]['max_nodes_to_shift'])
		priority        = settings[s]['priority']
		loop_free       = 'FALSE'
		trajectories    = int(settings[s]['trajectories'])
		out_path        = settings[s]['output_path']
		timestep        = 20000

		print('S_bounds', S_bounds)
		print('maxNodes', maxNodes)

		# load graph 
		G   = gp.load_graph_undirected (settings, s)
		DAG = gp.load_graph (settings, s)
	
		outdir = out_path + '/nparts/'
		connectivity_order = gp.rank_connectivity (DAG, primitive_only, outdir)
		lowest_connectivity = min([c for (n, c) in connectivity_order])
		highest_connectivity = max([c for (n, c) in connectivity_order])
		print('median degree of connectivity in subgraphs', connectivity_order)
		print('lowest median degree of connectivity in subgraph', lowest_connectivity)
		print('highest median degree of connectivity in subgraph', highest_connectivity)

		# if primitive only is true, input and output nodes will be removed from list of vertices to be partitioned
		if primitive_only == 'TRUE':
			in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (DAG)
			G_primitive = gp.get_G_primitive (G, nonprimitives)
		elif primitive_only == 'FALSE':
			G_primitive   = copy.deepcopy(G)
			nonprimitives = []

		actual_n = {}
		# optimize graph partition results to satisfy constraints
		nparts = list(range(int(len(G_primitive.nodes())/int(S_bounds[1]))+1, int(len(G_primitive.nodes())/int(S_bounds[0]))+1))
		for npart in nparts:
			part_sol = outdir + str(npart) + '/part_solns.txt'
			if os.path.exists(part_sol):
				cut, partDict = gp.load_metis_part_sol (part_sol)
				actual_n[npart] = [len(partDict.keys())]
			# median_connectivity = [c for (n, c) in connectivity_order if int(n) == npart][0]
		sorted_n = sorted(actual_n.items(), key=lambda kv: kv[1])
		sorted_n = [tn for (tn, n) in sorted_n]
		# sorted_n = [n for (n, c) in connectivity_order][:5]
		for npart in range(10, 22): 
			opt_path = outdir + str(npart) + '/optimized/'
			part_sol = outdir + str(npart) + '/part_solns.txt'
			cut, partDict = gp.load_metis_part_sol (part_sol)
			min_conn, max_conn, median_conn = gp.get_connectivity (DAG, primitive_only, partDict)
			print('min conn', min_conn)
			print('max conn', max_conn)
			if maxNodes + int(S_bounds[1]) >= max([len(v) for v in list(partDict.values())]): 
				print('TARGET N', npart)
				# if os.path.exists(opt_path):
				# 	shutil.rmtree(opt_path)
				# 	os.mkdir(opt_path)
				# else: 
				# 	os.mkdir(opt_path)
				# gp.visualize_assignment_graphviz (G, [gp.get_part(partDict, n) for n in G.nodes()], nonprimitives, primitive_only, outdir+str(npart), 1, [])
				
				if os.path.exists(opt_path) == False: 
					if os.path.exists(opt_path) == False: os.mkdir(opt_path)
				
				for conn in range(int(median_conn), 11):
					print ('target connectivity', conn)
					opt_path = outdir + str(npart) + '/optimized/' + str(conn) + '/'
					# if os.path.exists(opt_path):
					# 	shutil.rmtree(opt_path)
					# 	os.mkdir(opt_path)

					# 	has_solution = False
					# 	result = gp.choose_best_iteration (DAG, primitive_only, partDict, [conn], loop_free, opt_path)
						
					# 	if result != []: 
					# 		besti, order, solDict = result[0], result[1], result[2]
					# 		tot_i = len(solDict.keys())
					# 		sol_i = 0

					# 		if len(order) > 3: order = order[0:3]

					# 		for iteration in order:
					# 			print('continuing itegration', iteration)
					# 			part = solDict[iteration]['part'] 
					# 			cut = solDict[iteration]['cut']
					# 			part_opt = [gp.get_part(part, n) for n in G_primitive.nodes()]
					# 			matrix, partG = gp.partition_matrix (G_primitive, part_opt)
					# 			cell_unmet_const, cell_met_const = gp.get_cells_unmet_constraint (matrix, partG, [conn], loop_free)
					# 			# print('cell_unmet', cell_unmet_const)
					# 			if len(cell_unmet_const) == 0:
					# 				print('has solution') 
					# 				has_solution = True
					# 			else:
					# 				print('optimizing iteration', iteration)
					# 				soln_found, recorded_path = gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, part, maxNodes, 3, 3, True, [conn], loop_free, priority, timestep, 1, tot_i+sol_i, opt_path)
					# 				if soln_found == False: 
					# 					if len(cell_unmet_const) < 3: 
					# 						soln_found, recorded_path = gp.split_cells (DAG, primitive_only, S_bounds, cut, partDict, iteration, maxNodes, [conn], loop_free, priority, trajectories, opt_path)
					# 				sol_i += recorded_path
					# 				if soln_found:
					# 					has_solution = True
					# 					print('optimizing has solution')

					# 			# if found a solution, break
					# 			if has_solution: 
					# 				print('has solution, breaking')
					# 				break 
					# 		else:
					# 			continue
					# 		break
					
					# else: 

					if os.path.exists(opt_path) == False:	
						os.mkdir(opt_path)

						print('running graph partitioning')
						gp.optimize_signal_subnetwork_tmp (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, 3, 3, True, [str(conn)], loop_free, priority, timestep, trajectories, opt_path)
						# solDict = gp.load_opt_part_sol (opt_path + 'part_solns.txt')
						# for iteration in solDict.keys():
						# 	part = solDict[iteration]['part'] 
						# 	part_opt = [gp.get_part(part, n) for n in G_primitive.nodes()]
						# 	matrix, partG = gp.partition_matrix (G_primitive, part_opt)
						# 	cell_unmet_const, cell_met_const = gp.get_cells_unmet_constraint (matrix, partG, [conn], loop_free)
						# 	if cell_unmet_const == []: 
						# 		break
						# else: 
						# 	continue
						# break 
				
		# gp.determine_best_solution_bionetwork (DAG, primitive_only, str(median_connectivity), outdir)


			# for npart in os.listdir(outdir):
			# 	if npart.isdigit():
			# 		print('target npart', npart)
			# 		part_sol = outdir + npart + '/part_solns.txt'
			# 		cut, partDict = gp.load_metis_part_sol (part_sol)
			# 		if len(list(G.nodes()))>=15:
			# 			print('optimizing by subnetworks')
			# 			print('optimizing with low constraint')
			# 			gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, low_constraint, loop_free, priority, trajectories, outdir + npart +'/optimized_lc/')
			# 			print('execution time', datetime.datetime.now()-begin_time)
			# 			# print('optimizing with high constraint')
			# 			# gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, high_constraint, loop_free, priority, trajectories, outdir + npart +'/optimized_hc/')
			# 			# print('execution time', datetime.datetime.now()-begin_time)
			# 		else:                                
			# 			print('brute-force optimizing')
			# 			print('optimizing with low constraint')
			# 			gp.optimize_signal_bruteforce (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, low_constraint, outdir + npart +'/optimized_lc/')
			# 			print('execution time', datetime.datetime.now()-begin_time)
			# 			print('optimizing with high constraint')
			# 			gp.optimize_signal_bruteforce (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, high_constraint, outdir + npart +'/optimized_hc/')
			# 			print('execution time', datetime.datetime.now()-begin_time)
			# gp.determine_best_solution (DAG, primitive_only, high_constraint, low_constraint, outdir)


	print('TOTAL execution time', datetime.datetime.now()-begin_time)

if __name__ == "__main__":
	main()


	