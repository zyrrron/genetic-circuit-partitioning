import csv
import random
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
import networkx as nx
import metis
from collections import Counter
import numpy as np
import time
import numpy.linalg as la
import scipy.cluster.vq as vq
import itertools
import operator
import math
import copy
import collections
from mpmath import *
from itertools import chain
from itertools import product
from itertools import starmap
from functools import partial
import os
import seaborn as sns
import shutil
from networkx.drawing.nx_agraph import graphviz_layout
import re
import ujson
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
import numpy.linalg as la
import scipy.cluster.vq as vq
import scipy


def read_json(inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	ports, gates = {}, {}
	for idx, line in enumerate(lines):
		line = line.strip()
		if line.startswith('"ports"'):
			p_s = idx
			searchlines = lines[idx+1:]
			for i, sl in enumerate(searchlines, idx):
				if sl.strip().startswith('"cells"'):
					p_e = i+1
		if line.startswith('"cells"'):
			g_s = idx
			searchlines = lines[idx+1:]
			for i, sl in enumerate(searchlines, idx):
				if sl.strip().startswith('"netnames"'):
					g_e = i
	# get information of inputs and outputs
	spacer = [idx+p_s+1 for idx, line in enumerate(lines[p_s+1:p_e]) if ': {' in line.strip()]
	for i, v in enumerate(spacer):
		# get names
		s = lines[v].strip()
		s = re.search('"(.*)"', s)
		el = s.group(1)
		ports[el] = {}
		# get directions
		s = lines[v+1].strip()
		s = re.search('"direction": "(.*)"', s)
		direction = s.group(1)
		ports[el]['direction'] = direction
		# get offset if it exists
		s = lines[v + 2].strip()
		if s[-1] == ',':
			offset = s.split(':')[1].split(',')[0].strip()
			ports[el]['offset'] = int(offset)
			s = lines[v + 3].strip()
		# get bits
		bits = s.split('[')[1].split(']')[0].strip()
		if ',' in bits:
			del ports[el]
			continue
		else:
			if bits == ["0"]:
				print("sssssssssss")
			else:
				ports[el]['bits'] = int(bits)
	# get information of gates
	spacer = [idx+g_s+1 for idx, line in enumerate(lines[g_s+1:g_e]) if '$abc$' in line.strip()]
	for i, v in enumerate(spacer):
		# get names
		s = int(lines[v].strip().split('"')[1].split('$')[-1])
		gates[s] = {}
		gates[s]['input'] = {}
		gates[s]['output'] = {}
		# search for attributes of this gate
		if i != len(spacer)-1:
			searchlines = lines[v:spacer[i+1]]
		else:
			searchlines = lines[v:]
		for sl in searchlines:
			# get gate type
			if sl.strip().startswith('"type"'):
				gatetype = re.search('_(.*)_', sl.strip())
				gates[s]['type'] = gatetype.group(1)
			# get input(s)
			if sl.strip().startswith('"A": [') or sl.strip().startswith('"B": ['):
				port = re.search('"(.*)"', sl).group(1)
				bits = sl.split('[')[1].split(']')[0].strip()
				gates[s]['input'][port] = int(bits)
			# get output
			if sl.strip().startswith('"Y": ['):
				port = re.search('"(.*)"', sl).group(1)
				bits = sl.split('[')[1].split(']')[0].strip()
				gates[s]['output'][port] = int(bits)
	return ports, gates


read_json("sha256.json")
