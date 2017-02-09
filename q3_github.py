import os, random, sys


def build_debrujn(reads, kmer=30):
	"""
		Build a De Bruijn graph based on the sequenced reads.
		Default kmer is 30.
	"""
	kmer = 30
	nodes = set()
	edges = []
	for read in reads: # construct the graph from each read
		for i in range(len(read)-kmer+1):
			nodes.add(read[i:i+kmer-1])
			nodes.add(read[i+1:i+kmer])
			edges.append((read[i:i+kmer-1],read[i+1:i+kmer]))
	return nodes, edges

def pick_start_for_debruijn(nds,eds):
	"""
		Pick the node that is most likely to be the head of 
		the whole sequence.
		Return the node whose indegree is the fewest. 
		Randomly return one when there is a tie.
	"""
	nds = list(nds)
	res = [0] * len(nds) # record the tail number of eash node
	for n in range(len(nds)):
		for e in eds:
			if nds[n] == e[1]:
				res[n] += 1
	minnum = len(eds)
	tmp = []	
	for i in range(len(res)):
		if res[i] < minnum:
			minnum = res[i]
			tmp = [i]
		elif res[i] == minnum:
			tmp.append(i)
	return nds[random.choice(tmp)]

def dfs_rec(nds, eds, start, scs, visited):
	"""
		Recursively search the graph by depth.
		Return the set of the visited nodes and the superstring
	"""
	visited.add(start)
	if visited == nds: # all the nodes are visited
		return visited, scs

	noedge = True
	for e in eds:
		if e[0] == start and e[1] not in visited:
			noedge = False
	if noedge: # skip visited node although edge is not visited
		return visited, scs

	for edge in eds:
		if edge[0] == start and edge[1] not in visited:
			scs += edge[1][-1] # append the last character
			c = dfs_rec(nds, eds, edge[1], scs, visited)
			visited, scs = c[0],c[1]
	return visited, scs


def dfs(nds, eds):
	"""
		Initialize the parameters and invoke dfs_rec to search 
		the graph.
		Return the set of the visited nodes and the superstring
	"""	
	visited = set()
	start = pick_start_for_debruijn(nds-visited,eds)
	scs = start
	visited, scs = dfs_rec(nds, eds, start, scs, visited)
	while visited != nds: # until all the nodes are visited
		start = pick_start_for_debruijn(nds-visited,eds)
		visited, scs = dfs_rec(nds, eds, start, scs, visited)		
	return visited, scs

def check_overlap(scs, read):
	"""
		Check if the read overlaps with the superstring and 
		return the length of the overlapped part.
		If yes, return the shortest superstring of them.
		Else, return the concatenation of them.
	"""
	startl, startr = 0, 0
	fl, fr = True, True
	
	# the rightmost of read and the leftmost of scs
	while startl < len(read):
		if not scs.startswith(read[startl:]):
			startl += 1
		else:
			break
	
	# the rightmost of scs and the leftmost of read
	while startr < len(read):
		if not read.startswith(scs[-(len(read)-startr):]):
			startr += 1
		else:
			break
		
	if startl == len(read) and startr == len(read):
		return 0, random.choice([scs+read, read+scs])
	elif startl != len(read) and startr != len(read):
		if startl < startr:
			return len(read)-startl, read + scs[len(read)-startl:]
		else:
			return len(read)-startr, scs + read[len(read)-startr:]
	elif startr == len(read):
		return len(read)-startl, read + scs[len(read)-startl:]
	elif startl == len(read):
		return len(read)-startr, scs + read[len(read)-startr:]

def modify(scs, pick_read):
	"""
		Given scs and a list of all the unmapped reads.
		Return a new superstring containing all reads and scs.
	"""
	while pick_read != []:	
		maxovl = 0
		newscs = ''
		tmpread = ''
		for read in pick_read:
			tmp = check_overlap(scs, read)
			# pick one read which has the largest overlap with scs
			# one at a time
			if tmp[0] >= maxovl:
				maxovl = tmp[0]
				newscs = tmp[1]
				tmpread = read
		pick_read.pop(pick_read.index(tmpread))
		scs = newscs
	return scs

def unmapped_read(scs, reads):
	"""
		Given the superstring and all the reads.
		Return a list of the reads that cannot be mapped to scs
	"""
	pick_read = []
	for read in reads:
		if read in scs:
			continue
		pick_read.append(read)
	return pick_read

def check_existence(scs, reads):
	"""
		Given the superstring and all the reads.
		Return the shortest superstring of them
	"""
	pick_read = unmapped_read(scs, reads)
	scs = modify(scs, pick_read)
	return scs

	
def dfs_iter(nds, eds):
	"""
		Iteratively search the graph by depth.
		Return the set of the visited nodes and the superstring
	"""
	visited = set()
	# if a node's outdegree is bigger than 1, pick only one head
	# and store in the list stack
	# scs contains the contigs of the reads
	stack, scs = [], []
	while visited != nds:
		if stack == []:
			node = pick_start_for_debruijn(nds-visited,eds)	
			# treat node as the seed to grow to be a contig		
			scs.append(node)
		else:
			node = stack.pop()
			# append the last character of the node to the seed
			scs[-1] += node[-1]

		if node not in visited:
			visited.add(node)
			tmp = []
			# search for the next
			for edge in eds:
				if edge[0] == node and edge[1] not in visited:
					if edge[1] not in tmp:
						tmp.append(edge[1])
			# tmp stores all the heads from the node
			if len(tmp) >= 1: # pick one head randomly
				stack.append(random.choice(tmp))

	# check the contigs in the scs
	if len(scs) == 0:
		return visited, ''
	if len(scs) == 1:
		return visited, scs[0]

	while len(scs) != 1:
		maxlen = 0
		tmps1, tmps2, overlap = '','',''
		# find the overlap between each two contigs
		for s1 in scs:
			for s2 in scs:
				if s1 == s2:
					continue
				overlen = check_overlap(s1, s2)
				if overlen[0] >= maxlen:
					maxlen = overlen[0]
					overlap = overlen[1]
					tmps1, tmps2 = s1, s2
		# if no overlap between any two contigs
		if maxlen == 0 and len(scs) != 1:
			return visited, check_existence(scs[0], scs[1:])
		# combine the contigs with the largest overlap
		scs.append(overlap)
		scs.pop(scs.index(tmps1))
		scs.pop(scs.index(tmps2))
	return visited, scs[0]



def check_subsequence(s, t):
	"""
		s is the genome sequence and t is the assembled 
		quence.
		Return whether t is the subsequence of s
	"""
	idx = []
	i, j = 0, 0
	while i < len(s) and j < len(t):
		if s[i] == t[j]:
			i += 1
			j += 1
			idx.append(i)
		else:
			i += 1

	if j != len(t):
		return 0
	else:
		return 1

def recursive(genome, reads, nodes, edges):
	"""
		Find the superstring with recursive method.
	"""
	tmpscs = ''
	tmplen = len(reads)
	for _ in range(100):
		sup_scs =  dfs(nodes,edges)[1]
		true_scs = check_existence(sup_scs, reads)
		if check_subsequence(genome, true_scs):
			print 'True subsequence'
			return true_scs

def itereative(genome, reads, nodes, edges):
	"""
		Find the superstring with itereative method.
	"""
	tmpscs = ''
	tmplen = float('inf')
	for _ in range(100):
		sup_scs =  dfs_iter(nodes,edges)[1]
		true_scs = check_existence(sup_scs, reads)
		if check_subsequence(genome, true_scs):
			print 'True subsequence'
			return true_scs
	

def main():
	"""
		Input the name of the sequences file and the method.
		'recursive' of 'iterative'.
		Return the superstring.
	"""
	filename = sys.argv[1]
	methd = sys.argv[2]
	finput = open(os.path.join('tests',filename))
	genome = finput.readline()
	genome.rstrip()
	finput.readline()
	reads = finput.readlines()
	reads = [i[:-2] for i in reads]
	print 'done reading...'
	nodes, edges = build_debrujn(reads)
	print 'done building graph'
	print len(nodes), 'nodes'
	if methd == 'recursive':
		true_scs = recursive(genome, reads, nodes, edges)
	elif methd == 'iterative':
		true_scs = itereative(genome, reads, nodes, edges
	else:
		print "Please enter a valid method name: recursive or iterative!"
	print true_scs
	for read in reads:
		print true_scs.find(read, 0)+1

if __name__ == '__main__':
	main()







