'''
get closest neighbors ,the hamming distance between neighbor and pattern is just one.
function argument(s): pattern, return closest neighbors
'''
def Immediateneighbors(pattern): 
	neighborhood= [pattern]; NT=['A','T','C','G']
	for i in range(len(pattern)):
		symbol= pattern[i]
		for j in NT:
			if j!=symbol:
				npattern= list(pattern); npattern[i]= j
				seq=''.join(npattern); neighborhood.append(seq)
	return neighborhood


'''
get neighbors, the hamming distance between neighbor and pattern is no more than d
function argument(s): pattern d, return neighbors
'''
def iterativeNeighbors(pattern, d):
	neighborhood= [pattern]
	
	for i in range(int(d)):
		nneighborhood= neighborhood.copy() #赋值赋的是地址，所以这样赋值8行，记住！
#		print('{}  {}'.format(i, neighborhood))
		for j in nneighborhood:
#			print('j= ', j)
			neighborhood.extend(Immediateneighbors(j))
#			neighborhood.extend(str(Seq(str(Immediateneighbors(j))).reverse_complement()))
#			print('immediate once, the neighborhood now is :', neighborhood)
			neighborhood= list(set(neighborhood))
#	print('len of neighbor= ', len(neighborhood))
	return list(set(neighborhood))


# get all 4k_k_mers as pattern, argument(s): k, return 4k_k_mers
def build_4k_k_merset(k): return iterativeNeighbors('A'*int(k), int(k))


def distance(s1, s2):
	d=0
	for i in range(len(s1)):
		if s1[i]!=s2[i]: d+=1
	return d


def hammingd(j, i):
	d=[]
	for k in range(len(j)-len(i)+1):
		d.append(distance(j[k:k+len(i)], i))
	return min(d)


def multihammingd(l, i):
	d=[] 
	for j in l: d.append(hammingd(j, i))
	return max(d)


def MotifEnumeration(Dna, k, d):
	patterns= []; mers=build_4k_k_merset(k)
	for i in mers:
		if multihammingd(Dna, i)<= d: patterns.append(i)
	return ' '.join(list(set(patterns)))


# print(MotifEnumeration(['GCGGTAACTCTACTCGTCGACGATG','CGATGGCGACAAACGCCCCCGCATT','TCGAGCGACGAAGCGTAAATCTCTT','GCTCCATCTCCGAGGGATACCGAAT','GTGCTCTTTCTAATGAAGGTCGAAG','CGAGGGGTACCGGGCGTGTAACGAT'], 5, 1))
