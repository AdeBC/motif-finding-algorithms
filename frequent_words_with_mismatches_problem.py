from Bio.Seq import Seq
# calculate the neighbors of seq and reverse complement sequence seperately, remember!

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
def build_4k_k_merset(k):
	mer='A'*int(k)
	return(iterativeNeighbors('A'*int(k),int(k)))
	

'''
----------------------------------------------------------------------------------------------------------------------
four steps:
	1: get 4k_k_mers as pattern, initialization
	2: get neighbors of each pattern, and count the times their neithbors appear
	3: get c[sequence]+c[reverse_complement_sequence], c[reverse complement _sequence] is calculated in step 2
	4: search the max value, find and print result
---------------------------------------------------------------------------------------------------------------------	
'''
def Frewords_mismatchpro(string, k, d):
	patterns=build_4k_k_merset(k); c={}; neighbors={}; k=int(k); nc={}; print(len(patterns))
	for i in patterns:
		neighbors[i]=iterativeNeighbors(i, d)
		c[i]=0
		for j in range(len(string)-k+1):
			if string[j:j+k] in neighbors[i]: 
				c[i]+=1
	for i in c.keys():
			if str(Seq(str(i)).reverse_complement()) in c.keys():
				nc[i+' '+str(Seq(str(i)).reverse_complement())]=c[i]+c[str(Seq(str(i)).reverse_complement())]
	m=max(nc.values()); p=[] #改这里！分别算的结果相加最大，则取该二者。条件：互为反向互补链
	for i in nc.keys():
		if nc[i]==m:
			print('{}:{}'.format(i,nc[i])); p.append(i)
	r=set(' '.join(p).split()); print(' '.join(list(r)))


Frewords_mismatchpro('GGGGAATGGCAGGCAAATAATAATGGCACAAATGGGGAATAATAATAATGGCACAAATGGAATGGAATCAAATAATAATGGGGAATAATAATCACACACAGGGGAATGGGGGGGGAATAATGGGGCACAGGGGCAAATAATGGAATGGAATGGGGGGGGAATCAAATCAAATGGAATCAGGGGCACACAGGGGCAGGAATAATAATAATGGCAGGGGAATGGCAAATGGAATGGAATCAAAT',5, 2)
