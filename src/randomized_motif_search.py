#coding=utf-8

import time
import random
 
start = time.clock()


def scoring_motif(matrix, m):
	ntdict={'A':matrix[0], 'C':matrix[1],'G':matrix[2],'T':matrix[3]}; score=1
	for i in range(len(m)): score*=float(ntdict[m[i]][i])
	return score


def frequentletter(row):
	c=[]
	for i in row: c.append(row.count(i))
	for i in c: 
		if i== max(c): 
			j=row[c.index(i)] 
#	print('row= {}, max(c)={}\nfrequentletter={}'.format(row, max(c), j))
	return j


def build_pssm_matrix(motifs):
	length=4; width=len(motifs[0]); rows=[]; matrix={'A':[],'C':[],'G':[],'T':[]}
	for num in range(width): rows.append(''.join([i[num]for i in motifs if i != '']))
	for i in range(width): 
		matrix['A'].append(1+float(rows[i].count('A'))/len(rows[0]))
		matrix['C'].append(1+float(rows[i].count('C'))/len(rows[0]))
		matrix['G'].append(1+float(rows[i].count('G'))/len(rows[0]))
		matrix['T'].append(1+float(rows[i].count('T'))/len(rows[0]))
	return [matrix[i] for i in ['A','C','G','T']]


def profile_most(seq, k, matrix):
	score=0; index=0
	for i in range(len(seq)-k+1):
		newscore=scoring_motif(matrix, seq[i:i+k])
		if newscore >score: 
#			print('{} {}'.format(newscore, score))
			score=newscore
			index=i
	return seq[index:index+k]


def score(motifs):
	print('scoring motifs')
	print('the motifs is :', motifs)
	rows=[]; s=[]
	for num in range(len(motifs[0])): rows.append(''.join([i[num]for i in motifs if i != '']))
	for i in rows: 
		p= frequentletter(i); c=0
		for j in i: 
			if j !=p: c+=1
#		print('i:{}, s= {}'.format(i, c))
		s.append(c)
	return sum(s)


def randomized_motif_search(dnas, k, t):

	motifs=[]; c=0
	for i in dnas:
		r= random.randint(0, len(dnas[0])- k- 1); motifs.append(i[r:r+k])
	bestmotifs= motifs.copy(); print('		the original motifs is: ', motifs)
	while 1:
		pssm= build_pssm_matrix(motifs)
		for i in range(t): motifs[i]= profile_most(dnas[i], k, pssm)
		print('		the profile motifs is: ', motifs)
		if score(motifs) < score(bestmotifs): 
			bestmotifs= motifs.copy()
#			print('		loop once, the motifs is :{}',format(motifs))
		else : 
			return bestmotifs


def loop(dnas, k, t):
	motifs= randomized_motif_search(dnas, k, t); newmotifs=[]
	for i in range(1000):
		newmotifs= randomized_motif_search(dnas, k, t)
		if score(newmotifs) < score(motifs):
			motifs= newmotifs.copy()
	return motifs

f=open('dataset_161_5.txt'); lines=f.readlines(); seq=[]; f.close()
for i in range(1, len(lines)): seq.append(lines[i].strip())

# seq=['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG','TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC','AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
print(' '.join(loop(seq, 15, 20)))
elapsed = (time.clock() - start)
print("Time used:",elapsed)
