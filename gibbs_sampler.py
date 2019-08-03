#coding=utf-8

import datetime
import random
 
start = datetime.datetime.now()


def scoring_motif(matrix, m):
	ntdict={'A':matrix[0], 'C':matrix[1],'G':matrix[2],'T':matrix[3]}; score=1
#	print('len of matrix[0] is: {}'.format(len(matrix[0])))
#	print('len of motif is: {}'.format(len(m)))
	for i in range(len(m)): 
		score*=float(ntdict[m[i]][i])
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


def profile_random(seq, k, matrix):
	score= {}; l= []
	for i in range(len(seq)-k+1):
		score[seq[i:i+k]]= int(scoring_motif(matrix, seq[i:i+k])*100)
	for i in score.keys(): l.extend(score[i]*[i])
#	print(l)
	return random.choice(l)


def score(motifs):
#	print('scoring motifs')
#	print('the motifs is :', motifs)
	rows=[]; s=[]
	for num in range(len(motifs[0])): rows.append(''.join([i[num]for i in motifs if i != '']))
	for i in rows: 
		p= frequentletter(i); c=0
		for j in i: 
			if j !=p: c+=1
#		print('i:{}, s= {}'.format(i, c))
		s.append(c)
	return sum(s)


def gibbs_sampler(dnas, k, t, N):
	motifs=[]; c=0
	for i in dnas:
		r= random.randint(0, len(dnas[0])- k- 1); motifs.append(i[r:r+k])
	bestmotifs= motifs.copy()
	for j in range(N-1):
		i=random.randint(0, t-1)
		nmotifs= motifs.copy()
		nmotifs.pop(i)
		pssm= build_pssm_matrix(nmotifs)
		motifs[i]= profile_most(dnas[i], k, pssm)
		if score(motifs) < score(bestmotifs): 
			bestmotifs= motifs.copy()
	return bestmotifs


def loop(dnas, k, t, N):
	motifs= gibbs_sampler(dnas, k, t, N); newmotifs=[]
	for i in range(200):
#		print(i)
		newmotifs= gibbs_sampler(dnas, k, t, N)
		if score(newmotifs) < score(motifs):
			motifs= newmotifs.copy()
	print('score: ', score(motifs))
	return motifs

f=open('dataset_161_5.txt'); lines=f.readlines(); seq=[]; f.close()
for i in range(1, len(lines)): seq.append(lines[i].strip())

# seq=['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG','TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC','AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
print('\n'.join(loop(seq, 15, 20, 200)))
elapsed = (datetime.datetime.now() - start)
print("Time used:",elapsed)

