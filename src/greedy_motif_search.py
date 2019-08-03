
def scoring_motif(matrix, m):
#	print('{} {}'.format(len(matrix[0]), len(m)))
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
	rows=[]; s=[]
	for num in range(len(motifs[0])): rows.append(''.join([i[num]for i in motifs if i != '']))
	for i in rows: 
		p= frequentletter(i); c=0
		for j in i: 
			if j !=p: c+=1
#		print('i:{}, s= {}'.format(i, c))
		s.append(c)
	return sum(s)


def greedy_motif_search(dnas, k, t):
	bestmotifs=[i[0:k]for i in dnas]; motif=['']*t
	for num in range(len(dnas[0])-k+1):
#		print(num)
		nmotif=dnas[0][num:num+k]
		motif[0]=nmotif
		for i in range(1,t):
			matrix= build_pssm_matrix(motif[0:i])
#			for line in matrix: print(line) 
			motif[i]= profile_most(dnas[i], k, matrix)
#			print(motif[i])
#		for line in matrix: print(line)
		motifs=motif[:]
		sm=score(motifs); sb=score(bestmotifs)
#		print('motif:    {}, score:{}'.format(motifs, sm))
#		print('bestmotif:{}, score:{}'.format(bestmotifs, sb)); print('\n')
		if sm< sb:
			bestmotifs= motifs[:]
	return bestmotifs

f=open('dataset_160_9.txt'); lines=f.readlines(); seq=[]; f.close()
for i in range(1, len(lines)): seq.append(lines[i].strip())
print(' '.join(greedy_motif_search(seq, 12, 25)))
