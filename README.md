# Motif-finding

Motif finding problem is a classical bioinformatics problem, aiming to quickly find a series of motifs on genes with the same enzyme (DNA replicase, etc.) binding site.  
Motif：  
> 1. Does not have an independent tertiary structure.  
> 2. Has specific biological functions: binding, modification, cell sublocalization, maintenance of structures, etc.  
> 3. The length is generally several to several tens of amino acids / base.      


If you want to know more about motif and motif finding problem, please see [*序列模式识别 by Prof. Xue*](http://xue.biocuckoo.org/course.html 'http://xue.biocuckoo.org/course.html "Chinese"') and [*Sequence motif - Wikipedia*](https://en.wikipedia.org/wiki/Sequence_motif 'https://en.wikipedia.org/wiki/Sequence_motif "English"')     

---
模体发现是一个经典的生物信息学问题，致力于在被同样的蛋白质结合的基因集上游序列里找到模体。  
模体（Motif）：  
> 1. 不具有独立的三级结构；  
> 2. 具有特定的生物学功能：结合、修饰、细胞亚定位和结构维持等；  
> 3. 长度通常在几个到几十个碱基/氨基酸。  


关于模体和模体发现问题，如果你想了解更多，请参阅[*序列模式识别 by Prof. Xue*](http://xue.biocuckoo.org/course.html 'http://xue.biocuckoo.org/course.html "Chinese"') 和 [*Sequence motif - Wikipedia*](https://en.wikipedia.org/wiki/Sequence_motif 'https://en.wikipedia.org/wiki/Sequence_motif "English"') 。
## Algorithm implementation

In this toturial project, we have three algorithms to solve the Motif finding problem.  
在这个教程项目里，我们将会使用三种算法来解决模体发现问题。  


| Algorithm | Code  |  
| :------- | :---------- |  
| 枚举算法——Enumeration algorithm | [*src/median_string.py*](https://github.com/ChongHui-007/Motif-finding/blob/master/src/median_string.py 'view median_string.py') [*src/motif_enumeration.py*](https://github.com/ChongHui-007/Motif-finding/blob/master/src/motif_enumeration.py 'view motif_enumeration.py') |  
| 贪心算法——Greedy algorithm | [*src/greedy_motif_search.py*](https://github.com/ChongHui-007/Motif-finding/blob/master/src/greedy_motif_search.py 'view greedy_motif_search.py') |  
| 随机算法——Randomized algorithm | [*src/gibbs_sampler.py*](https://github.com/ChongHui-007/Motif-finding/blob/master/src/gibbs_sampler.py 'view gibbs_sampler.py') [*src/randomized_motif_search.py*](https://github.com/ChongHui-007/Motif-finding/blob/master/src/randomized_motif_search.py 'view randomized_motif_search.py') |    

## Pseudocode
### MotifEnumeration
```
    MotifEnumeration(Dna, k, d)
        Patterns ← an empty set
        for each k-mer Pattern in the first string in Dna
            for each k-mer Pattern' differing from Pattern by at most d mismatches
                if Pattern' appears in each string from Dna with at most d mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns
        
        
    New_MotifEnumeration(Dna, k, d)
        Patterns ← an empty set
        for each k-mer Pattern' from AA…AA to TT…TT
            if Pattern' appears in each string from Dna with at most d mismatches
                add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns
 ```
 ```py
    def New_MotifEnumeration(Dna, k, d):
	    patterns= []                              
        mers=build_4k_k_merset(k)
		for i in mers:
	    	if multihammingd(Dna, i)<= d:
        		patterns.append(i)
        patterns= list(set(patterns))
    	return ' '.join(patterns)
```
### MedianString
```
    MedianString(Dna, k)
        distance ← ∞
        for each k-mer Pattern from AA…AA to TT…TT
            if distance > d(Pattern, Dna)
                 distance ← d(Pattern, Dna)
                 Median ← Pattern
        return Median
        
        
    New_MedianString(Dna, k)
        scores ← an empty set
        median ← an empty set
        merset ← a k-mer set of sequence from AA…AA to TT…TT
        for each k-mer Pattern from merset
            sumdistance ← sum of distances between k-mer and each sequence from Dna
            add sumdistance to scores
        for each number from 0 to length of scores
            if scores[number] equals to minimum value of scores
                add merset[numer] to median
        return median
```
```py
    def New_MedianString(dnas, k):
        merset=build_4k_k_merset(k)
        scores=[]
        median=[]
        for i in merset: 
            scores.append(sumdistance(dnas, i))
        for i in range(len(scores)):
            if scores[i]==min(scores): 
                median.append(merset[i])
        return ' '.join(median)
```
### GreedyMotifSearch  
```
    GreedyMotifSearch(Dna, k, t)
        BestMotifs ← motif matrix formed by first k-mers in each string from Dna
        for each k-mer Motif in the first string from Dna
            Motif1 ← Motif
            for i = 2 to t
                form Profile from motifs Motif1, …, Motifi - 1
                Motifi ← Profile-most probable k-mer in the i-th string in Dna
            Motifs ← (Motif1, …, Motift)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs
```
```py
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
```
### RandomizedMotifSearch  
```
    RandomizedMotifSearch(Dna, k, t)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
        BestMotifs ← Motifs
        while forever
            Profile ← Profile(Motifs)
            Motifs ← Motifs(Profile, Dna)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
            else
                return BestMotifs
```
```py
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
```
### Gibbs sampler  
```
    Gibbs sampler(Dna, k, t, N)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
        BestMotifs ← Motifs
        for j ← 1 to N
            i ← Random(t)
            Profile ← profile matrix constructed from all strings in Motifs except for Motifi
            Motifi ← profile most(Dnai, k, Profile)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs
```
```py
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
```
## Further more

* If you want to fully understand these Pseudocode, try  
    1. Learn [*Finding Hidden Messages in DNA*](https://www.coursera.org/learn/dna-analysis/home/welcome 'https://www.coursera.org/learn/dna-analysis/home/welcome')      
    2. Read [*Bioinformatics Algorithms: an Active Learning Approach*](http://bioinformaticsalgorithms.com/index.htm 'http://bioinformaticsalgorithms.com/index.htm')  
* 如果你想要完全理解这些伪代码，请尝试：  
    1. 学习在线课程[*Finding Hidden Messages in DNA*](https://www.coursera.org/learn/dna-analysis/home/welcome 'https://www.coursera.org/learn/dna-analysis/home/welcome')；   
    2. 阅读[*Bioinformatics Algorithms: an Active Learning Approach*](http://bioinformaticsalgorithms.com/index.htm 'http://bioinformaticsalgorithms.com/index.htm')。   
## To-do
1. Bug fixes
2. Add more code comments
3. Usage description
