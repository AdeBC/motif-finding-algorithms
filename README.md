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

```
    MotifEnumeration(Dna, k, d)
        Patterns ← an empty set
        for each k-mer Pattern in the first string in Dna
            for each k-mer Pattern’ differing from Pattern by at most d mismatches
                if Pattern' appears in each string from Dna with at most d mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns
```

```
    MedianString(Dna, k)
        distance ← ∞
        for each k-mer Pattern from AA…AA to TT…TT
            if distance > d(Pattern, Dna)
                 distance ← d(Pattern, Dna)
                 Median ← Pattern
        return Median
```
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
## Further more

If you want to fully understand these Pseudocode, try  
    1. Learn [*Finding Hidden Messages in DNA*](https://www.coursera.org/learn/dna-analysis/home/welcome 'https://www.coursera.org/learn/dna-analysis/home/welcome')      
    2. Read [*Bioinformatics Algorithms: an Active Learning Approach*](http://bioinformaticsalgorithms.com/index.htm 'http://bioinformaticsalgorithms.com/index.htm')  
如果你想要完全理解这些伪代码，请尝试：  
    1. 学习在线课程[*Finding Hidden Messages in DNA*](https://www.coursera.org/learn/dna-analysis/home/welcome 'https://www.coursera.org/learn/dna-analysis/home/welcome')；   
    2. 阅读[*Bioinformatics Algorithms: an Active Learning Approach*](http://bioinformaticsalgorithms.com/index.htm 'http://bioinformaticsalgorithms.com/index.htm')。   
## To-do
1. Bug fixes
2. Add more code comments
3. Usage description
