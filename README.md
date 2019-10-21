# Motif-finding
---

Motif finding problem is a classical bioinformatics problem, aiming to quickly find a series of motifs on genes with the same enzyme (DNA replicase, etc.) binding site.  
Motif：  
> 1. Does not have an independent tertiary structure.  
> 2. Has specific biological functions: binding, modification, cell sublocalization, maintenance of structures, etc.  
> 3. The length is generally several to several tens of amino acids / base.      


If you want to know more about motif and motif finding problem, please see [*序列模式识别 by Prof. Xue*](http://xue.biocuckoo.org/course.html 'Chinese') and [*Sequence motif - Wikipedia*](https://en.wikipedia.org/wiki/Sequence_motif 'English')   

## Algorithm implementation  

In this toturial project, we have three algorithms to solve the Motif finding problem.  


| Algorithm | Code  |  
| :------- | :---------- |  
| Enumeration algorithm | [*src/median_string.py*](https://github.com/ChongHui-007/Motif-finding/blob/master/src/median_string.py 'view median_string.py') [*src/motif_enumeration.py*](https://github.com/ChongHui-007/Motif-finding/blob/master/src/motif_enumeration.py 'view motif_enumeration.py') |  
| Greedy algorithm | [*src/greedy_motif_search.py*](https://github.com/ChongHui-007/Motif-finding/blob/master/src/greedy_motif_search.py 'view greedy_motif_search.py') |  
| Randomized algorithm | [*src/gibbs_sampler.py*](https://github.com/ChongHui-007/Motif-finding/blob/master/src/gibbs_sampler.py 'view gibbs_sampler.py') [*src/randomized_motif_search.py*](https://github.com/ChongHui-007/Motif-finding/blob/master/src/randomized_motif_search.py 'view randomized_motif_search.py') |    

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
    1. Learn [*Finding Hidden Messages in DNA*](https://www.coursera.org/learn/dna-analysis/home/welcome)      
    2. Read [*Bioinformatics Algorithms: an Active Learning Approach*](http://bioinformaticsalgorithms.com/index.htm)  
## To-do list
1. Bug fixes
2. Add more code comments
3. Usage description
