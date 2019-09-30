# Motif-finding
## 模体发现
模体发现是经典的生物信息学问题之一，目的在于快速地找到一系列具有相同酶（DNA复制酶等）作用位点的基因上的模体。   
以下关于模体的解释来自华中科技大学生命科学与技术学院教授薛宇老师制作的[生物信息学课件](http://xue.biocuckoo.org/course.html)第六章：序列模式识别  
![Motif](https://github.com/ChongHui-007/Motif-finding/blob/master/src/motif.png)
## 算法实现
总共有三种算法解决Motif finding问题：
>枚举算法-Enumeration algorithm  
贪心算法-Greedy algorithm  
随机算法-Randomized algorithm    

**项目文件中的src/median_string.py和src/motif_enumeration.py即是通过枚举算法解决了模体发现问题**  
**而src/greedy_motif_search.py是通过贪心算法解决了模体发现问题**  
**最后，src/gibbs_sampler.py和src/randomized_motif_search.py都是通过随机算法解决了模体发现问题**  

## 伪代码
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
如果你想要理解伪代码本身，建议你：  
学习在线课程[Finding Hidden Messages in DNA (Bioinformatics I)](https://www.coursera.org/learn/dna-analysis/home/welcome)  
或阅读书籍[《Bioinformatics Algorithms: an Active Learning Approach》](http://bioinformaticsalgorithms.com/index.htm)
