from motif_enumeration import *


def sumdistance(l, i):
	d=[] 
	for j in l: d.append(hammingd(j, i))
	return sum(d)



def medianstring(dnas, k):
	d=len(dnas[0]); merset=build_4k_k_merset(k); s=[]; median=[]
	for i in merset: 
		s.append(sumdistance(dnas, i))
	for i in range(len(s)):
		if s[i]==min(s): median.append(merset[i])
	return ' '.join(median)


print(medianstring(['CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC','GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC','GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG'], 7))	
