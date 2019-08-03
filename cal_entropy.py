import math



def calentropy(dict):
	entropy={}
	for i in dict.keys():
		entropy[i]=0
		for j in dict[i]:
			if j!=0: entropy[i]-=(j*(math.log(j, 2)))
	print(entropy)


calentropy({'A':[0.5, 0, 0, 0.5],'B':[0.25, 0.25, 0.25, 0.25],'C':[0, 0, 0, 1],'D':[0.25, 0, 0.5, 0.25]})
