def hammingd(seq1, seq2):
	c=0
	for i in range(len(seq1)):
		if seq1[i]!=seq2[i]: c+=1
	print(c)


hammingd('CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT','CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG')
