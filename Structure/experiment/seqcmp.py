
#count number of bases that are paired
def complementarity(seq1, seq2):
	cnt = 0
	length = len(seq1)
	for i in range(0, length):
		if ((seq1[i] == 'G' and seq2[i] == 'C') or (seq1[i] == 'C' and seq2[i] == 'G') or (seq1[i] == 'T' and seq2[i] == 'A') or (seq1[i] == 'A' and seq2[i] == 'T')):
			cnt += 1
	return cnt*2

#trim sequence for complementarity calculation
def trim(seq1, seq2):

	if (len(seq1) > len(seq2)):
		seq1, seq2 = seq2, seq1

	len1 = len(seq1)
	len2 = len(seq2)

	maximum = 0
	state = 0
	pos = 0

	#state 1:ATCGATCG  				~ ATCGATCG
	#		 		ATCGATCGATCG   	  ATCGATCGATCG
	for i in range(1, len1+1):
		sub1 = seq1[-i:]
		sub2 = seq2[:i]
		comp = complementarity(sub1, sub2)
		if (comp > maximum):
			maximum = comp
			state = 1
			pos = i

	#state 2: ATCGATCG    ~     ATCGATCG
	#		 ATCGATCGATCG	ATCGATCGATCG
	for i in range(1,len2-len1+1):
		sub2 = seq2[i:i+len1]
		comp = complementarity(sub1, sub2)
		if (comp > maximum):
			maximum = comp
			state = 2
			pos = i

	#state 3:	  ATCGATCG  ~             ATCGATCG
	#		 ATCGATCGATCG	   ATCGATCGATCG
	for i in range(1,len1):
		sub1 = seq1[:-i]
		sub2 = seq2[len2-len1+i:]
		comp = complementarity(sub1, sub2)
		if (comp > maximum):
			maximum = comp
			state = 3
			pos = i

	return maximum, maximum / (len1+len2), state, pos

def main():

	Right_name = ["Rb-A'", "Rb-D_E_F", "R-D'_C_F'_B", "R-D_C'_F_B'_A'", "Rl-B_A", "Rl-X_E'"]
	Right = ["AAAAAAAACCGGTCGCTG", "AAAAAAAACACCGTACAGCCTCGTTCC", "GCCAGTGGAAGGTGCAGCTACGGTG", "CACCGTAGCTGCACCTTCCACTGGCCCGGTCGCTG", "AAAAAAAACAGCGACCGGGCCAGT", "AAAAAAAACGAGGCTGTACACGT"]

	Left_name = ["Lb-A'", "Lb-D_E_F", "L-D'_C_F'_B", "L-D_C'_F_B'_A'", "Ll-B_A", "Ll-X_E'"]
	Left = ["GCGACCTCCGAAAAAAAA", "AGCTGCGGCAGCTGAGTCCAAAAAAAA", "GGACTCCGACTGCGTAGCTCGTCAG", "GCGACCTCCGCTGACGAGCTACGCAGTCGGAGTCC", "CGTCAGCGGAGGTCGCAAAAAAAA", "ACTATCAGCTGCCGCAAAAAAAA"]

	hinge_name = ["Rb-Hinge1", "Rl-Hinge1", "Rb-Hinge2", "Rl-Hinge2", "Lb-Hinge1", "Ll-Hinge1", "Lb-Hinge2", "Ll-Hinge2"]
	hinge = ["GTCCGTGTCACAAAAAAAA", "CACTGTGCCTGAAAAAAAA", "GCATGAACACGAAAAAAAA", "GCACAAGTACGAAAAAAAA","AAAAAAAAAACAGCCTGCT", "AAAAAAAAAATCGTCCGAC", "AAAAAAAAAGACCAGTGAC", "AAAAAAAAACAGTGACCAG"]

	link_name = ["R-link1", "R-link2", "L-link1", "L-link2"]
	link = ["CAGGCACAGTGAAAAAAAAAAAAGTGACACGGAC", "CGTACTTGTGCAAAAAAAAAAAACGTGTTCATGC", "AGCAGGCTGAAAAAAAAAAAAGTCGGACGA", "GTCACTGGTCAAAAAAAAAAAACTGGTCACTG"]
	
	seq =  Right+ Left + hinge + link
	name = Right_name + Left_name + hinge_name + link_name

	#output notation
	print("[Notation]")
	print("\tl : lid")
	print("\tb : base")	
	print("\tL : left-hand side")
	print("\tR : right-hand side")
	print("\tComplementary level ＝ # of paired bases / total # of bases")
	print("==============================================================")

	for i in range(0,24):
		for j in range(i+1,24):
			maximum, percentage, state, pos = trim(seq[i], seq[j][::-1])
			
			#threshold
			if (percentage < 0.50):
				continue

			len1 = len(seq[i])
			len2 = len(seq[j])
			space1 = 0
			space2 = 0
			seq1 = {"name": name[i], "space": 0*' ', "sequence": '5\'-'+seq[i]+'-3\''}
			seq2 = {"name": name[j], "space": 0*' ', "sequence": '3\'-'+seq[j][::-1]+'-5\''}

			#restore alignment from state and pos
			if state == 1 and len1 <= len2:		seq2["space"] = (len1-pos)*' '
			elif state == 1:					seq1["space"] = (len2-pos)*' '	
			elif state == 2 and len1 <= len2:	seq1["space"] = pos*' '		
			elif state == 2:					seq2["space"] = pos*' '
			elif state == 3 and len1 <= len2:	seq1["space"] = (len2-len1+pos)*' '
			elif state == 3:					seq2["space"] = (len1-len2+pos)*' '

			#output result
			print('{name: <22}: {space}{sequence}'.format(**seq1))
			print('{name: <22}: {space}{sequence}'.format(**seq2))
			print("# of paired bases     : {}".format(maximum))
			print("Total # of bases      : {}".format(len(seq[i]+seq[j])))
			print("Complementary level  : {}".format(percentage))
			print("==============================================================")

main()