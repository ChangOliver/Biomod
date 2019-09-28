
def affinity(seq1, seq2):
	cnt = 0
	length = len(seq1)
	for i in range(0, length):
		if ((seq1[i] == 'G' and seq2[i] == 'C') or (seq1[i] == 'C' and seq2[i] == 'G') or (seq1[i] == 'T' and seq2[i] == 'A') or (seq1[i] == 'A' and seq2[i] == 'T')):
			cnt += 1
	return cnt*2

def trim(seq1, seq2):

	if (len(seq1) > len(seq2)):
		seq1, seq2 = seq2, seq1

	len1 = len(seq1)
	len2 = len(seq2)

	maximum = 0
	state = 0
	pos = 0
	for i in range(1, len1+1):
		sub1 = seq1[-i:]
		sub2 = seq2[:i]
		affi = affinity(sub1, sub2)
		if (affi > maximum):
			maximum = affi
			state = 1
			pos = i

	for i in range(1,len2-len1+1):
		sub2 = seq2[i:i+len1]
		affi = affinity(sub1, sub2)
		if (affi > maximum):
			maximum = affi
			state = 2
			pos = i

	for i in range(1,len1):
		sub1 = seq1[:-i]
		sub2 = seq2[len2-len1+i:]
		affi = affinity(sub1, sub2)
		if (affi > maximum):
			maximum = affi
			state = 3
			pos = i

	return maximum, maximum / (len1+len2), state, pos

def main():

	Top_name = ["Tb-A'", "Tb-D_E_F", "T-D'_C_F'_B", "T-D_C'_F_B'_A'", "Tl-B_A", "Tl-X_E'"]
	Top = ["AAAAAAAACCGGTCGCTG", "AAAAAAAACACCGTACAGCCTCGTTCC", "GCCAGTGGAAGGTGCAGCTACGGTG", "CACCGTAGCTGCACCTTCCACTGGCCCGGTCGCTG", "AAAAAAAACAGCGACCGGGCCAGT", "AAAAAAAACGAGGCTGTACACGT"]

	Bottom_name = ["Bb-A'", "Bb-D_E_F", "B-D'_C_F'_B", "B-D_C'_F_B'_A'", "Bl-B_A", "Bl-X_E'"]
	Bottom = ["GCGACCTCCGAAAAAAAA", "AGCTGCGGCAGCTGAGTCCAAAAAAAA", "GGACTCCGACTGCGTAGCTCGTCAG", "GCGACCTCCGCTGACGAGCTACGCAGTCGGAGTCC", "CGTCAGCGGAGGTCGCAAAAAAAA", "ACTATCAGCTGCCGCAAAAAAAA"]

	hinge_name = ["TbR-Hinge", "TlR-Hinge", "TbL-Hinge", "TlL-Hinge", "BbL-Hinge", "BbR-Hinge", "BlL-Hinge", "BlR-Hinge"]
	hinge = ["GTCCGTGTCACAAAAAAAA", "CACTGTGCCTGAAAAAAAA", "GCATGAACACGAAAAAAAA", "GCACAAGTACGAAAAAAAA", "AAAAAAAAAGACCAGTGAC", "AAAAAAAAAACAGCCTGCT", "AAAAAAAAACAGTGACCAG", "AAAAAAAAAATCGTCCGAC"]

	link_name = ["TR-link", "TL-link", "BR-link", "BL-link"]
	link = ["CAGGCACAGTGAAAAAAAAAAAAGTGACACGGAC", "CGTACTTGTGCAAAAAAAAAAAACGTGTTCATGC", "AGCAGGCTGAAAAAAAAAAAAGTCGGACGAA", "GTCACTGGTCAAAAAAAAAAAACTGGTCACTG"]
	
	precip_name = ["Tl-precip", "Tl-precip'", "Bl-precip", "Bl-precip'"]
	precip = ["AAACAACACAACACAACAC", "GTGTTGTGTTGTGTTG", "AAAGAAAGAGAAAGAGAAAG", "CTTTCTCTTTCTCTTTC"]
	
	seq = Top + Bottom + hinge + link + precip
	name = Top_name + Bottom_name + hinge_name + link_name + precip_name

	print("[Notation]")
	print("\tl : lid")
	print("\tb : base")	
	print("\tT : top")
	print("\tB : bottom")
	print("\tL : left-hand side")
	print("\tR : right-hand side")
	print("\tAffinity level ＝ # of paired bases / total # of bases")
	print("==============================================================")

	for i in range(0,28):
		for j in range(i+1,28):
			maximum, percentage, state, pos = trim(seq[i], seq[j][::-1])
			#if (percentage < 0.50):
			#	continue
			len1 = len(seq[i])
			len2 = len(seq[j])
			space1 = 0
			space2 = 0
			seq1 = {"name": name[i], "space": 0*' ', "sequence": '5\'-'+seq[i]+'-3\''}
			seq2 = {"name": name[j], "space": 0*' ', "sequence": '3\'-'+seq[j][::-1]+'-5\''}

			if state == 1 and len1 <= len2:		seq2["space"] = (len1-pos)*' '
			elif state == 1:					seq1["space"] = (len2-pos)*' '	
			elif state == 2 and len1 <= len2:	seq1["space"] = pos*' '		
			elif state == 2:					seq2["space"] = pos*' '
			elif state == 3 and len1 <= len2:	seq1["space"] = (len2-len1+pos)*' '
			elif state == 3:					seq2["space"] = (len1-len2+pos)*' '

			print('{name: <22}: {space}{sequence}'.format(**seq1))
			print('{name: <22}: {space}{sequence}'.format(**seq2))
			print("# of paired bases     : {}".format(maximum))
			print("Total # of bases      : {}".format(len(seq[i]+seq[j])))
			print("Affinity level  : {}".format(percentage))
			print("==============================================================")

main()