from scipy.special import binom

def is_admissible(tup):
	res = True
	for i in range(len(tup) - 1):
		if tup[i] < 2 * tup[i + 1]:
			res = False
			break
	return res


def make_admissible(tup):
	# if is_admissible(tup):
		# return {tup:1}
	admissible = True
	for i in range(len(tup) - 1):
		if tup[i] < 2 * tup[i + 1]:
			#split
			admissible = False
			break
	if admissible:
		return {tup : 1}
	res = {}
	left, right = tup[:i], tup[i+2:]
	#alter
	for c in range(int(tup[i] / 2) + 1):
		coef = binom(tup[i + 1] - c - 1, tup[i] - 2 * c)
		temp_tup = left + (tup[i + 1] + tup[i] - c, c) + right  
		# temp_tup may be not admissible
		res_temp = make_admissible(temp_tup)
		for m in res_temp:
			if m in res:
				res[m] = res[m] + coef * res_temp[m]
			else:
				res[m] = coef * res_temp[m]
		# res.update(res_temp)
		# else:
			# res[temp_tup] += res_temp[temp_tup]
	return res


def main():
	tup = input('Input sequence: ')
	print make_admissible(tup)

if __name__ == '__main__':
	main()