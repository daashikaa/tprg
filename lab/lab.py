import random
import numpy as np
import math
import scipy.stats as stats
import argparse
import matplotlib.pyplot as plt

def gen10000(file_path):
	with open(file_path, 'w') as file:
		random_floats = [str(random.random()) for _ in range(10000)]
		file.write(','.join(random_floats))

def make_lst(file_path):
	with open(file_path, 'r', encoding='utf-8-sig') as f:
		data = f.read(); lst = list(map(int, data.split(","))); return lst

def read_all_files():
	ress = []
	nms = ['5p.txt', 'add.txt', 'bbs.txt', 'lc.txt', 'lfsr.txt', 'mt.txt', 'nfsr.txt', 'rc4.txt', 'rsa.txt']
	for i, nm in enumerate(nms):
		lst = make_lst(nm)
		mx = max(lst)
		lst = list(map(lambda x: x / mx, lst))
		ress.append(lst)
		print(i)
	return ress

# gen10000('random_floats.txt')
# lst = make_lst('random_floats.txt')
# lst1 = [0.0021, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

# хи-квадрат
def my_chi2(lst, a_=0.05, lst_=None, ex=None, x=None):
	if x is None: x = len(np.unique(lst))
	if lst_ is None: _, lst_ = np.unique(lst, return_counts=True)
	if ex is None: ex = np.full(x, len(lst) / x)
	chi_stat = np.sum(pow((lst_ - ex), 2) / ex); critical_value = stats.chi2.ppf(1 - a_, x - 1)
	return chi_stat <= critical_value

# print(my_chi2(lst))

# серий
def my_series(lst, a_=0.05, d=64):
	x = pow(d, 2); half = len(lst) // 2; res = np.zeros(x, dtype=int)
	for j in range(half):
		q_ = int(lst[2 * j] * d); r_ = int(lst[2 * j + 1] * d)
		q_ = min(q_, d - 1); r_ = min(r_, d - 1); res[q_ * d + r_] += 1
	exp = np.full(x, len(lst) / (2 * x))
	return my_chi2(res, a_, res, exp)

# print(my_series(lst))

# интервалов
def my_intervals(lst, n=1000, a_=0.05, a=0, b=0.5, t=100):
	lst = np.array(lst)
	j = -1; s = 0; r = 0; crs = np.zeros(t + 1, dtype=int); length = lst.shape[0]
	while s < n and j < length:
		r = 0
		while j < length and a <= lst[j] <= b: r += 1; j += 1
		crs[min(r, t)] += 1; s += 1; j += 1 
	p = b - a
	tmp = [n * p * pow((1.0 - p), r) for r in range(t)]
	exp = tmp + [n * pow((1.0 - p), t)]
	return my_chi2(lst, a_, crs, np.array(exp), t + 1)

# разбиений
def my_parts(num, a_=0.05, n=1000):
	segment_length = 10000 // n
	r = np.zeros(segment_length + 1, dtype=int)
	for i in range(n):
		segment = num[segment_length * i : segment_length * (i + 1)]
		unique_count = len(np.unique(segment))
		r[unique_count] += 1
	p = []
	scale_factor = 1.0 / 1000.0
	factorial_1000 = math.factorial(1000)
	for count in range(segment_length + 1):
		if count == 0:
				prs = 1.0
		else:
				prs = factorial_1000 / (math.factorial(1000 - count) * pow(1000,count))
		p.append(prs * scale_factor)
	lst_ = np.array([math.comb(segment_length + index - 1, index) / pow(1000, segment_length) for index in range(segment_length + 1)])
	return my_chi2(num, a_, lst_[1:], np.array(p[1:]), segment_length)

# перестановок
def my_perms(lst, a_=0.05, t=20):
	n = len(lst); x = math.factorial(t); group_counts = {}
	for i in range(0, n, t):
		group = tuple(sorted(lst[i:i + t]))
		if group in group_counts: group_counts[group] += 1
		else: group_counts[group] = 1
	lst_ = list(group_counts.values()); exp = np.array([n / x] * len(lst_))
	return my_chi2(lst, a_, lst_, exp, x)

a_mtrx = [[4529.4, 9044.9, 13568, 22615,  22615,  27892 ],
					[9044.9, 18097,  27139, 36187,  452344, 55789 ],
					[13568,  27139,  40721, 54281,  67582,  83685 ],
					[18091,  36187,  54281, 72414,  90470,  111580],
					[22615,  45234,  67852, 90470,  113262, 139476],
					[27892,  55789,  83685, 111580, 139476, 172860]]

b_vec = [1 / 6, 5 / 24, 11 / 120, 19 / 720, 29 / 5040, 1 / 840]


# монотонности
def get_series_lens(lst):
	res = []
	i = 0
	while i < len(lst):
		lenl = 1
		while i + lenl < len(lst) and lst[i + lenl - 1] <= lst[i + lenl]:
			lenl += 1
		res.append(lenl)
		i += lenl
	return res

def calculate_expected_values(lst, series_lens):
	n = len(lst); res = []
	tmp = 0
	for lenl in series_lens:
		m = 1/6; min_val = min(lenl, 6)
		for i in range(min_val):
			for j in range(min_val):
				m += (lst[i + tmp] - n * b_vec[i]) * (lst[j + tmp] - n * b_vec[j]) * a_mtrx[i][j]
		tmp += lenl; res.append(m)
	return res

def my_mono(lst, a_=0.05):
	series_lengths = get_series_lens(lst)
	res = calculate_expected_values(lst, series_lengths)
	return my_chi2(res, a_)

def plot_exp_vs_sample_size(exps):
	sample_sizes = [len(exps[:i]) for i in range(len(exps)) ]
	if len(sample_sizes) != len(exps):
		raise ValueError("Длины списков sample_sizes и expectations должны быть равны.")
	plt.plot(sample_sizes, exps, marker='o', markersize=1, linestyle='-', linewidth=1, color='red', markerfacecolor='red', markeredgecolor='red')
	plt.xlabel('Объем выборки (x)')
	plt.ylabel('Математическое ожидание (y)')
	plt.title('График зависимости математического ожидания от объема выборки')
	plt.grid(True)
	plt.show()

def plot_std_vs_sample_size(std):
	sample_sizes = [len(std[:i]) for i in range(len(std)) ]
	if len(sample_sizes) != len(std):
		raise ValueError("Длины списков sample_sizes и expectations должны быть равны.")
	plt.plot(sample_sizes, std, marker='o', markersize=1, linestyle='-', linewidth=1)
	plt.xlabel('Объем выборки (x)')
	plt.ylabel('Cр-кв. отклонение (y)')
	plt.title('График зависимости ср-кв. отклонения от объема выборки')
	plt.grid(True)
	plt.show()

# конфликтов
def my_confs(lst):
	n = len(lst); m = pow(2, 20); avr = n / m; x = 1 - n / m + math.factorial(n) / (2 * math.factorial(n - 2) * pow(m, 2)); _x_ = n / m - 1 + x
	return not(_x_ < avr + 10 or _x_ > avr - 10)

# матожидание и среднеквадратическое отклонение
def math_exp_and_std(num_lst):
	math_exp = num_lst.mean(); std = num_lst.std(); return math_exp, std

def absolute(lst_):
	math_const = 0.4993; std_const = 0.2881; math_exp, std = math_exp_and_std(lst_)
	print("Отн. погрешность мат. ожидания: ", abs(math_const - math_exp))
	print("Отн. погрешность ср-кв. отклонения: ", abs(std_const - std))

def lsts_for_dia(lst_):
	res_m = []
	res_s = []
	for i in range(len(lst_)):
		if len(lst_[:i]) > 0:
			res_m.append(np.mean(lst_[:i]))
			res_s.append(np.std(lst_[:i]))
	return res_m, res_s

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	help1 = "Имя файла с входной последовательностью"
	help2 = """ Критерии для проверки заданной последовательности на равномерное распределение:
								->  chi - Критерий хи-квадрат,
								->  ser - Критерий серий,
								->  int - Критерий интервалов,
								->  par - Критерий разбиений,
								->  per - Критерий перестановок,
								->  mon - Критерий монотонности,
								->  con - Критерий конфликтов,
								->  all - Все критерии
				"""
	types_crit = ["chi", "ser", "int", "par", "per", "mon", "con", "all"]

	parser.add_argument("-t", nargs=1, choices=types_crit, required=True, default=[None], help=help2)
	parser.add_argument("-f", nargs=1, required=False, default=[None], help=help1)

	args = parser.parse_args()
	type = args.t[0]
	path = args.f[0]
	# llsstt = read_all_files()
	# for k in range(len(llsstt)):
	# 	print("Критерий хи-квадрат:")
	# 	print("-> ", my_chi2(llsstt[k]))
	# 	print("Критерий серий:")
	# 	print("-> ", my_series(llsstt[k]))
	# 	print("Критерий интервалов:")
	# 	print("-> ", my_intervals(llsstt[k]))
	# 	print("Критерий разбиений:")
	# 	print("-> ", my_parts(llsstt[k]))
	# 	print("Критерий перестановок:")
	# 	print("-> ", my_perms(llsstt[k]))
	# 	print("Критерий монотонности:")
	# 	print("-> ", my_mono(llsstt[k]))
	# 	print("Критерий конфликтов:")
	# 	print("-> ", my_confs(llsstt[k]))
	if path == None:
		gen10000("filegen10000.txt")
		lstt = np.array(make_lst("filegen10000.txt"))
	else:
		lst1 = make_lst(path)
		mx = max(lst1)
		lstt = np.array(list(map(lambda x: x / mx, lst1)))
	v1, v2 = math_exp_and_std(lstt); print("Мат. ожидание: ", v1); print("Cр-кв. отклонение: ", v2); absolute(lstt)

	try:
		match type:
				case 'chi':
						print("Критерий хи-квадрат:")
						print("-> ", my_chi2(lstt))
				case 'ser':
						print("Критерий серийт:")
						print("-> ", my_series(lstt))
				case 'int':
						print("Критерий интервалов:")
						print("-> ", my_intervals(lstt))
				case 'par':
						print("Критерий разбиений:")
						print("-> ", my_parts(lstt))
				case 'per':
						print("Критерий перестановок:")
						print("-> ", my_perms(lstt))
				case 'mon':
						print("Критерий монотонности:")
						print("-> ", my_mono(lstt))
				case 'con':
						print("Критерий конфликтов:")
						print("-> ", my_confs(lstt))
				case 'all':
						print("Критерий хи-квадрат:")
						print("-> ", my_chi2(lstt))
						print("Критерий серий:")
						print("-> ", my_series(lstt))
						print("Критерий интервалов:")
						print("-> ", my_intervals(lstt))
						print("Критерий разбиений:")
						print("-> ", my_parts(lstt))
						print("Критерий перестановок:")
						print("-> ", my_perms(lstt))
						print("Критерий монотонности:")
						print("-> ", my_mono(lstt))
						print("Критерий конфликтов:")
						print("-> ", my_confs(lstt))
	except Exception as err:
			print("В процессе проверки произошла ошибка!")
			print(" -> " + str(err))

	mm, st = lsts_for_dia(lstt)
	plot_exp_vs_sample_size(mm); plot_std_vs_sample_size(st)