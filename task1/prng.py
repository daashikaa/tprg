import sys
"""линейный конгруэнтный метод"""
def lc (m, a, c, x0, n_ = 10000):
	x = []
	x.append(x0);
	for i in range(1, n_ - 1):
		x.append(((a * x[i - 1] + c) % m) % 2**10)
	return x

"""аддитивный метод"""
def add(m, low_i, up_i, start_seq, n_ = 10000):
	x = start_seq
	seq = start_seq.copy()
	n = len(seq)
	l = n_ - n
	for _ in range(l):
		xn = (seq[n - low_i] + seq[n - up_i]) % m
		x.append(xn % 2**10)
		seq.append(xn)
		del seq[0]
	return x

"""5p"""
def _5p(p, q1, q2, q3, w, x0, n_ = 10000):
	x = []
	x0 = int(x0, 2)
	for _ in range(n_):
		cur = 0
		for _ in range(w):
			bit_q1 = (x0 >> p - q1) & 1
			bit_q2 = (x0 >> p - q2) & 1
			bit_q3 = (x0 >> p - q3) & 1
			bit_x0 = x0 & 1
			xor = bit_q1 ^ bit_q2 ^ bit_q3 ^ bit_x0
			cur = (cur << 1) | xor
			x0 = (x0 >> 1) | (xor << p - 1)
		x.append(cur % 2**10)
	return x

"""рслос"""
def lfsr(vec_bin, start_reg, n_= 10000):
	x = []
	reg_len = len(start_reg)
	vec_bin = int(vec_bin, 2)
	start_reg = int(start_reg, 2)
	shift = reg_len - 1
	for _ in range(n_):
		new = (start_reg ^ vec_bin >> shift) & 1
		start_reg = (start_reg >> 1) | (new << (shift))
		x.append(start_reg % 2**10)
	return x

"""нелинейная комбинация рслос"""
def nfsr(R1, R2, R3, w, x1, x2, x3, n_ = 10000):
	x = []
	lR1 = len(R1)
	lR2 = len(R2)
	lR3 = len(R3)
	for _ in range(n_):
		cur = 0
		for _ in range(w):
			xorR1 = (x1 ^ (x1 >> lR1 - 1)) 
			xorR2 = (x2 ^ (x2 >> lR1 - 1))
			xorR3 = (x3 ^ (x3 >> lR1 - 1))
			res = ((xorR1 ^ xorR2) + (xorR2 ^ xorR3) + xorR3) & 1
			x1 = (x1 >> 1) | (res << lR1 - 1)
			x2 = (x2 >> 1) | (res << lR2 - 1)
			x3 = (x3 >> 1) | (res << lR3 - 1)
			cur = (cur << 1) | res
		x.append(cur % 2**10)
	return x

"""вихрь мерсенна"""
def mt(x0, n = 624, n_= 10000):
	w = 32
	r = 31
	m = 397
	a = 0x9908B0DF
	u = 11
	d = 0xFFFFFFFF
	s = 7
	t = 15
	l = 18
	b = 0x9D2C5680
	c = 0xEFC60000
	f = 1812433253
	lst = [0] * n
	ind = n + 1

	def f1(core):
		lst[0] = core
		for i in range(1, n):
			tmp = f * (lst[i - 1] ^ (lst[i - 1] >> (w - 2))) + i
			lst[i] = tmp & d

	def f2():
		nonlocal ind
		if ind >= n:
			f3()
			ind = 0
		y = lst[ind]
		y = y ^ ((y >> u) & d)
		y = y ^ ((y << s) & b)
		y = y ^ ((y << t) & c)
		y = y ^ (y >> l)
		ind = ind + 1
		return y & d

	def f3():
		for i in range(n):
			x = (lst[i] >> r) + (lst[(i + 1) % n] & ((1 << r) - 1))
			x_a = x >> 1
			lst[i] = lst[(i + m) % n] ^ x_a
			if (x % 2) != 0:
				lst[i] ^= a

	f1(x0)
	res = []
	for _ in range(n_):
		o = f2()
		res.append(o % 2**10)
	return res

"""rc4"""
def rc4(k, count_num = 10000):
	len_k = len(k)
	S = [i for i in range(256)]
	j = 0
	for i in range(256):
			j = (j + S[i] + k[i % len_k]) % 256
			S[i], S[j] = S[j], S[i]
	i = 0
	j = 0
	Ks = []
	for _ in range(1, count_num + 1):
		i = (i + 1) % 256
		j = (j + S[i]) % 256
		S[i], S[j] = S[j], S[i]
		cur = S[(S[i] + S[j]) % 256]
		Ks.append(cur % 2**10)
	return Ks

""""гпсч на основе rsa"""
def rsa(n, e, w, x0, n_ = 10000):
	x = []
	for _ in range(n_):
		cur = 0
		for _ in range(w):
			x0 = pow(x0, e, n)
			cur = (cur << 1) | (x0 & 1)
		x.append(cur % 2**10)
	return x

"""алгоритм блюма-блюма-шуба"""
def bbs(x0, l, n_ = 10000):
	p = 127
	q = 131
	n = p * q
	x = []
	for _ in range(n_):
		cur = 0
		for _ in range(l):
			x0 = pow(x0, 2, n)
			cur = (cur << 1) | (x0 & 1)
		x.append(cur % 2**10)
	return x

def generate_pseudo_random(method, args, n):
	if method == 'lc':
		int_values = [int(x) for x in args.split(',')]
		if len(int_values) != 4:
			print('Было передано неверное число аргументов!')
			return []
		m, a, c, x0 = int_values
		return lc(m, a, c, x0, n)

	elif method == 'add':
		args_split = args.split(',')
		if len(args_split) < 4:
			print('Было передано неверное число аргументов!')
			return []
		m = int(args_split[0])
		low_i = int(args_split[1])
		up_i = int(args_split[2])
		start_seq = list(map(int, args_split[3:]))
		return add(m, low_i, up_i, start_seq, n)

	elif method == '5p':
		args_split = args.split(',')
		if len(args_split) != 6:
			print('Было передано неверное число аргументов!')
			return []
		p = int(args_split[0])
		q1 = int(args_split[1])
		q2 = int(args_split[2])
		q3 = int(args_split[3])
		w = int(args_split[4])
		x0 = args_split[5]
		return _5p(p, q1, q2, q3, w, x0, n)

	elif method == 'lfsr':
		args_split = args.split(',')
		if len(args_split) != 2:
			print('Было передано неверное число аргументов!')
			return []
		vec_bin = args_split[0]
		start_reg = args_split[1]
		return lfsr(vec_bin, start_reg, n)

	elif method == 'nfsr':
		args_split = args.split(',')
		if len(args_split) != 7:
			print('Было передано неверное число аргументов!')
			return []
		R1 = args_split[0]
		R2 = args_split[1]
		R3 = args_split[2] 
		w = int(args_split[3])
		x1 = int(args_split[4])
		x2 = int(args_split[5])
		x3 = int(args_split[6])
		return nfsr(R1, R2, R3, w, x1, x2, x3, n)

	elif method == 'mt':
		args_split = args.split(',')
		if len(args_split) != 1:
			print('Было передано неверное число аргументов!')
			return []
		x0 = int(args_split[0])
		return mt(x0, 624, n)

	elif method == 'rc4':
		args_split = args.split(',')
		if len(args_split) != 256:
			print('Было передано неверное число аргументов!')
			return []
		k = list(map(int, args_split))
		return rc4(k, n)

	elif method == 'rsa':
		args_split = args.split(',')
		if len(args_split) != 4:
			print('Было передано неверное число аргументов!')
			return []
		m, e, w, x0 = map(int, args_split)
		return rsa(m, e, w, x0, n)

	elif method == 'bbs':
		args_split = args.split(',')
		if len(args_split) != 2:
			print('Было передано неверное число аргументов!')
			return []
		x0 = int(args_split[0])
		l = int(args_split[1])
		return bbs(x0, l, n)

	else:
		print("Введенный метод не предусмотрен или не существует!")
		return []

def print_usage():
	print('Применение: prng.exe /g:<generator> [/n:<count>] [/f:<filename>] [/i:<params>]')
	print('       prng.exe /h for help')

def print_help():
	print_usage()
	print('Генерация последовательности псевдослучайных чисел по выбранному методу')
	print()
	print('Возможные аргументы:')
	print('  /g: <generator>  параметр указывает на метод генерации ПСЧ:')
	print('                      lc  – линейный конгруэнтный метод;')
	print('                      add – аддитивный метод;')
	print('                      5p  – пятипараметрический метод;')
	print('                      lfsr – регистр сдвига с обратной связью (РСЛОС);')
	print('                      nfsr – нелинейная комбинация РСЛОС;')
	print('                      mt  – вихрь Мерсенна;')
	print('                      rc4 – RC4;')
	print('                      rsa – ГПСЧ на основе RSA;')
	print('                      bbs – алгоритм Блюма-Блюма-Шуба;')
	print('  /n: <count>      количество генерируемых чисел (по умолчанию 10000)')
	print('  /f: <filename>   полное имя файла, в который будут выводиться данные (по умолчанию rnd.dat)')
	print('  /i: <params>     перечисление параметров для выбранного генератора')

def main():

	args = sys.argv[1:]
	method = None
	n = 10000
	filename = 'rnd.dat'

	for arg in args:
		if arg.startswith('/g:'):
			method = arg[3:]
		elif arg.startswith('/n:'):
			n = int(arg[3:])
		elif arg.startswith('/f:'):
			filename = arg[3:]
		elif arg.startswith("/i:"):
			input_args = arg[3:]
		elif arg == '/h':
			if arg != args[0]:
				print('Ошибка: /h не может быть использован одновременно с другими аргументами!')
				print_usage()
				return
			else:
				print_help()
				return

	if method is None:
		print('Метод не указан!')
		print_usage()
		return

	elif input_args is None:
		print('Неуказаны аргументы для функций!')
		print_usage()
		return

	random_numbers = generate_pseudo_random(method, input_args, n)

	if random_numbers == []:
		return
	else:
		with open(filename, 'w', encoding = 'UTF-8') as f:
			f.write(','.join(map(str, random_numbers)))
		print(f'Числа сгенерированы и сохранены в файл {filename}!')


if __name__ == '__main__':
	main()