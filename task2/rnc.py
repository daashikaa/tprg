import sys
import math

"""стандартное равномерное распределение"""
def st(a, b, lst):
	res = []
	X_ = max(lst)
	m = X_ + 1
	for x in lst:
		U = x / m
		Y = b * U + a
		res.append(Y)
	return res

"""треугольное распределение"""
def tr(a, b, lst):
	res = []
	X_ = max(lst)
	m = X_ + 1
	lst = [x / m for x in lst]
	lst.append(lst[-1])
	for i in range(len(lst) - 1):
		U1 =lst[i]
		U2 = lst[i + 1]
		Y = a + b * (U1 + U2 - 1)
		res.append(Y)
	return res

"""общее экспоненциальное распределение"""
def ex(a, b, lst):
	res = []
	X_ = max(lst)
	m = X_ + 1
	for x in lst:
		U = x / m
		Y = -b * math.log(U) + a
		res.append(Y)
	return res

"""нормальное распределение"""
def nr(mu, sigma, lst):
	res = []
	X_ = max(lst)
	m = X_ + 1
	lst = [x / m for x in lst]
	lst.append(lst[-1])
	for i in range(0, len(lst) - 1, 2):
		U1 =lst[i]
		U2 = lst[i + 1]
		Z1 = mu + sigma * math.sqrt(-2 * math.log(1 - U1)) * math.cos(2 * math.pi * U2)
		Z2 = mu + sigma * math.sqrt(-2 * math.log(1 - U1)) * math.sin(2 * math.pi * U2)
		res.append(Z1)
		res.append(Z2)
	return res

"""гамма распределение"""
def gm(a, b, k, lst):
	res = []
	X_ = max(lst)
	m = X_ + 1
	lst = [x / m for x in lst]
	uk = []
	for i in range(0, len(lst), k):
		uk.append(lst[i : i + k])
	if len(uk[-1]) != k: uk.pop()
	for luk in uk:
		for_log = 1
		for u in luk:
			for_log = for_log * (1 - u)
		Y = a - b * math.log(for_log)
		res.append(Y)
	return res

"""логнормальное распределение"""
def ln(a, b, lst):
	res = []
	lst = nr(0, 1, lst)
	for Z in lst:
		Y = a + math.exp(b - Z)
		res.append(Y)
	return res

"""логистическое распределение"""
def ls(a, b, lst):
	res = []
	X_ = max(lst)
	m = X_ + 1
	for x in lst:
		U = x / m
		Y = a + b * math.log(U / (1 - U))
		res.append(Y)
	return res

"""биномиальное распределение"""
def bi(p, p2, lst):
	res = []
	X_ = max(lst)
	m = X_ + 1
	lst = [x / m for x in lst]
	n = len(lst)
	for i in range(0, n):
		sum = 0
		y = lst[i]
		while sum < lst[i]:
			k = 0
			while y > k:
				sum += math.comb(n, k) + pow(p, k) + pow(1 - p, n - k)
				k = k + 1
			y = y + 1
		t = 1
		for _ in range(p2):
			t *= 10
		res.append(int(y *  10**23) % t)
	return res

def print_usage():
	print("Применение: rnc.exe [/f:<filename>] /d:<distribution> /p1:<param1> /p2:<param2> [/p3:<param3>]")
	print('       rnc.exe /h for help')

def print_help():
	print_usage()
	print('Преобразование ПСЧ к заданному распределению')
	print()
	print('Возможные аргументы:')
	print('  /f: <filename>   полное имя файла, из которого будет браться начальная последовательность (по умолчанию lst.dat)')
	print('  /d: <distribution>  код распределения для преобразования последовательности:')
	print('                      st – стандартное равномерное с заданным интервалом;')
	print('                      tr – треугольное распределение;')
	print('                      ex – общее экспоненциальное распределение;')
	print('                      nr – нормальное распределение;')
	print('                      gm – гамма распределение;')
	print('                      ln – логнормальное распределение;')
	print('                      ls – логистическое распределение;')
	print('                      rsa – ГПСЧ на основе RSA;')
	print('                      bi – биномиальное распределение.')
	print('  /p1:<parameter1>      1-й параметр, необходимый, для генерации ПСЧ заданного распределения')
	print('  /p1:<parameter2>      2-й параметр, необходимый, для генерации ПСЧ заданного распределения')
	print('  /p1:<parameter3>      3-й параметр, необходимый, для генерации ПСЧ заданного распределения')


def main():
	args = sys.argv[1:]
	filename = 'lst.dat'
	distribution = None
	params = {}

	for arg in args:
		if arg.startswith("/f:"):
			filename = arg[3:]
		elif arg.startswith("/d:"):
			distribution = arg[3:]
		elif arg.startswith("/p"):
			param_name, param_value = arg[1:].split(":")
			params[param_name] = int(param_value)
		elif arg == '/h':
			if arg != args[0]:
				print('Ошибка: /h не может быть использован одновременно с другими аргументами!')
				print_usage()
				return
			else:
				print_help()
				return
			
	if distribution is None:
		print('Метод не указан!')
		print_usage()
		return
	
	if params is None:
		print('Неуказаны аргументы для функций!')
		print_usage()
		return

	if filename is None or distribution is None:
		print("Применение: rnc.exe /f:<filename> /d:<distribution> /p1:<param1> /p2:<param2> /p3:<param3>")
		return

	try:
		with open(filename, "r") as file:
			numbers = file.readline().strip().split(',')
			lst = [int(num) for num in numbers]
	except FileNotFoundError:
		print(f"Файл '{filename}' не найден!")
		return
 
	transformed_numbers = []
	# Apply distribution transformation
	if distribution == "st":
		transformed_numbers = st(params["p1"], params["p2"], lst)
	elif distribution == "tr":
			transformed_numbers = tr(params["p1"], params["p2"], lst)
	elif distribution == "ex":
			transformed_numbers = ex(params["p1"], params["p2"], lst)
	elif distribution == "nr":
			transformed_numbers = nr(params["p1"], params["p2"], lst)
	elif distribution == "gm":
			transformed_numbers = gm(params["p1"], params["p2"], params["p3"], lst)
	elif distribution == "ln":
			transformed_numbers = ln(params["p1"], params["p2"], lst)
	elif distribution == "ls":
			transformed_numbers = ls(params["p1"], params["p2"], lst)
	elif distribution == "bi":
			transformed_numbers = bi(params["p1"], params["p2"], lst)
	else:
		print('Введенный метод не предусмотрен или не существует!')
		return

	output_filename = f"distr-{distribution}.dat"

	with open(output_filename, 'w', encoding = 'UTF-8') as f:
		f.write(','.join(map(str, transformed_numbers)))
		print(f'Числа сгенерированы и сохранены в файл {output_filename}!')


if __name__ == '__main__':
	main()
