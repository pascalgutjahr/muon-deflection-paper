from multiprocessing import Pool
import time 


def my_function(arg):
	print('arg: ', arg)
	return arg

def main():
	start_time = time.time()
	args = [1, 2]

	pool = Pool(2)
	results = pool.map(my_function, args)
	pool.close()
	pool.join()
	end_time = time.time()
	print('results: ', results)
	print('duration: {}s'.format(end_time-start_time))

if __name__ == '__main__':
	main()
