import time
import multiprocessing

def poolTest():
    for _ in range(5):
        starttime = time.time()
        pool = multiprocessing.Pool()
        pool.map(multiprocessing_func, range(0, 10))
        pool.close()
        print('That took {} seconds'.format(time.time() - starttime))

def basic_func(x):
    if x == 0:
        return 'zero'
    elif x%2 == 0:
        return 'even'
    else:
        return 'odd'

def multiprocessing_func(x):
    y = x*x
    time.sleep(0.2)
    print('{} squared results in a/an {} number'.format(x, basic_func(y)))