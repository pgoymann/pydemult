def my_data_processor(x):
    # some expensive computation
    print(x)
    x+=1
    return x

def my_data_generator():
    for i in range(20):
        yield i
from mputil import lazy_map

gen = my_data_generator()
print(gen)
print(lazy_map(data_processor=my_data_processor, 
               data_generator=gen, 
               n_cpus=5))