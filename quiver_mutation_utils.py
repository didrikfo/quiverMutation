import itertools

def list_intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def sublist_exists(list, sublist):
    for i in range(len(list)-len(sublist)+1):
        if sublist == list[i:i+len(sublist)]:
            return True #return position (i) if you wish
    return False

def divisors(n):
    # get factors and their counts
    factors = {}
    nn = n
    i = 2
    while i*i <= nn:
        while nn % i == 0:
            factors[i] = factors.get(i, 0) + 1
            nn //= i
        i += 1
    if nn > 1:
        factors[nn] = factors.get(nn, 0) + 1
    primes = list(factors.keys())
    # generates factors from primes[k:] subset
    def generate(k):
        if k == len(primes):
            yield 1
        else:
            rest = generate(k+1)
            prime = primes[k]
            for factor in rest:
                prime_to_i = 1
                # prime_to_i iterates prime**i values, i being all possible exponents
                for _ in range(factors[prime] + 1):
                    yield factor * prime_to_i
                    prime_to_i *= prime
    # in python3, `yield from generate(0)` would also work
    for factor in generate(0):
        yield factor

def powerset(iterable):
    "list(powerset([1,2,3])) --> [(), (1,), (2,), (3,), (1,2), (1,3), (2,3), (1,2,3)]"
    power_list = []
    s = list(iterable)
    for tup in itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1)):
        power_list.append(list(tup))
    return power_list

def rel_set_to_string(rel_set):
    string_list = []
    for i in range(len(rel_set)):
        string_ints = [str(int) for int in rel_set[i][0]]
        string_of_ints = ";".join(string_ints)
        string_list.append(string_of_ints)
    joined_string = "|".join(string_list)
    return joined_string

def get_vertex_numbering_key_from_value(vertex_numbering, vertex):
    key_list = list(vertex_numbering.keys())
    value_list = list(vertex_numbering.values())
    if vertex > 0:
        position = value_list.index(vertex)
        key = key_list[position]
    else:
        position = value_list.index(-vertex)
        key = -key_list[position]
    return key
