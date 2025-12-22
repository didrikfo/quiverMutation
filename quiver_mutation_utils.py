import itertools

def listIntersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def sublistExists(list, sublist):
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
    powerList = []
    s = list(iterable)
    for tup in itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1)):
        powerList.append(list(tup))
    return powerList

def relSetToString(relSet):
    stringList = []
    for i in range(len(relSet)):
        stringInts = [str(int) for int in relSet[i][0]]
        stringOfInts = ";".join(stringInts)
        stringList.append(stringOfInts)
    joinedString = "|".join(stringList)
    return joinedString

def getVertexNumberingKeyFromValue(vertexNumbering, vertex):
    keyList = list(vertexNumbering.keys())
    valueList = list(vertexNumbering.values())
    if vertex > 0:
        position = valueList.index(vertex)
        key = keyList[position]
    else:
        position = valueList.index(-vertex)
        key = -keyList[position]
    return key
