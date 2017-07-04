def find_relations(q, g, F):
    p = 2*q + 1
    if p.is_prime() == False:
        return False
    if q.is_prime() == False:
        return False
    L = GF(q)
    A = Matrix(L, len(F), len(F) + 1, range(0))
    r = A.rank()
    count = 0
    i = 1
    while count < len(F):
        x, a, b = relation_step(g, i, p, F)
        if b:
            A[count] = vector(L, x)
            #lineare Unabhängigeit testen über rang von A
            if A.rank() > r:
                r = A.rank()
                count += 1
        i += 1
    B = A[0:len(F),0:len(F)]
    v = A.column(len(F))
    w = B.solve_right(v)
    w = [ ZZ(comp) for comp in list(w) ]
    return power_step(g, p, q, F, w)

def index_calculus(q, g, F, h):
    w = find_relations(q, g, F)
    s = find_discrete_log(q, g, F, h, w)
    p = 2*q + 1
    if power_mod(g, s, p) == h:
        return s
    else:
        return "index calculus failed"

def test_run():
    return index_calculus(7901, 80, [2,3,5,7], 200)

#load('index_ZZ-cython.spyx'); load('index_ZZ-cython.sage'); test_run()
#timeit('test_run()', precision = 4, number = 2^7, repeat = 5)
