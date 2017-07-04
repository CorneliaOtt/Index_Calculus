def find_relations(q, g, F):
    p = 2*q + 1
    if p.is_prime() == False:
        return False
    if q.is_prime() == False:
        return False
    K = GF(p)
    L = GF(q)
    #Matrix in die die Relationen zeilenweise gespeichert werden
    A = Matrix(L, len(F), len(F) + 1, range(0))
    r = A.rank()
    count = 0
    i = 1
    while count < (len(F)):
        temp = power_mod(g, i, p)
        x = [ temp.valuation(pr) for pr in F ]
        x.append(i)
        y = vector(x)
        a = [ pf[0] for pf in factor(temp) ]
        b = True
        #y wird nur aufgenommen falls die Faktoren ausschließlich aus der Faktorbasis stammen
        for h in a:
            if h not in F:
                b = False
        if b:
            A[count] = y
            #lineare Unabhängigeit testen über rang von A
            if A.rank() > r:
                r = A.rank()
                count += 1
        i += 1
    B = A[0:len(F),0:len(F)]
    v = A.column(len(F))
    w = B.solve_right(v)
    w = vector(IntegerModRing(2*q),w)
    a = [ K(g)^(ZZ(comp)) for comp in w ]
    b = [ K(g)^(ZZ(comp) + q) for comp in w ]
    for k in range(len(F)):
        if a[k] > b[k]:
            w[k] += q
    return w

def find_discrete_log(q, g, F, h, w):
    p = 2*q + 1
    i = 1
    b = True
    while b:
        temp = (power_mod(g, i, p) * h) % p
        x = [ temp.valuation(pr) for pr in F ]
        x.append(i)
        y = vector(x)
        a = [ pf[0] for pf in factor(temp) ]
        b = False
        for j in a:
            if j not in F:
                b = True
        i += 1
        s = 0
    for j in range(len(F)):
        s += y[j]*w[j]
    s = (s - y[len(F)]) % (p - 1)
    return ZZ(s)

def index_calculus(q, g, F, h):
    w = find_relations(q, g, F)
    #w = [ ZZ(comp) for comp in w ]
    s = find_discrete_log(q, g, F, h, w)
    p = 2*q + 1
    if power_mod(g, s, p) == h:
        return s
    else:
        return "index calculus failed"

def test_run():
    return index_calculus(7901, 80, [2,3,5,7], 200)

#load('index_ZZ.sage'); test_run()
#timeit('test_run()', precision = 4, number = 2^7, repeat = 5)
