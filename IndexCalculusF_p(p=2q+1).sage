#The algorithm solves the discrete logarithm problem g^x=h mod p
#The Algorithm only works in F_p with prime numbers p,q with p=2*q+1!!!
#Input: q, g, F (the factorbase) and h
#Output: x
#index_calculus(q,g,F,h) return the dicrete logarithm: log_g(h)mod p =x

def find_relations( q, g, F):
    p=2*q+1
    if p.is_prime()==false:
        print p
        print "p is not prime"
        return false
    if q.is_prime()==false:
        print q
        print "q is not prime"
        return false
    K=GF(p)
    L=GF(q)
    #matrix A to save the relations row-wise
    A=matrix(L, len(F), len(F)+1, range(0))
    r=A.rank()
    count = 0
    i=1
    while count < (len(F)) :
        temp=g^i %p
        x=[temp.valuation(pr) for pr in F]
        x.append(i)
        y=vector(x)
        a= [pf[0] for pf in factor(temp)]
        b= true
        #save y if and only if all factors are in the facorbase F
        for h in a:
            if h not in F:
                b=false
        if b==true:
            y
            A[count]=y
            # test the linear independence with the rank of A
            if A.rank()> r:
                r=A.rank()
                count=count+1
        i=i+1
    B=A[0:len(F),0:len(F)]
    v=A.column(len(F))
    w = B.solve_right(v)
    w=vector(IntegerModRing(2*q),w)
    a=[ K(g)^(ZZ(bla)) for bla in w ]
    b=[ K(g)^(ZZ(bla) + q) for bla in w ]
    for k in range(len(F)):
        if a[k]>b[k]:
            w[k]=w[k]+q
    return w



def find_discrete_log(q,g,F,h,w):
    p=2*q+1
    K=GF(p)
    L=GF(q)
    i=1
    b=true
    while b :
        temp=(g^i * h) %p
        x=[temp.valuation(pr) for pr in F]
        x.append(i)
        y=vector(x)
        a= [pf[0] for pf in factor(temp)]
        b= false
        for j in a:
            if j not in F:
                b=true
        i=i+1
        s=0
    for j in range(len(F)):
        s=s+y[j]*w[j]
    s=(s-y[len(F)])%(p-1)
    return s



def index_calculus(q,g,F,h):
    p=2*q+1
    w = find_relations(q,g,F)
    print w
    s=find_discrete_log(q,g,F,h,w)
    if (g^s)%p ==h:
        return s
    else:
        return "index calculus failed"



#an example
#q=576500911530019

#K = FiniteField(p)
#g=4949832
#h=3708399
#index_calculus(q, ZZ(g), [2,3,5,7,11,13,17,19,23], ZZ(h))