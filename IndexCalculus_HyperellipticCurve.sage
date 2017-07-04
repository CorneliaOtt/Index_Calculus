def build_factorbase( C, q):
    #only put one of the opposite divisors into the factorbasis
    #so for finite fields Fq we test if y is in [0,(q+1)/2]
    bound=Integer((q+1)/2)
    divlist=[]
    J=C.jacobian()
    for i in range(len(C.points())):
        P=C.points()[i]
        #if P is a finite point
        if(P[2]!=0):
            y_i=P[1]
            if(y_i in range(bound)):
                D=J(P)
                divlist.append([D,P])
    #returns a list of pairs[D,P] with a Divisor D and the according point P
    return divlist


def init_random_walk( C, q, r, D_1, D_2, n):
    T=[]
    #j in [0, r-1]
    for j in range(r):
        #alpha and beta are random integers in [0,n-1]
        alpha=GF(n).random_element()
        beta=GF(n).random_element()
        alpha=Integer(alpha)
        beta=Integer(beta)
        T.append([alpha*D_1+beta*D_2, alpha, beta])
    return T


def is_smooth_divisor(D):
    a=D[0]
    b=D[1]
    d1=count_roots_with_multiplicity(a)
    d2=a.degree()
    #if d1=d2 a splits completly over F_q
    if(d1==d2):
        return true
    else:
        return false
    
    
def count_roots_with_multiplicity(a):
    list=a.roots()
    n=len(list)
    sum =0
    for i in range(n):
        sum = sum+list[i][1]
    return sum


def hash_divisor(D,r):
    #hash function h returns values in [0,r-1]
    a=D[0]
    b=D[1]
    if(b!=0):
        h= hash(tuple([a,b]))
    else:
         h= hash(a)
    h=h%r
    return h


def main_loop(D1, D2, r, T, G, n ):
    alpha_0=GF(n).random_element()
    beta_0=GF(n).random_element()
    alpha=[]
    beta=[]
    M=matrix(GF(n), len(G), len(G),range(0))
    R_0=Integer(alpha_0)*D_1+Integer(beta_0)*D_2
    k=0
    count =0
    while k <len(G):
        #step a: find a smooth divisor R_0
        j=hash_divisor(R_0,r)
        alpha_0=(alpha_0+T[j][1]) %n
        beta_0=(beta_0+T[j][2])%n
        R_0=Integer(alpha_0)*D1+Integer(beta_0)*D2
        while(is_smooth_divisor(R_0)==false):
            j=hash_divisor(R_0,r)
            alpha_0=alpha_0+T[j][1] %n
            beta_0=beta_0+T[j][2]%n
            R_0=Integer(alpha_0)*D1+Integer(beta_0)*D2
        #step b: express R_0 on the factorbasis G
        a=R_0[0]
        if(a!=1):
            roots=a.roots()
            m=len(roots)
            #mk is the k-th row in the matrix M
            mk=[]
            Rk=0
            for j in range(len(G)):
                root_found=false
                gj=G[j][0][0]
                gj_root=G[j][1][0]
                for i in range(m):
                    if(roots[i][0]== gj_root):
                        mk.append(roots[i][1])
                        root_found=true
                        D=roots[i][1]*G[j][0]
                        Rk=Rk+D
                        break
                if(root_found==false):
                    mk.append(0)
            # print 'Rk='
#             print Rk
#             print "R_0 ="
#             print R_0
#             print '='
#             print R_0[0].factor()
            if((Rk[0]==R_0[0]) and (Rk[1]==R_0[1])):
                # print 'mk'
#                 print mk
                M[k]=mk
                alpha.append(alpha_0)
                beta.append(beta_0)
                k=k+1
        count=count+1
        if count>(len(G)*10000):
            return [-1,-1,-1]
    return [M, alpha, beta]



def find_kernel(M,n,q,G):
    #find a non zero vector of the kernel of the matrix
    Ker=kernel(M)
    Ker=Ker.basis_matrix()
    #v zero-vector
    v=vector(GF(n), len(G))
    for i in range(Ker.dimensions()[0]):
        if(Ker[i]!=v):
            gamma=Ker[i]
            break
    return gamma



def solution(alpha, beta, gamma, n, G):
    counter=0
    denominator=0
    for k in range(len(G)):
        counter = (counter+alpha[k]*gamma[k])%n
        denominator=(denominator+beta[k]*gamma[k])%n
    if (denominator==0):
        return -1
    else:
        factor=inverse_mod(Integer(denominator),n)
        x=-counter*factor
        x=x%n
        return x
    
    
def index_calculus(C,q,r,D1,D2,n):
    G=build_factorbase( C, q)
    # print G
#     print ""
    count=0
    x=0
    while(count<20):
        print count
        count=count+1
        T=init_random_walk( C, q, r, D1, D2, n)
        [M,alpha,beta]=main_loop(D1, D2, r, T, G, n)
        # print [M,alpha,beta]
#         print ""
        if(M!=-1):
            gamma=find_kernel(M,n,q,G)
            x=solution(alpha,beta,gamma, n, G)
            if(x!=-1):
                break
    return x

#an example
q=53
n=4451
r=log(n,2).round()
K.<x>=GF(q)[]
f=40*x^7 + 31*x^6 + 9*x^5 + 22*x^4 + 34*x^3 + 47*x^2 + 7*x + 48
h=0
C=HyperellipticCurve(f,h)
J=C.jacobian()
P=C([-6,-43,1])
D_1=J(P)
D_1
D_2=(n-5)*D_1
D_2


%time
index_calculus(C,q,r,D_1,D_2,n)