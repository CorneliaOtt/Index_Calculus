#The algorithm solves the discrete logarithm problem alpha^x=beta mod n
#The Algorithm only works in F_n with prime numbers n,q with n=2*q+1!!!
#Input: n, alpha and beta
#Output: x
#pollard_rho(n,alpha,beta) return the dicrete logarithm: log_alpha(beta)mod n =x

def new_x_a_b(x,a,b,n,alpha, beta):
    x=Integer(x)
    if x%3 ==0:
        x=(x^2) % n
        a =  (a*2)  % (n-1)
        b =  (b*2)  % (n-1)
    if x%3 ==1:
        x=x*alpha % n
        a = (a+1) % (n-1)
    if x%3 ==2:
        x=(x*beta) % n
        b = (b+1) % (n-1)
    return [x,a,b]


def pollard_rho(n,alpha,beta):
    p=Integer((n-1)/2)
    a=GF(n).random_element()
    b=GF(n).random_element()
    a=Integer(a)
    b=Integer(b)
    x=Integer(alpha^a*beta^b)
    X=x
    A=a
    B=b
    gamma=0
    for i in range (1,n-1):
        [x,a,b]=new_x_a_b( x, a, b, n,alpha, beta)
        [X,A,B]=new_x_a_b( X, A, B, n,alpha, beta)
        [X,A,B]=new_x_a_b( X, A, B, n,alpha, beta)
        if x == X :
            gamma=((a-A)*inverse_mod((B-b),p))%(p)
            if(alpha^gamma%n==beta%n):
                break
            else:
                gamma=gamma+p
                break
    return gamma

#an example
#n=19532699
#alpha=3561329
#beta=14204345
#pollard_rho(n,alpha,beta)