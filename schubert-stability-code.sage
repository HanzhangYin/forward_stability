import re
import sys
import random
from multipolynomial_bases import *
A.<x,y> = MultivariatePolynomialAlgebra(QQ)
Schub = A.schubert_basis()

def lehmer_to_permutation(code):
    out = []
    included = [False] * int(len(code) + max(code) + 1)
    for i in range(len(code)):
        num_falses_seen = 0
        j = 0
        while j < len(included):
            if not included[j]:
                num_falses_seen += 1
            if num_falses_seen == code[i] + 1:
                out.append(j + 1)
                included[j] = True
                break
            j += 1
    for i in range(1, max(out)):
        if i not in out:
            out.append(i)
    return out

def add_codes(c, d):
    out = []
    for i in range(max(len(c),len(d))):
        ci = c[i] if i<len(c) else 0
        di = d[i] if i<len(d) else 0
        out.append(ci+di)
    return out

def dual_lehmer_to_permutation(code):
    rev_code = [code[i] for i in range(len(code)-1,-1,-1)]
    n = len(lehmer_to_permutation(rev_code))
    return w0(n)*lehmer_to_permutation(rev_code)*w0(n)

#sum of nonzero rows, as defined in H-W
def ordinary_lambda(w,i,j):
    lehmer_theta = []
    for k in range(i,j+1):
        to_append = 0
        for l in range(k+1,len(w)+1):
            if w[l-1]<w[k-1]:
                to_append = 1
        lehmer_theta.append(to_append)
    return sum(lehmer_theta)

def dual_lambda(w,i,j):
    dual_lehmer_theta = []
    for k in range(i,min(len(w)+1,j+1)):
        to_append = 0
        for l in range(1,k):
            if w[l-1]>w[k-1]:
                to_append = 1
        dual_lehmer_theta.append(to_append)
    return sum(dual_lehmer_theta)

def dual_lambda_sum(w,i):
    dual_lehmer_theta = []
    for j in range(i-1,len(w)):
        to_append = 0
        for k in range(j):
            if w[k]>w[j]:
                to_append = 1
        dual_lehmer_theta.append(to_append)
    return sum(dual_lehmer_theta)

def w0(n):
    return Permutation(range(n,0,-1))

def iota(w):
    n = len(w)
    return w0(n)*w*w0(n)

def push(w, n=1):
    lc = w.to_lehmer_code()
    for i in range(n):
        lc.insert(0,0)
    return Permutation(lehmer_to_permutation(lc))

def largest_moved(w):
    to_return = len(w)
    while(to_return > 1 and w[to_return-1] == to_return):
        to_return = to_return - 1
    return to_return

def stability_number(u,v):
    starting_point = max(largest_moved(u),largest_moved(v))
    to_return = starting_point
    for i in range(starting_point,0,-1):
        term = dual_lambda_sum(u,i) + dual_lambda_sum(v,i) + i - 1
        to_return = max(to_return, term)
    return to_return


def get_support(u,v):
    uA = Schub[u.to_lehmer_code()]
    vA = Schub[v.to_lehmer_code()]
    prod_list = list(uA*vA)
    return [lehmer_to_permutation(s[0]) for s in prod_list]


def schubert_product(u,v,push_forward=0):
    print("Schubert product:")
    #print("S_%s * S_%s = " %(u,v))
    u = push(u,push_forward)
    v = push(v,push_forward)
    print("S_%s * S_%s = " %(u,v))
    prod = Schub[u.to_lehmer_code()]*Schub[v.to_lehmer_code()]
    prefix_plus = "  "
    for perm,coeff in list(prod):
        print("%s%s * S_%s" %(prefix_plus, extract_coeff(coeff), lehmer_to_permutation(perm)))
        prefix_plus = "+ "
    f_stab = max([largest_moved(lehmer_to_permutation(w[0])) for w in prod])
    print("Largest Number Moved: %s" %(f_stab))

def push(w, n=1):
    lc = w.to_lehmer_code()
    for i in range(n):
        lc.insert(0,0)
    return Permutation(lehmer_to_permutation(lc))

def extract_coeff(y):
    y = str(y)
    if y=="y[0]":
        return 1
    if y=="-y[0]":
        return -1
    else:
        return int(y[:-5])

def stability_number_from_product(u,v,cohom="Sch"):
    to_return = 0
    support = get_support(u,v) if cohom=="Sch" else getGrothPerms(u,v)
    for w in support:
        to_return = max(to_return, largest_moved(w))
    return to_return

def check_stability_number(n, cohom="Sch", verbose=False):
    for u in Permutations(n):
        for v in Permutations(n):
               formula_num = stability_number(u,v)
               actual_num = stability_number_from_product(u,v, cohom)
               if verbose or random.random() < 0.01:
                   print(u,v,formula_num,actual_num)
               if formula_num != actual_num:
                   print("ERROR!")
                   return
    print("ALL GOOD!")

def w0(n):
    return Permutation(range(n,0,-1))

def print_support(u,v):
    uA = Schub[u.to_lehmer_code()]
    vA = Schub[v.to_lehmer_code()]
    prod_list = list(uA*vA)
    print(u.to_lehmer_code(),v.to_lehmer_code(),"--",[s[0] for s in prod_list])

def check_lehmer_code_sum(n, verbose=False):
    for u in Permutations(n):
        for v in Permutations(n):
            formula_num = stability_number(u,v)
            code_sum_num = largest_moved(lehmer_to_permutation(add_codes(u.to_lehmer_code(), v.to_lehmer_code())))
            inverse_code_sum_num = largest_moved(lehmer_to_permutation(add_codes(u.inverse().to_lehmer_code(), v.inverse().to_lehmer_code())))
            dual_code_sum_num = largest_moved(dual_lehmer_to_permutation(add_codes(u.to_lehmer_cocode(), v.to_lehmer_cocode())))
            inverse_dual_code_sum_num = largest_moved(dual_lehmer_to_permutation(add_codes(u.inverse().to_lehmer_cocode(), v.inverse().to_lehmer_cocode())))
            if verbose or random.random() < 0.01:
                print(u,v,formula_num, "--", code_sum_num, inverse_code_sum_num, dual_code_sum_num, inverse_dual_code_sum_num)
            maximum = max(code_sum_num, inverse_code_sum_num, dual_code_sum_num, inverse_dual_code_sum_num)
            if formula_num != maximum:
                print("ERROR!", formula_num < maximum)
    print("ALL GOOD!")

def integer_support(w):
    to_return = set()
    for i in range(1,len(w)):
        for j in range(1,i+1):
            if w[j-1]>i:
                to_return.add(i)
    return to_return

def integer_support_from_product(u,v, cohom="Sch"):
    to_return = set()
    perm_support = get_support(u,v) if cohom=="Sch" else getGrothPerms(u,v)
    for w in perm_support:
        to_return = to_return | integer_support(w)
    return to_return

def integer_support_from_formula(u,v):
    to_return = set()
    for j in range(1,len(u)+len(v)+1):
        for i in range(1,len(u)+len(v)+1):
           if j in integer_support(u) or j in integer_support(v):
               to_return.add(j)
           if ordinary_lambda(u,j,i) + ordinary_lambda(v,j,i) > abs(j-i):
               #print("Adding %s from %s, Ordinary-u: %s,  Ordinary-v: %s, distance: %s" %(j, i, ordinary_lambda(u,j,i), ordinary_lambda(v,j,i), abs(j-i)))
               to_return.add(j)
           if dual_lambda(u,i,j) + dual_lambda(v,i,j)-1 > abs(j-i):
               #print("Adding %s from %s, Dual-u: %s,  Dual-v: %s, distance: %s" %(j, i, dual_lambda(u,j,i), dual_lambda(v,j,i), abs(j-i)))
               to_return.add(j)
    return to_return


def check_integer_support(n, cohom="Sch", verbose=False, quit_on_error=True):
    for u in Permutations(n):
        for v in Permutations(n):
               formula_support = integer_support_from_formula(u,v)
               actual_support = integer_support_from_product(u,v, cohom)
               if verbose or random.random() < 0.01:
                   print(u,v,formula_support,actual_support)
               if formula_support != actual_support:
                   print("ERROR!")
                   if quit_on_error:
                       return
    print("ALL GOOD!")

def random_integer_support_check(n, k, push_num=0, cohom="Sch", verbose=False, quit_on_error=True):
    perms = Permutations(n)
    for i in range(k):
        u = push(random.choice(perms),push_num)
        v = push(random.choice(perms),push_num)
        formula_support = integer_support_from_formula(u,v)
        actual_support = integer_support_from_product(u,v, cohom)
        if verbose or random.random() < 0.01:
            print(i+1,u,v,formula_support,actual_support)
        if formula_support != actual_support:
            print("ERROR!")
            if quit_on_error:
                return
    print("ALL GOOD!")
