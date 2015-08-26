def add(a,b):
    return int(a+b)

def multiply(a,b):
    return int(a*b)

def divide(a,b):
    return int(a/b)

def find(s, ch):
    return [(i+1,) for i, ltr in enumerate(s) if ltr == ch]


def find_subword(word,subword):

    arr=[]
    for i in subword:
        arr.append(find(word,i))

    current=arr[0][:]
    temp=[]
    for k in range(1,len(arr)):
        for j in arr[k]:
            for i in current:
                if(i[-1]<j[-1]):
                    temp.append(i+j)
        current=temp[:]
        temp=[]
    return current


def find_subwords(word,subwords):
    arr=[]
    for i in subwords:
        arr.extend(find_subword(word,i))
    return arr


def k(j,t,n):
    if j==0:
        return 0
    elif j>len(t):
        return n*(n-1)/2+1
    else:
        return t[j-1]


def omega(i,n):
    om=[]
    for k in range(0,i):
        om.append(1)
    for k in range(i,n):
        om.append(0)
    return om


def alpha(i,n):
    al=[]
    for k in range(0,i-1):
        al.append(0)
    al.append(1)
    al.append(-1)
    for k in range(i+1,n):
        al.append(0)
    return al


def innerprod(a,b,n):
    inner=0
    for i in range(0,n):
        inner=inner+a[i]*b[i]
    return inner


def permute(h,t,n):
    if h.order()>1:
        perm=[]
        for i in range(0,n):
            perm.append(t[h.dict()[i+1]-1])
        return perm
    else:
        return t


def addlambdainequality(lamb,Z,l,word,n):#should now be correct
    P=PolynomialRing(ZZ,l+1,"t")
    t=P.gens()
    sum=0
    for k in range(0,l):
        for i in range(k+1,l):
            sum=sum-innerprod(alpha(word[i],n),alpha(word[k],n),n)*t[i+1]
        Z.append(innerprod(lamb,alpha(word[k],n),n)*t[0]+sum-t[k+1])
        sum=0
    return Z

def lambdainequality(lamb,l):
    P=PolynomialRing(ZZ,l+1,"t")
    t=P.gens()
    Z=[]     
    sum=0
    for k in range(0,l):
        for i in range(k+1,l):
            sum=sum-innerprod(alpha(word[i],n),alpha(word[k],n),n)*t[i+1]
        Z.append(innerprod(lamb,alpha(word[k],n),n)*t[0]+sum-t[k+1])
        sum=0
    return Z

def addfulllambdainequality(Z,l,word,n):
    P=PolynomialRing(ZZ,l+n+1,"t")
    t=P.gens()
    sum=0
    for k in range(0,l):
        for i in range(k+1,l):
            sum=sum-innerprod(alpha(word[i],n),alpha(word[k],n),n)*t[i+1]
        Z.append(t[l+word[k]]+sum-t[k+1])
        sum=0
    return Z

def Matrixlist(ineq,l):
    Mlist=[]
    for i in range(0,len(ineq)):
        Mlist.append([])
        for k in range(0,l+1):
            ev=[]
            for s in range(0,l+1):
                ev.append(kronecker_delta(s,k))
            Mlist[i].append((ineq[i])(ev))
    return Mlist

def Totalconematrixlist(ineq,l,n):
    Mlist=[]
    for i in range(0,len(ineq)):
        Mlist.append([])
        for k in range(0,add(l,n)+1):
            ev=[]
            for s in range(0,add(l,n)+1):
                ev.append(kronecker_delta(s,k))
            Mlist[i].append((ineq[i])(ev))
    return Mlist

def macaulay2matrixlist(ineq,l):
    Mlist=[]
    for i in range(0,len(ineq)):
        Mlist.append([])
        for k in range(1,l+1):
            ev=[]
            for s in range(0,l+1):
                ev.append(-kronecker_delta(s,k))
            Mlist[i].append((ineq[i])(ev))
    Mlist1=str(Mlist)
    Mlist2=Mlist1.replace('[','{')
    Mlist3=Mlist2.replace(']','}')
    return Mlist3


def Vectorlist(ineq,l):
    Vlist=[]
    for i in range(0,len(ineq)):
        k=0
        ev=[]
        for s in range(0,l+1):
            ev.append(kronecker_delta(s,k))
        Vlist.append([(ineq[i])(ev)])
    return Vlist

def macaulay2vectorlist(ineq,l):
    Vlist=[]
    for i in range(0,len(ineq)):
        k=0
        ev=[]
        for s in range(0,l+1):
            ev.append(kronecker_delta(s,k))
        Vlist.append([(ineq[i])(ev)])
    Vlist1=str(Vlist)
    Vlist2=Vlist1.replace('[','{')
    Vlist3=Vlist2.replace(']','}')
    return Vlist3

def findpolytope(n,lam,number):
    G=SymmetricGroup(n)    #S3
    w0=G.w0  #as a group element
    w0_words = G.w0.reduced_words()
    gens=[]
    W=[]
    for i in range(1, n):
        gens.append(G.simple_reflection(i))
    #Finding the Wi
    W=[]
    for i in gens:
        b=gens[:]
        b.remove(i)
        W.append(G.subgroup(b))
    #Finding the Ui
    U=[]
    i=1
    for j in W:
        elmt = gens[i-1]*w0 #need new name, i-1 because list index, i corresponds to s_1
        rightCosets=G.cosets(j,side="right")
        for sublist in rightCosets:
            if elmt in sublist:
                coset=sublist
                break
        u=min(coset)
        U.append(u.reduced_words()) #we don't want the group element, we want it expressed as s_i
        i+=1
    #print "U:", U
    #print "\n\n"
    word=list(w0_words[number])
    #word=list(w0_words[12]) #for pablo's choice in s4
    l=len(word)
    print "Longest Word:", word
    print "\n\n"
    #finding U(i) as a subword of w0
    subwords=[]
    for u in U:
        subwords.append(find_subwords(word,u)) #get permutation group element as a list, list of list of  tuples of indices of u in word
    #print "subwords are \n\n", subwords
    #for each i
    P=PolynomialRing(ZZ,l+1,"t")
    t=P.gens()
    inequalities=[]
    for i in range(1,n):
        #print "i=",i,"\n"
        for current_subword in subwords[i-1]:
            #print "Current subword="
            #current_subword
            sum=0
            #the i_j are precisely the indices in the tuple current_subword... i.e. U(1)=S_1=(1,) means S_i_1 = 1.  If U(1)=S_1=(3,),     then S_i_1=(3).
            for j in range(0,len(current_subword)+1):
                perm_elmt=PermutationGroupElement(()) #product of s_i_j's
                for s in range(1,j+1):
                    perm_elmt=perm_elmt*G.simple_reflection(word[k(j,current_subword,n)-1])
                #print "J=",j,"\n" 
                for kk in range(n)+1,k(j+1,current_subword,n)):
                    #print "K=",kk,"\n"
                    #print "hi",perm_elmt
                    sum=sum+innerprod(permute(perm_elmt,alpha(word[kk-1],n),n),omega(i,n),n)*t[kk]
            inequalities.append(sum)    
    inequalities
    ineq=addlambdainequality(lam,inequalities,l,word,n) #t0 will act as a placeholder for the constant term making the next step nice
    ineq
    return ineq


def findstringcone(n,number):
    G=SymmetricGroup(n)    #S3
    w0=G.w0  #as a group element
    w0_words = G.w0.reduced_words()
    gens=[]
    W=[]
    for i in range(1, n):
        gens.append(G.simple_reflection(i))
    #Finding the Wi
    W=[]
    for i in gens:
        b=gens[:]
        b.remove(i)
        W.append(G.subgroup(b))
    #Finding the Ui
    U=[]
    i=1
    for j in W:
        elmt = gens[i-1]*w0 #need new name, i-1 because list index, i corresponds to s_1
        rightCosets=G.cosets(j,side="right")
        for sublist in rightCosets:
            if elmt in sublist:
                coset=sublist
                break
        u=min(coset)
        U.append(u.reduced_words()) #we don't want the group element, we want it expressed as s_i
        i+=1
    #print "U:", U
    #print "\n\n"
    word=list(w0_words[number])
    #word=list(w0_words[12]) #for pablo's choice in s4
    l=len(word)
    #print "Longest Word:", word
    #print "\n\n"
    #finding U(i) as a subword of w0
    subwords=[]
    for u in U:
        subwords.append(find_subwords(word,u)) #get permutation group element as a list, list of list of  tuples of indices of u in word
    #print "subwords are \n\n", subwords
    #for each i
    P=PolynomialRing(ZZ,l+1,"t")
    t=P.gens()
    inequalities=[]
    for i in range(1,n):
        #print "i=",i,"\n"
        for current_subword in subwords[i-1]:
            #print "Current subword="
            #current_subword
            sum=0
            #the i_j are precisely the indices in the tuple current_subword... i.e. U(1)=S_1=(1,) means S_i_1 = 1.  If U(1)=S_1=(3,),     then S_i_1=(3).
            for j in range(0,len(current_subword)+1):
                perm_elmt=PermutationGroupElement(())
                for s in range(1,j+1): #this line is wrong
                    perm_elmt=perm_elmt*G.simple_reflection(word[k(s,current_subword,n)-1])
                for kk in range(k(j,current_subword,n)+1,k(j+1,current_subword,n)):
                    #print "K=",kk,"\n"
                    #print "hi",perm_elmt
                    #print current_subword
                    sum=sum+innerprod(permute(perm_elmt,alpha(word[kk-1],n),n),omega(i,n),n)*t[kk]
            inequalities.append(sum)    
    return inequalities

def findstringconemorevariables(n,number):
    G=SymmetricGroup(n)    #S3
    w0=G.w0  #as a group element
    w0_words = G.w0.reduced_words()
    gens=[]
    W=[]
    for i in range(1, n):
        gens.append(G.simple_reflection(i))
    #Finding the Wi
    W=[]
    for i in gens:
        b=gens[:]
        b.remove(i)
        W.append(G.subgroup(b))
    #Finding the Ui
    U=[]
    i=1
    for j in W:
        elmt = gens[i-1]*w0 #need new name, i-1 because list index, i corresponds to s_1
        rightCosets=G.cosets(j,side="right")
        for sublist in rightCosets:
            if elmt in sublist:
                coset=sublist
                break
        u=min(coset)
        U.append(u.reduced_words()) #we don't want the group element, we want it expressed as s_i
        i+=1
    #print "U:", U
    #print "\n\n"
    word=list(w0_words[number])
    #word=list(w0_words[12]) #for pablo's choice in s4
    l=len(word)
    #print "Longest Word:", word
    #print "\n\n"
    #finding U(i) as a subword of w0
    subwords=[]
    for u in U:
        subwords.append(find_subwords(word,u)) #get permutation group element as a list, list of list of  tuples of indices of u in word
    #print "subwords are \n\n", subwords
    #for each i
    P=PolynomialRing(ZZ,l+n+1,"t")
    t=P.gens()
    inequalities=[]
    for i in range(1,n):
        #print "i=",i,"\n"
        for current_subword in subwords[i-1]:
            #print "Current subword="
            #current_subword
            sum=0
            #the i_j are precisely the indices in the tuple current_subword... i.e. U(1)=S_1=(1,) means S_i_1 = 1.  If U(1)=S_1=(3,),     then S_i_1=(3).
            for j in range(0,len(current_subword)+1):
                perm_elmt=PermutationGroupElement(()) #product of s_i_j's
                for s in range(1,j+1):
                    perm_elmt=perm_elmt*G.simple_reflection(word[k(j,current_subword,n)-1])
                #print "J=",j,"\n"
                for kk in range(n)+1,k(j+1,current_subword,n)):
                    #print "K=",kk,"\n"
                    #print current_subword
                    sum=sum+innerprod(permute(perm_elmt,alpha(word[kk-1],n),n),omega(i,n),n)*t[kk]
            inequalities.append(sum)    
    return inequalities


def findstringconematrix(n,number): 
    G=SymmetricGroup(n)    #S3
    w0=G.w0  #as a group element
    w0_words = G.w0.reduced_words()
    gens=[]
    W=[]
    for i in range(1, n):
        gens.append(G.simple_reflection(i))
    #Finding the Wi
    W=[]
    for i in gens:
        b=gens[:]
        b.remove(i)
        W.append(G.subgroup(b))
    #Finding the Ui
    U=[]
    i=1
    for j in W:
        elmt = gens[i-1]*w0 #need new name, i-1 because list index, i corresponds to s_1
        rightCosets=G.cosets(j,side="right")
        for sublist in rightCosets:
            if elmt in sublist:
                coset=sublist
                break
        u=min(coset)
        U.append(u.reduced_words()) #we don't want the group element, we want it expressed as s_i
        i+=1
    #print "U:", U
    #print "\n\n"
    word=list(w0_words[number])
    #word=list(w0_words[12]) #for pablo's choice in s4
    l=len(word)
    print "Longest Word:", word
    print "\n\n"
    #finding U(i) as a subword of w0
    subwords=[]
    for u in U:
        subwords.append(find_subwords(word,u)) #get permutation group element as a list, list of list of  tuples of indices of u in word
    #print "subwords are \n\n", subwords
    #for each i
    P=PolynomialRing(ZZ,l+1,"t")
    t=P.gens()
    inequalities=[]
    for i in range(1,n):
        #print "i=",i,"\n"
        for current_subword in subwords[i-1]:
            #print "Current subword="
            #current_subword
            sum=0
            #the i_j are precisely the indices in the tuple current_subword... i.e. U(1)=S_1=(1,) means S_i_1 = 1.  If U(1)=S_1=(3,),     then S_i_1=(3).
            for j in range(0,len(current_subword)+1):
                perm_elmt=PermutationGroupElement(()) #product of s_i_j's
                for s in range(1,j+1):
                    perm_elmt=perm_elmt*G.simple_reflection(word[k(j,current_subword,n)-1])
                #print "J=",j,"\n"
                for kk in range(n)+1,k(j+1,current_subword,n)):
                    #print "K=",kk,"\n"
                    #print "hi",perm_elmt
                    sum=sum+innerprod(permute(perm_elmt,alpha(word[kk-1],n),n),omega(i,n),n)*t[kk]
            inequalities.append(sum)    
    A=matrix(Matrixlist(inequalities,n*(n-1)/2))
    return A
    

def findstringconeonly(n,number):
    G=SymmetricGroup(n)    #S3
    w0=G.w0  #as a group element  
    w0_words = G.w0.reduced_words()
    gens=[]
    W=[]
    for i in range(1, n):
        gens.append(G.simple_reflection(i))
    #Finding the Wi
    W=[]
    for i in gens:
        b=gens[:]
        b.remove(i)
        W.append(G.subgroup(b))
    #Finding the Ui
    U=[]
    i=1
    for j in W:
        elmt = gens[i-1]*w0 #need new name, i-1 because list index, i corresponds to s_1
        rightCosets=G.cosets(j,side="right")
        for sublist in rightCosets:
            if elmt in sublist:
                coset=sublist
                break
        u=min(coset)
        U.append(u.reduced_words()) #we don't want the group element, we want it expressed as s_i
        i+=1
    #print "U:", U
    #print "\n\n"
    word=list(w0_words[number])
    #word=list(w0_words[12]) #for pablo's choice in s4
    l=len(word)
    #print "Longest Word:", word
    #print "\n\n"
    #finding U(i) as a subword of w0
    subwords=[]
    for u in U:
        subwords.append(find_subwords(word,u)) #get permutation group element as a list, list of list of  tuples of indices of u in word
    #print "subwords are \n\n", subwords
    #for each i
    P=PolynomialRing(ZZ,l+1,"t")
    t=P.gens()
    inequalities=[]
    for i in range(1,n):
        #print "i=",i,"\n"
        for current_subword in subwords[i-1]:
            #print "Current subword="
            #current_subword
            sum=0
            #the i_j are precisely the indices in the tuple current_subword... i.e. U(1)=S_1=(1,) means S_i_1 = 1.  If U(1)=S_1=(3,),     then S_i_1=(3).
            for j in range(0,len(current_subword)+1):
                perm_elmt=PermutationGroupElement(()) #product of s_i_j's
                for s in range(1,j+1):
                    perm_elmt=perm_elmt*G.simple_reflection(word[k(j,current_subword,n)-1])
                #print "J=",j,"\n"
                for kk in range(n)+1,k(j+1,current_subword,n)):
                    #print "K=",kk,"\n"
                    #print "hi",perm_elmt
                    sum=sum+innerprod(permute(perm_elmt,alpha(word[kk-1],n),n),omega(i,n),n)*t[kk]
            inequalities.append(sum)    
    return inequalities

def partition(n):
    G=SymmetricGroup(n)
    part=[]
    for i in range(0,len(G.w0.reduced_words())):
        S=findstringconeonly(n,i)
        for ss in range(0,len(part)):
            if set(S)==set(part[ss][0]):
                part[ss].append(S)
                break
        else:
            part.append([])
            part[len(part)-1].append(S)
    for i in range(0,len(part)):
        print part[i]
        print len(part[i])
        
def numericalpartition(n):
    G=SymmetricGroup(n)
    part=[]
    for i in range(0,len(G.w0.reduced_words())):
        S=[findstringconeonly(n,i),i]
        for ss in range(0,len(part)):
            if set(S[0])==set(part[ss][0][0]):
                part[ss].append(S)
                break
        else:
            part.append([])
            part[len(part)-1].append(S)
    for i in range(0,len(part)):
        group=[]
        for ss in range(0,len(part[i])):
            group.append(part[i][ss][1])
        print group
        print len(part[i])
        

def numericalpartitionlattice(n,lamb):
    G=SymmetricGroup(n)
    part=[]
    for i in range(0,len(G.w0.reduced_words())):
        S=[fullfindpolytope(n,lamb,i),i]
        for ss in range(0,len(part)):
            if (S[0].face_lattice()).is_isomorphic(part[ss][0][0].face_lattice()):
                part[ss].append(S)
                break
        else:
            part.append([])
            part[len(part)-1].append(S)
    for i in range(0,len(part)):
        group=[]
        for ss in range(0,len(part[i])):
            group.append(part[i][ss][1])
        print group 
        print len(part[i])
        
        

def fullfindpolytope(n,lam,number):
    G=SymmetricGroup(n)    #S3
    w0=G.w0  #as a group element
    w0_words = G.w0.reduced_words()
    gens=[]
    W=[]
    for i in range(1, n):
        gens.append(G.simple_reflection(i))
    #Finding the Wi
    W=[]
    for i in gens:
        b=gens[:]
        b.remove(i)
        W.append(G.subgroup(b))
    #Finding the Ui
    U=[]
    i=1
    for j in W:
        elmt = gens[i-1]*w0 #need new name, i-1 because list index, i corresponds to s_1
        rightCosets=G.cosets(j,side="right")
        for sublist in rightCosets:
            if elmt in sublist:
                coset=sublist
                break
        u=min(coset)
        U.append(u.reduced_words()) #we don't want the group element, we want it expressed as s_i
        i+=1
    #print "U:", U
    #print "\n\n"
    word=list(w0_words[number])
    #word=list(w0_words[12]) #for pablo's choice in s4
    l=len(word)
    #print "Longest Word:", word
    #print "\n\n"
    #finding U(i) as a subword of w0
    subwords=[]
    for u in U:
        subwords.append(find_subwords(word,u)) #get permutation group element as a list, list of list of  tuples of indices of u in word
    #print "subwords are \n\n", subwords
    #for each i
    P=PolynomialRing(ZZ,l+1,"t")
    t=P.gens()
    inequalities=[]
    for i in range(1,n):
        #print "i=",i,"\n"
        for current_subword in subwords[i-1]:
            #print "Current subword="
            #current_subword
            sum=0
            #the i_j are precisely the indices in the tuple current_subword... i.e. U(1)=S_1=(1,) means S_i_1 = 1.  If U(1)=S_1=(3,),     then S_i_1=(3).
            for j in range(0,len(current_subword)+1):
                perm_elmt=PermutationGroupElement(()) #product of s_i_j's
                for s in range(1,j+1):
                    perm_elmt=perm_elmt*G.simple_reflection(word[k(j,current_subword,n)-1])
                #print "J=",j,"\n"
                for kk in range(n)+1,k(j+1,current_subword,n)):
                    #print "K=",kk,"\n"
                    #print "hi",perm_elmt 
                    sum=sum+innerprod(permute(perm_elmt,alpha(word[kk-1],n),n),omega(i,n),n)*t[kk]
            inequalities.append(sum)    
    inequalities
    ineq=addlambdainequality(lam,inequalities,l,word,n) #t0 will act as a placeholder for the constant term making the next step nice
    ineq
    A=Matrixlist(ineq,n*(n-1)/2)
    B=Polyhedron(ieqs=A)
    return B

def findtotalcone(n,i):
    C=findstringcone(n,i)
    G=SymmetricGroup(n)
    R=G.w0.reduced_words()
    l=(n*(n-1))/2
    l=int(l)
    addfulllambdainequality(C,l,R[i],n)
    return C

def wrongfullfindpolytope(n,lam,number):
    G=SymmetricGroup(n)    #S3
    w0=G.w0  #as a group element
    w0_words = G.w0.reduced_words()
    gens=[]
    W=[]
    for i in range(1, n):
        gens.append(G.simple_reflection(i))
    #Finding the Wi
    W=[]
    for i in gens:
        b=gens[:]
        b.remove(i)
        W.append(G.subgroup(b))
    #Finding the Ui
    U=[]
    i=1
    for j in W:
        elmt = gens[i-1]*w0 #need new name, i-1 because list index, i corresponds to s_1
        rightCosets=G.cosets(j,side="right")
        for sublist in rightCosets:
            if elmt in sublist:
                coset=sublist
                break
        u=min(coset)
        U.append(u.reduced_words()) #we don't want the group element, we want it expressed as s_i
        i+=1
    #print "U:", U
    #print "\n\n"
    word=list(w0_words[number])
    #word=list(w0_words[12]) #for pablo's choice in s4
    l=len(word)
    #print "Longest Word:", word
    #print "\n\n"
    #finding U(i) as a subword of w0
    subwords=[]
    for u in U:
        subwords.append(find_subwords(word,u)) #get permutation group element as a list, list of list of  tuples of indices of u in word
    #print "subwords are \n\n", subwords
    #for each i
    P=PolynomialRing(ZZ,l+1,"t")
    t=P.gens()
    inequalities=[]
    for i in range(1,n):
        #print "i=",i,"\n"
        for current_subword in subwords[i-1]:
            #print "Current subword="
            #current_subword
            sum=0
            #the i_j are precisely the indices in the tuple current_subword... i.e. U(1)=S_1=(1,) means S_i_1 = 1.  If U(1)=S_1=(3,),     then S_i_1=(3).
            for j in range(0,len(current_subword)+1):
                perm_elmt=PermutationGroupElement(()) #product of s_i_j's
                for s in range(1,j+1):
                    perm_elmt=perm_elmt*G.simple_reflection(word[k(n,current_subword,n)-1])
                #print "J=",j,"\n"
                for kk in range(k(j,current_subword,n)+1,k(j+1,current_subword,n)):
                    #print "K=",kk,"\n"
                    #print "hi",perm_elmt 
                    sum=sum+innerprod(permute(perm_elmt,alpha(word[kk-1],n),n),omega(i,n),n)*t[kk]
            inequalities.append(sum)    
    inequalities
    ineq=addlambdainequality(lam,inequalities,l,word,n) #t0 will act as a placeholder for the constant term making the next step nice
    ineq
    A=Matrixlist(ineq,n*(n-1)/2)
    B=Polyhedron(ieqs=A)
    return B

