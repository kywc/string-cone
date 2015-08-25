import cPickle as pickle
import os.path

class ITrail:
    def __init__(self, i, start, c=[], has_end=false, has_d=false):
        if len(c) > len(i):
            raise TypeError("Sequence too long")
        self.c = c
        self.i = i
        self.start = start
        if has_end==false or has_d==false:
            self.end = self.find_end()[0]
            self.d = self.find_end()[1]
        else:
            self.end = has_end
            self.d = has_d

    def lengthen(self, c_k):
        k = len(self.c)         # this is actually k-1, but index of i list starts at 0, so would need to subtract one anyway
        gamma_k = self.end - c_k*self.i[k]
        d_k = c_k + gamma_k.scalar(self.i[k].associated_coroot())
        new_c = self.c + [c_k]
        new_d = self.d + [d_k]
        return ITrail(self.i, self.start, new_c, gamma_k, new_d)

    def find_end(self):                       #find end of trail, compute d's
        gamma_k = self.start
        new_d = []
        for k in range(len(self.c)):
            gamma_k = gamma_k - self.c[k]*self.i[k]
            new_d = new_d + [self.c[k] + gamma_k.scalar(self.i[k].associated_coroot())]
        return [gamma_k, new_d]

    def check(self, allowed_vertices):   #check if trail goes outside allowed vertices (weight space)
        trail_end = start                #or if induced map is otherwise zero (todo)
        for k in range(len(self.c)):
            try:
                allowed_vertices.index(trail_end)
            except ValueError:
                return false
            trail_end = trail_end - self.c[k]*self.i[k]
        try:
            allowed_vertices.index(trail_end)
        except ValueError:
            return false
        return true

    def visualize(self):
        m = self.start
        m2 = self.start
        lines = []
        for k in range(len(self.c)):
            m2 = m - self.c[k]*self.i[k]
            col = (self.i[k][0] + 1 / 2, self.i[k][1] + 1 / 2, self.i[k][2] + 1 / 2)
            lines.append(line([(m[0], m[1], m[2]), (m2[0], m2[1], m2[2])], corner_cutoff=0, radius = 0.014, color = "orange", arrow_head = true))
            m = m2
        return sum(lines)



def show_weights(V = WeylCharacterRing("C3")(2,1,0)):   #returns graphics object encoding weight multiplicities of representation V
    mult = V.weight_multiplicities()
    weights = mult.keys()
    if V.parent().rank() == 3:
        return sum([point3d((pt[0],pt[1],pt[2]), size=5*(mult[pt]^(1/3))) for pt in weights])
    elif V.parent().rank() < 3:
        return sum([point((pt[0],pt[1]), size=5*(mult[pt]^(1/2))) for pt in weights])
    else:
        return 0

def recompute_trails(trail, allowed_vertices):
    k = len(trail.c)
    index = trail.i[k]
    trail_list = []
    new_c = 0
    while true:
        new_trail = trail.lengthen(new_c)
        if new_trail.end not in allowed_vertices:
            break
        else:
            trail_list.append(new_trail)
            new_c = new_c + 1
    return trail_list


def compute_trails(V, indices, start, end, prune=false, show_incomplete=false):
    g = V.parent()  #parent character ring
    r = g.rank()    #rank of associated Lie algebra
    dim = r         #dimension of ambient space (same as rank except in type A, where Sage models lie alg as gl_n)
    if g.cartan_type()[0] == 'A':
        dim = r + 1
    for i in indices:
        if i < 1 or i > r:
            raise TypeError("Index greater than rank")
    alpha = g.simple_roots()

    #for i in alpha: #testing: display simple roots
    #    print i,
    #print ""

    mult = V.weight_multiplicities()
    wts = mult.keys()

    #allowed_verts = []
    #for lamb in wts:
    #    dom = true
    #    for a in alpha:
    #        if (lamb - end).scalar(a.associated_coroot()) < 0:
    #            dom = false
    #    if dom:
    #        allowed_verts = allowed_verts + [lamb]
    #print wts
    #print allowed_verts

    #wts = map(lambda x: [x[i] for i in range(dim)], wts)
    indices = map(lambda x: alpha[x], indices)             #replace list of integers with corresponding simple roots
    trails = []
    final_trails = []
    trails.append(ITrail(indices, start))

    for i in range(len(indices)):
        trails = map(lambda t: recompute_trails(t, wts), trails)    #replace each trail in list with list of valid trails going in proper direction
        #trails = map(lambda t: recompute_trails(t, allowed_verts), trails)    #replace each trail in list with list of valid trails going in proper direction
        trails = (lambda z: [x for y in z for x in y])(trails)      #flatten list

    for trail in trails:
        if trail.end == end:
            final_trails.append(trail)
    if show_incomplete:
        return [trails, final_trails]
    else:
        return final_trails

def compute_inequalities(cartan_type, word, debug=false):
    ct = CartanType(cartan_type)
    dualct = ct.dual()
    dualring = WeylCharacterRing(dualct)

    weyl = dualring.space().weyl_group()
    w0 = weyl.w0

    dualfund = dualring.fundamental_weights()

    ts = []
    ineqs = []
    for i in range(len(dualfund)):
        #wicheck = fund[i+1].associated_coroot()
        wicheck = dualfund[i+1]
        if debug:
            print wicheck
        final = wicheck.simple_reflection(i+1).weyl_action(w0)
        V = dualring(wicheck)
        ts = compute_trails(V, word, wicheck, final)
        for x in ts:
            if debug:
                print x.start, x.end, x.c, x.i
                print x.d
            ineqs = ineqs + [[0] + x.d]
        if debug:
            print ""
    return ineqs

def compute_Qw0(cartan_type, word, weight):
    ineqs = compute_inequalities(cartan_type, word)
    alpha = WeylCharacterRing(cartan_type).simple_roots()
    n = len(ineqs[0]) - 1
    #print n
    for i in range(n):
        #k = i+1
        new_ineq = [0] * (n+1)
        new_ineq[0] = weight.scalar(alpha[word[i]])
        new_ineq[i+1] = -1
        for j in range(n - i - 1):
            new_ineq[i+j+2] = -alpha[word[i+j+1]].scalar(alpha[word[i]])
            #print i+j+2
        ineqs = ineqs + [new_ineq]
    return ineqs

def compute_Cw0(cartan_type, word):
    ineqs = compute_inequalities(cartan_type, word)
    n = len(ineqs[0]) - 1
    ineqs = map(lambda ineq: ([0] * n) + ineq, ineqs)
    print ineqs
    

def compute_barCw0(cartan_type, word):
    ineqs = compute_inequalities(cartan_type, word)
    alpha = WeylCharacterRing(cartan_type).simple_roots()
    n = len(ineqs[0]) - 1
    m=len(alpha)
    for h in range(0,len(ineqs)):
        for s in range(0,m):
            ineqs[h].append(0)
    #print n
    for i in range(n):
        #k = i+1
        new_ineq = [0] * (n+m+1)
        new_ineq[0] =0
        new_ineq[i+1] = -1
        for j in range(n - i - 1):
            new_ineq[i+j+2] = -alpha[word[i+j+1]].scalar(alpha[word[i]])
            #print i+j+2
        for r in range(1,m+1):
            new_ineq[n+r]=kronecker_delta(word[i],r)
        ineqs = ineqs + [new_ineq]
    return ineqs

def Qw0_all_words(cartan_type, weight, save_to=false):
    ct = CartanType(cartan_type)
    dualct = ct.dual()
    dualring = WeylCharacterRing(dualct)
    weyl = dualring.space().weyl_group()
    w0 = weyl.w0
    R=w0.reduced_words()
    polyhedra = []
    if save_to:
        if os.path.exists(save_to):
            print "File exists"
            return []
        f = open(save_to, 'a')
    for current_word in R:
        P=Polyhedron(ieqs = compute_Qw0(cartan_type, current_word, weight))
        if save_to:
            #pickle.write(current_word + "\n")
            #pickle.dump([P, current_word], f)
            P_dump = P.dumps()
            pickle.dump([P_dump, current_word], f)
        polyhedra.append([P, current_word])
    return polyhedra

def Polytope_Partition(cartan_type, weight, from_file=false):
    polyhedra = []
    if not from_file:
        polyhedra = Qw0_all_words(cartan_type, weight)
    if from_file:
        f = open(from_file, 'r')
        while True:
            try:
                poly = pickle.load(f)
                polyhedra.append([loads(poly[0]), poly[1]])
            except EOFError:
                break
    partition=[]
    for couple in polyhedra:
        P = couple[0]
        current_word = couple[1]
        for i in range(0,len(partition)):
            if (P.face_lattice()).is_isomorphic(partition[i][0].face_lattice()):
                partition[i].append(current_word)
                break
        else:
            partition.append([])
            partition[-1].append(P)
            partition[-1].append(current_word)
    return partition

def Polytope_Partition_facelattice(cartan_type, weight, from_file=false):
    polyhedra = []
    if not from_file:
        polyhedra = Qw0_all_words(cartan_type, weight)
    if from_file:
        f = open(from_file, 'r')
        while True:
            try:
                poly = pickle.load(f)
                polyhedra.append([loads(poly[0]), poly[1]])
            except EOFError:
                break
    partition=[]
    for couple in polyhedra:
        P = couple[0]
        current_word = couple[1]
        for i in range(0,len(partition)):
            if (P.face_lattice()).is_isomorphic(partition[i][0].face_lattice()):
                partition[i].append(current_word)
                break
        else:
            partition.append([])
            partition[-1].append(P)
            partition[-1].append(current_word)
    return partition

#def Polytope_Partition(cartan_type, weight):
#    ct = CartanType(cartan_type)
#    dualct = ct.dual()
#    dualring = WeylCharacterRing(dualct)
#    weyl = dualring.space().weyl_group()
#    w0 = weyl.w0
#    R=w0.reduced_words()
#    partition=[]
#    for current_word in R:
#        P=Polyhedron(ieqs = compute_Qw0(cartan_type, current_word, weight))
#        for i in range(0,len(partition)):
#            if (P.face_lattice()).is_isomorphic(partition[i][0].face_lattice()):
#                partition[i].append(current_word)
#                break
#        else:
#            partition.append([])
#            partition[-1].append(P)
#            partition[-1].append(current_word)
#    return partition


    

def sim(A,B):
    print A==B
    print (A.face_lattice()).is_isomorphic(B.face_lattice())
    
    
