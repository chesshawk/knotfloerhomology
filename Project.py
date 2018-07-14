import itertools
sequence = []
#import a signed sequence
# get the max of the sequence in absolute value
def get_max(list):
    maximum = 0
    for i in list:
        if abs(i) > maximum:
            maximum = abs(i)
    return maximum
#/ over \ is a + in the sequence, \ over / is a - in the sequence

sequence = [1,4,1] #This is where you tell the program what crossings are present in the knot's tangle decomposition
#If any element in the above array is zero I will find you and I will kill you. 

#supposing notation is that i corresponds to a short crossing between i and i + 1

#if max = i then number of boundary points is 2(i+1)
boundary_points = 2*(get_max(sequence)+1)
#print(boundary_points)

tangles = len(sequence) + boundary_points
print(tangles)

#further tangle decomposition

crossing_tangles = len(sequence) #number of tangles with a cross

caps = get_max(sequence) + 1

coordinates = []
def make_grid(list): 
    for i in range(tangles-1):
        for j in range(boundary_points):
            list.append((i,j))
make_grid(coordinates)
# number of coordinates is tangles times boundaries
_matrix = [0] * len(coordinates)
for i in range(len(coordinates)):
    _matrix[i] = [0] * len(coordinates)
 
def get_entry(x, y): #gets entry in the list of coordinates for the knot
    if 0 <= x < tangles-1 and 0 <= y < boundary_points:
        x = boundary_points * x
        return x + y
    else:
        print("Out of Bounds") #TODO come back and figure this out with exceptions

def get_back(z): #inverse of get_entry
    if z > boundary_points * (tangles - 1):
        print("Too high, bruh")
    x = z//boundary_points  #TODO come back and figure this out with exceptions
    y = z%boundary_points
    return (x,y)
def connected(x,y): #Input is two coordinates
    if _matrix[get_entry(list(x)[0], list(x)[1])][get_entry(list(y)[0],list(y)[1])] != 0:
        return True
    else:
        return False

for i in range(caps): #coding in the caps data
    _matrix[get_entry(i,i)][get_entry(i, boundary_points - i - 1)] = 1
    _matrix[get_entry(i, boundary_points - i - 1)][get_entry(i,i)] = -1
    _matrix[get_entry(tangles - 2 - i, boundary_points - i - 1)][get_entry(tangles - 2 - i, i)] = 1
    _matrix[get_entry(tangles - 2 - i,i)][get_entry(tangles - 2 - i, boundary_points - i - 1)] = -1
    for j in range(i): 
        _matrix[get_entry(i,j)][get_entry(i-1,j)] = 1
        _matrix[get_entry(i-1,j)][get_entry(i,j)] = -1
        _matrix[get_entry(i-1,boundary_points-j-1)][get_entry(i,boundary_points-j-1)] = 1
        _matrix[get_entry(i,boundary_points-j-1)][get_entry(i-1,boundary_points-j-1)] = -1

        _matrix[get_entry(tangles-2-i,j)][get_entry(tangles-1-i,j)] = -1
        _matrix[get_entry(tangles-1-i,j)][get_entry(tangles-2-i,j)] = 1
        _matrix[get_entry(tangles-1-i,boundary_points-j-1)][get_entry(tangles-2-i,boundary_points-j-1)] = -1
        _matrix[get_entry(tangles-2-i,boundary_points-j-1)][get_entry(tangles-1-i,boundary_points-j-1)] = 1

#Horizontal lines in the bottom half of the crossing tangles:
for i in range(crossing_tangles):
    for j in range(boundary_points//2):
        _matrix[get_entry(i+caps,j)][get_entry(i+caps-1,j)] = 1
        _matrix[get_entry(i+caps-1,j)][get_entry(i+caps,j)] = -1

#The crossings themselves
for i in range(crossing_tangles):
    for j in range(boundary_points//2):
        if j == abs(sequence[i]):
            _matrix[get_entry(i+caps,boundary_points-j-1)][get_entry(i+caps-1,boundary_points-j)] = -1
            _matrix[get_entry(i+caps-1,boundary_points-j)][get_entry(i+caps,boundary_points-j-1)] = 1
            _matrix[get_entry(i+caps,boundary_points-j)][get_entry(i+caps-1,boundary_points-j-1)] = -1
            _matrix[get_entry(i+caps-1,boundary_points-j-1)][get_entry(i+caps,boundary_points-j)] = 1
        elif j+1 != abs(sequence[i]):
            _matrix[get_entry(i+caps,boundary_points-j-1)][get_entry(i+caps-1,boundary_points-j-1)] = -1
            _matrix[get_entry(i+caps-1,boundary_points-j-1)][get_entry(i+caps,boundary_points-j-1)] = 1



'''
#Print the matrix
for i in range((tangles - 1) * boundary_points):
    for j in range((tangles - 1) * boundary_points):
        if _matrix[i][j] != 0:
            print (get_back(i), get_back(j), _matrix[i][j])  

'''

#Need to see the convention for the crossings but here we will give each crossing a coordinate. Check the convention...
crossing_coord = []
for i in range(len(sequence)):
    if sequence[i] > 0:
        crossing_coord.append((caps - 1 + i + 0.75, boundary_points - sequence[i] - 0.5))
    else:
        crossing_coord.append((caps - 1 + i + 0.25, boundary_points + sequence[i] - 0.5))
#print (crossing_coord)
#Tells if crossing is over or under
def is_positive(x):
    if sequence[x - 1] > 0:
        return True
    else:
        return False

s = []
for i in range(boundary_points+1):
    s.append(i-0.5)
def findsubsets(S,m):
    return set(itertools.permutations(S,m))
a = []
for i in range(len(s)):
    #print("Finding subsets of size", i)
    if 2*i < len(s):
        c = list(findsubsets(s,i))
    else:
        c = list(findsubsets(s,len(s)-i))
        for j in range(len(c)):
            c[j] = set(s) - set(c[j])

    a.append(c)
#print("Finding subsets of siiize", len(s))
b = []
b.append(set(s))
a.append(b)


def sub_bijections_list(b):
    bijlist = []
    for i in range(len(b)):
        for j in range(len(b)):
            bij = list(itertools.izip_longest(b[i],b[j]))
            bijlist.append(bij)
    return bijlist
def modified_sub_bijections_list(b,c):
    bijlist = []
    for i in range(len(b)):
        for j in range(len(c)):
            bij = list(itertools.izip_longest(b[i],c[j]))
            bijlist.append(bij)
    return bijlist
 
def all_bijections (a):
    '''This method is tricky. It takes the input a, which this program automatically generates given a tangle. a is just a list of list of sets, where each list of sets is a list of all subsets of T (where the algebra A(eL) is generated by partial bijections of T). Anytime this method is called, call it as all_bijections(a). It will return for you all the partial bijections of T.'''
    all_bij = []
    for i in range(len(a)):
        all_bij.append(sub_bijections_list(a[i]))
    return all_bij
# The ordered pairs have left coordinates domains and right coordinates images.

def get_index(element): #Inputs an element in the algebra (list of tuples) and outputs the index in the list all bijections. Returns a tuple.
    return (len(element),all_bijections(a)[len(element)].index(element))
def get_element(x, y): #Given a tuple of positive integers returns the corresponding bijection.
    if x >= len(all_bijections(a)) or y >= len(all_bijections(a)[x]):
        print ("nope")
    else:
        return all_bijections(a)[x][y]

# Coding the multiplcation of algebra elements. I.e. the concatenations.

def domain (a): #For domain need to input a list of tuples or all_bijections then 2 indices. Returns a set. 
    domain = []
    for i in range(len(a)):
        domain.append(list(a[i])[0])
    domain.sort()
    return domain
def image (a): #Similar to domain method
    image = []
    for i in range(len(a)):
        image.append(list(a[i])[1])
    image.sort()
    return image

def bsc(x,y): #Short for black strand crossing. This tells us if there is a crossing between two black strands in the algebra strand diagrams.
    if (list(x)[0] < list(y)[0] and list(x)[1] > list(y)[1]) or (list(x)[0] > list(y)[0] and list(x)[1] < list(y)[1]):
        return True
def dbsc(a,b): #Tells if there is a double crossing between two bijections. Input is two bijections. 
    for i in range(len(a)):
        for j in range(len(a)):
            for k in range(len(b)):
                for l in range(len(b)):
                    if bsc(a[i], a[j]) and bsc(b[k], b[l]) and a[i][1] == b[k][0] and a[j][1] == b[l][0]:
                        return True
def asc(a,b): 
    for i in range(len(a)):
        if a[i][0] < a[i][1]:
            if a[i][1] in domain(b):
                k = domain(b).index(a[i][1])
                if b[k][0] > b[k][1]:
                    return True
        if a[i][0] > a[i][1]:
            if a[i][1] in domain(b):
                l = domain(b).index(a[i][1])
                if b[l][0] < b[l][1]:
                    return True

def f2in_permutation_finder(elem, asshole): #This method searches for if any permutation of the list elem is in the list of lists asshole. If it is, the program returns which permutation of elem is in asshole by returning the index of this permutation in itertools.permutations(elem), along with the permutation of elem itself.
    for k in range(len(asshole)):
        if sorted(elem) == sorted(asshole[k]):
            return [1,asshole[k]]
    return [-1,"nope not in here"]

def f2in(elem, asshole): #It returns true if the list elem or any permutations of elem is in the list of lists asshole.
    if f2in_permutation_finder(elem,asshole)[0] == 1:
        return True
    else:
        return False

def f2add(d, diff):
    if not(f2in(d, diff)):
        diff.append(d)
    else:
        diff.remove(f2in_permutation_finder(d,diff)[1])


def f3in_permutation_finder(elem, asshole): #This method searches for if any permutation of the two elements of elem, a list of size two, is in the list of lists asshole. If it is, the program returns which permutation of elem is in asshole by returning the index of this permutation in itertools.permutations(elem), along with the permutation of elem itself.
    for k in range(len(asshole)):
        if [sorted(elem[0]),sorted(elem[1])] == [sorted(list(asshole[k])[0]),sorted(list(asshole[k]))[1]]:
            return [1,asshole[k]]
    return [-1,"nope not in here"]

def f3in(elem, asshole): #It returns true if the list elem or any permutations of elem is in the list of lists asshole.
    if f3in_permutation_finder(elem,asshole)[0] == 1:
        return True
    else:
        return False

def f3add(d, diff):
    if not(f3in(d, diff)):
        diff.append(d)
    else:
        diff.remove(f3in_permutation_finder(d,diff)[1])




def alg_mult (a, b): #imput is two elements of all_bij a[number of elements][bijection]
    product = []
    if image(a) != domain(b):
        return 0
    elif dbsc(a,b):
        return 0
    elif asc(a,b):
        return 0
    else:
        for i in range(len(a)):
            if not( f2in((list(a[i])[0],list(b[i])[1]), product) ):
                product.append((list(a[i])[0],list(b[i])[1]))
            else:
                product.remove(f2in_permutation_finder((list(a[i])[0],list(b[i])[1]), product)[1])
    return product            


#We consider a sum of non-zero elements in the algebera to be given by a list/array of bijections. Here is a method to multiply sums.
# Remember, each bijection is a list of tuples! So the arguments here are lists of lists of tuples! How fun!

def mul_sum(a,b):
    prod_sum = []
    for i in range(len(a)):
        for j in range(len(b)):
            if alg_mult(a[i], b[j]) != 0 and not(f2in(alg_mult(a[i],b[j]) , prod_sum)):
                prod_sum.append(alg_mult(a[i],b[j]))
            elif alg_mult(a[i], b[j]) != 0 and f2in(alg_mult(a[i],b[j]) , prod_sum):
                prod_sum.remove(f2in_permutation_finder(alg_mult(a[i],b[j]) , prod_sum)[1])                
    return prod_sum


"""
What needs to be done is to continue to work with the alg_mult method in order to take into account the modular relationships as discussed.
"""


def is_idempotent_generator(a): #for a algebra generator, returns true if it's idempotent
    for i in range(len(a)):    
        if list(a[i])[0] != list(a[i])[1]:
            return False
    return True

def is_idempotent(a): #for a list of bijections, i.e. general linear combination of algebra generators, returns true if its idempotent
    for i in range(len(a)) :
        if not is_idempotent_generator(a[i]):
            return False
    return True

''' 
def grid_state(b):
    grid_state = []
    for i in range(boundary_points + 1):
        for j in range(len(all_bijections(a)[i])):
            print(i,j)
            if set(image(b)).isdisjoint(set(domain(all_bijections(a)[i][j]))):
                grid_state.append(all_bijections(a)[i][j])
    return grid_state
    '''

#b is any bijection. Returns an array of all the possible bijections that it can be next to in a grid state

#gets the sign sequence at the left edge of tangle i
#if you give me some stupid shitty knot with no crossings I WILL BE ANGRY

def sign_sequence(i): #Returns an array representing the sign sequence along the ith range of alpha curves (the ith red line in the strand diagram). The array starts from the bottom and goes up telling you where it's - (indicated by -1) and where it's + (indicated by 1).
    signs = [0] * boundary_points
    for j in range(boundary_points):
        #print("j",j)
        if 0 <= i-1 < tangles-1 and 0 <= j-1 < boundary_points and signs[j] == 0:
            signs[j] = -1 * _matrix[get_entry(i,j)][get_entry(i-1,j-1)] 
        if 0 <= i-1 < tangles-1 and 0 <= j < boundary_points and signs[j] == 0:        
            signs[j] = -1 * _matrix[get_entry(i,j)][get_entry(i-1,j)] 
        if 0 <= i-1 < tangles-1 and 0 <= j+1 < boundary_points and signs[j] == 0:        
            signs[j] = -1 * _matrix[get_entry(i,j)][get_entry(i-1,j+1)] 

        if 0 <= i+1 < tangles-1 and 0 <= j-1 < boundary_points and signs[j] == 0:
            signs[j] = 1 * _matrix[get_entry(i,j)][get_entry(i+1,j-1)] 
        if 0 <= i+1 < tangles-1 and 0 <= j < boundary_points and signs[j] == 0:
            signs[j] = 1 * _matrix[get_entry(i,j)][get_entry(i+1,j)] 
        if 0 <= i+1 < tangles-1 and 0 <= j+1 < boundary_points and signs[j] == 0:        
            signs[j] = 1 * _matrix[get_entry(i,j)][get_entry(i+1,j+1)] 
            
    
    return signs


def allowable_grid_state(b): #Given b a list of bijections in the form [(-0.5,0.5),(0.5,2.5),etc] that represents one half of a tangle, this program returns an array chock full of all the possible grid states for the entire tangle.
    #TODO Remove points that are empty
    grid_state = []
    points = []
    for i in range(boundary_points+1):
        points.append(i-0.5)
    new_domain = set(points).difference(image(b))
    new_ranges = list(itertools.permutations(points,len(new_domain)))
    for i in range(len(new_ranges)):
        bij = list(itertools.zip_longest(new_domain,new_ranges[i]))
        grid_state.append(bij)

    return grid_state

def alg_diff_generator_modulo(b,i,j): #The mod relations for the differential on the algebra. Returns true if for a list of bijections b, the ith and jth strands cross in such a way that the differential would spit out a zero for this term.
    c = domain(b)
    if b[i][0] < b[j][0] and b[i][1] > b[j][1]:
        for k in range(int(b[i][0]+0.5),int(b[j][0]-0.5)):
            p = k + 0.5
            if p in c:
                l = c.index(p)
                if b[j][1] < b[l][1] < b[i][1]:
                    return True
    elif b[i][0] > b[j][0] and b[i][1] < b[j][1]:
        for k in range(int(b[j][0]+0.5),int(b[i][0]-0.5)):
            p = k + 0.5
            if p in c:
                l = c.index(p)
                if b[i][1] < b[l][1] < b[j][1]:
                    return True

def alg_diff_generator(b): #Differential for a single element (b a bijection). Returns a list of bijections, each element in the list is one of the terms. 
    diff = []
    for i in range(len(b)):
        j = i + 1
        while j < len(b):
            d = []
            for k in range(len(b)):
                d.append(b[k])
            if  bsc(b[i], b[j]):
                num1 = (list(b[i])[0], list(b[j])[1])
                num2 = (list(b[j])[0], list(b[i])[1])
                d.remove(b[i])
                d.append(num1)
                d.remove(b[j])
                d.append(num2)
                if not alg_diff_generator_modulo(b,i,j) and not f2in(d , diff):
                    diff.append(d)
                elif not alg_diff_generator_modulo(b,i,j) and f2in(d , diff):
                    diff.remove(f2in_permutation_finder(d,diff)[1])
            j += 1
    return diff


def alg_diff(b): #For a linear combination of algebra basis elements (b a list of bijection). Returns a list of bijections, each element in the list is one of the terms. .
    diff = []
    for i in range(len(b)):
        war = alg_diff_generator(b[i])
        for j in range(len(war)):
            if not f2in(war[j], diff):
                diff.append(war[j])
            else:
                diff.remove(f2in_permutation_finder(war[j], diff)[1])
    return diff



#Here are all the gradings we need. Might need to take a look again later...
def gradings_alg(b, c): #Takes a bijection and a tangle. And returns an ordered pair of the maslov and alexander gradings.
    num_cross = 0
    right_cross = 0
    left_cross = 0
    for i in range(len(b)):
        j = i + 1
        while j < len(b):
            if bsc(b[i],b[j]):
                num_cross += 1
            j +=1 
        for k in range(boundary_points//2, boundary_points):  
            if min(list(b[i])[0], list(b[i])[1]) < k < max(list(b[i])[0], list(b[i])[1]) and sign_sequence(c)[k] == 1 :
                right_cross += 1
        for l in range(boundary_points//2):
            if min(list(b[i])[0], list(b[i])[1]) < l < max(list(b[i])[0], list(b[i])[1]) and sign_sequence(c)[l] == -1:
                left_cross += 1
    return (num_cross - right_cross,(left_cross - right_cross)/2)
def left_maz(b, c): #Left Mazlov for the grid states.
    lmaz = 0
    for i in range(len(b)):
        j = i + 1
        while j < len(b):
            if bsc(b[i],b[j]):
                lmaz += 1
            j +=1
    for l in range(boundary_points//2):
            if min(list(b[i])[0], list(b[i])[1]) < l < max(list(b[i])[0], list(b[i])[1]) and sign_sequence(c)[l] == -1:
                lmaz += 1
    for k in range(boundary_points):
        if sign_sequence(c)[k] == -1:
            lmaz += 1
    return lmaz
        
def right_maz(b,c): #Right Mazlov for the grid states.
    rmaz = list(gradings_alg(b,c))[0]
    if caps <= c < tangles - caps:
        rmaz += 1
    return rmaz
def grid_state_grading(b, c): #Takes a single partial bijection and a tangle and returns the gradings. Right now only Alexander grading. 
    alex = list(gradings_alg(b,c))[1]
    if caps <= c < tangles - caps: #Check to see if this is correct...
        alex += 1
    for i in range(boundary_points):
        if sign_sequence(c)[i] == -1:
            alex += 1
    return alex
def full_grid_grading(b):
    alex = 0
    if len(b) != tangles:
        print("Not a full state")
    else:
        for i in range(boundary_points):
            alex += grid_state_grading(b[i],i)
    return alex

 # takes a half-integer and gives the alpha or beta curves at that point.
 # u = 0 corresponds to A^0 and will give you the 0th alpha curves on the line containing points with coordinates (0,0) (0,1), etc.
 # u = 0.5 corresponds to B^0
 # u = 1 corresponds to A^1
 # and so on
def alpha_betas_helper(u):
    ev = []
    if u == -1:
        ev.append((u, 0.5*(boundary_points-1)))
        return ev
    elif u == tangles-1:
        ev.append((u, 0.5*(boundary_points-1)))
        return ev
    elif u == tangles-1.5:
        ev.append((u,-0.5))
        ev.append((u,boundary_points-0.5))
        return ev
    elif u == -0.5:
        ev.append((u,-0.5))
        ev.append((u,boundary_points-0.5))
        return ev
    elif 0 <= u < caps-1:
        i = 0
        if (2*u)%2 == 0:
            lim = 2+u
        else:
            lim = 2+u-0.5
        while i < lim:
            ev.append((u,-0.5+i))
            i = i+1
        if (2*u)%2 == 0:
            i = 2+u-1
        else:
            i = 2+u-1-0.5
        while i > -1:
            ev.append((u,boundary_points-0.5-i))
            i = i-1
        return ev
    elif tangles-1-caps < u <= tangles-2:
        wall = tangles-2-u
        i = 0
        if (2*wall)%2 == 0:
            lim = 2+wall
        else:
            lim = 2+wall-0.5

        while i < lim:
            ev.append((u,-0.5+i))
            i = i+1
        if (2*wall)%2 == 0:
            i = 2+wall-1
        else:
            i = 2+wall-1-0.5
        while i > -1:
            ev.append((u,boundary_points-0.5-i))
            i = i-1
        return ev

    elif caps-1 <= u <= tangles-1-caps:
        for i in range(boundary_points+1):
            ev.append((u,i-0.5))
        return ev
    else:
        #print("wtf no ew. index out of bounds. no alpha no beta here.")
        return ev

#for grid states...
def alpha_betas(t): # takes a strand and gives the coordinates of the alpha and beta curves in the strand version of the grid diagram. A tangle is given by the right most coordinate.
    return [alpha_betas_helper(t-1),alpha_betas_helper(t-0.5),alpha_betas_helper(t)]


'''
#for grid states...
def alpha_betas(t): # takes a strand and gives the coordinates of the alpha and beta curves in the strand version of the grid diagram. A tangle is given by the right most coordinate.
    left_coords = []
    mid_coords = []
    right_coords = []
    if caps <= t < tangles - caps:
        for i in range(boundary_points + 1):
            left_coords.append((t - 2, i-0.5))
            mid_coords.append((t-1.5, i- 0.5))
            right_coords.append((t-1, i- 0.5))
        return [left_coords, mid_coords, right_coords]
    elif t == 1:
        mid_coords.append((-0.5,-0.5))
        mid_coords.append((-0.5, boundary_points - 0.5))
        right_coords.append((0,-0.5))
        right_coords.append((0, (boundary_points - 1)/2))
        right_coords.append((0, boundary_points - 0.5))
        left_coords.append((-1,(boundary_points)/2))
        return [left_coords, mid_coords, right_coords]
    elif t == tangles:
        mid_coords.append((tangles - 1.5, -0.5))
        mid_coords.append((tangles - 1.5, boundary_points - 0.5))
        left_coords.append((tangles - 2, -0.5))
        left_coords.append((tangles - 2, (boundary_points- 1)/2))
        left_coords.append((tangles - 2, boundary_points - 0.5))
        return [left_coords, mid_coords]
    elif t < caps:
        if t ==2:
            left_coords = alpha_betas(t-1)[1]
        else:
            left_coords = alpha_betas(t-1)[2]
            for i in range(t):
                mid_coords.append((t-1.5,-0.5 + i)) 
                mid_coords.append((t-1.5, boundary_points - 0.5 - i))
                right_coords.append((t-1,-0.5 + i))
                right_coords.append((t-1, boundary_points hw- 0.5 - i))
        right_coords.append((t-1, (boundary_points-1)/2))
        return [left_coords, mid_coords, right_coords]
    elif t >= tangles - caps:
        right_coords = alpha_betas(t+1)[0]
        for i in range(t):
            mid_coords.append((t-1.5,-0.5 + i))
            mid_coords.append((t-1.5, boundary_points - 0.5 - i))
            left_coords.append((t-1,-0.5 + i))
            left_coords.append((t-1, boundary_points - 0.5 - i))
        left_coords.append((t-1, (boundary_points-1)/2))
        return [left_coords, mid_coords, right_coords]
'''

#HOW TO USE HEEGARD METHODS

#NOTE: The array generated by Heegard diagram methods uses indexing different from the convention. (0)(1) is on the bottom row, unlike on the left column
#looks like
#10 11 12 13
#00 01 02 03
#so that 01=X, 11=O corresponds to the cap
#'' 'O' '' ''
#'' 'X' '' ''
#gives the X's and O's as a grid, [0][0] is the bottom-left, [0][1] is immediately to the right of [0][0].
#gives the left Heegard diagram immediately to the right of i on the grid. For example, left_heegard(0) gives us the first left Heegard diagram immediately following the first cap
def left_heegard(i):
    ss = sign_sequence(i)
    left = [""] * boundary_points
    crosses = False
    crosspoint = -100
    for j in range(len(crossing_coord)):
        if list(crossing_coord[j])[0] == i+0.25:
            crosses = True
            crosspoint = list(crossing_coord[j])[1]
    for ass in range(boundary_points):
        left[ass] = [""] * boundary_points

    for hewler in range(boundary_points):
        if i == tangles-2:
            if ss[hewler] == -1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, hewler)][get_entry(i,k)] != 0:
                        cappair = k
                left[hewler][min(cappair,hewler)] = "O"
            elif ss[hewler] == 1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, hewler)][get_entry(i,k)] != 0:
                        cappair = k
                left[hewler][min(cappair,hewler)] = "X"
        else:
            if sign_sequence(i+1)[hewler] == 0 and ss[hewler] == -1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, hewler)][get_entry(i,k)] != 0:
                        cappair = k
                left[hewler][min(cappair,hewler)] = "O"
            elif sign_sequence(i+1)[hewler] == 0 and ss[hewler] == 1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, hewler)][get_entry(i,k)] != 0:
                        cappair = k
                left[hewler][min(cappair,hewler)] = "X"

            elif ss[hewler] == 1 and abs(crosspoint - hewler) > 0.5:
                left[hewler][hewler] = "X"
            elif ss[hewler] == -1 and abs(crosspoint - hewler) > 0.5:
                left[hewler][hewler] = "O"
            elif ss[hewler] == 1 and crosspoint - hewler == 0.5:
                left[hewler][hewler+1] = "X"
            elif ss[hewler] == -1 and crosspoint - hewler == 0.5:
                left[hewler][hewler+1] = "O"
            elif ss[hewler] == 1 and crosspoint - hewler == -0.5:
                left[hewler][hewler-1] = "X"
            elif ss[hewler] == -1 and crosspoint - hewler == -0.5:
                left[hewler][hewler-1] = "O"
                


    return left

#unlike left_heegard, starts from the BOTTOM RIGHT as (0,0).
#looks like
#13 12 11 10
#03 02 01 00
#so that 01=X, 11=O corresponds to the cap
#'' '' 'O' ''
#'' '' 'X' ''

#gives the right Heegard diagram immediately to the left of i on the grid. For example, right_heegard(0) gives us the first cap
def right_heegard(i):
    ss = sign_sequence(i)



    left = [""] * boundary_points
    crosses = False
    crosspoint = -100
    for j in range(len(crossing_coord)):
        if list(crossing_coord[j])[0] == i-0.25:
            crosses = True
            crosspoint = list(crossing_coord[j])[1]
    for ass in range(boundary_points):
        left[ass] = [""] * boundary_points

    for hewler in range(boundary_points):
        if i != 0:
            if sign_sequence(i-1)[hewler] == 0 and ss[hewler] == -1:
                cappair = -1
                for k in range(boundary_points):
                    #print("k=",k)
                    if _matrix[get_entry(i, hewler)][get_entry(i,k)] != 0:
                        cappair = k
                left[hewler][min(cappair,hewler)] = "X"
            elif sign_sequence(i-1)[hewler] == 0 and ss[hewler] == 1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, hewler)][get_entry(i,k)] != 0:
                        cappair = k
                left[hewler][min(cappair,hewler)] = "O"

            elif ss[hewler] == 1 and abs(crosspoint - hewler) > 0.5:
                left[hewler][hewler] = "O"
            elif ss[hewler] == -1 and abs(crosspoint - hewler) > 0.5:
                left[hewler][hewler] = "X"
            elif ss[hewler] == 1 and crosspoint - hewler == 0.5:
                left[hewler][hewler+1] = "O"
            elif ss[hewler] == -1 and crosspoint - hewler == 0.5:
                left[hewler][hewler+1] = "X"
            elif ss[hewler] == 1 and crosspoint - hewler == -0.5:
                left[hewler][hewler-1] = "O"
            elif ss[hewler] == -1 and crosspoint - hewler == -0.5:
                left[hewler][hewler-1] = "X"
        else:
            if ss[hewler] == -1:
                cappair = -1
                for k in range(boundary_points):
                    #print("k=",k)
                    if _matrix[get_entry(i, hewler)][get_entry(i,k)] != 0:
                        cappair = k
                left[hewler][min(cappair,hewler)] = "X"
            elif ss[hewler] == 1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, hewler)][get_entry(i,k)] != 0:
                        cappair = k
                left[hewler][min(cappair,hewler)] = "O"

    return left

#Prints all the Heegard diagrams from cap to cap. The first list of lists printed corresponds to the right grid on the Heegard diagram for the first cap
def print_heegard():
    print("see code to understand how to interpret output. each tangle has a left grid and a right grid ")
    print("First Cap:")
    for i in range(tangles-1):
        print(right_heegard(i))
        if i != tangles -2:
            print("Tangle between", i, "and", i+1)
        else:
            print("Last Cap:")
        print(left_heegard(i))

#Next thing to do is to generate all the gridstates as bijections. Can do it given an individual tangle?
def halves(t): #Inputs a tangle and returns the halves of each tangle and their bijections.
    left_half = []
    right_half = []
    if t != tangles:
        for i in range(min(len(alpha_betas(t)[0]), len(alpha_betas(t)[1])) + 1):
            if not(t == 0 and i == 0):
                left_half.append(list(modified_sub_bijections_list(list(findsubsets(alpha_betas(t)[0],i)),list(findsubsets(alpha_betas(t)[1],i)))))
        for k in range(min(len(alpha_betas(t)[1]), len(alpha_betas(t)[2])) + 1):
            if not(t == tangles-1 and i == 0):
                right_half.append(modified_sub_bijections_list(list(findsubsets(alpha_betas(t)[1],k)),list(findsubsets(alpha_betas(t)[2],k))))
        return[left_half, right_half]
    else:
        for i in range(min(len(alpha_betas(t)[0]), len(alpha_betas(t)[1])) + 1):
            left_half.append(modified_sub_bijections_list(list(findsubsets(alpha_betas(t)[0],i)),list(findsubsets(alpha_betas(t)[1],i))))
        return[left_half]
#The indexing is thus... First left or right, sie of bijection and then the bijection itself then the elements of the bijection... Need to list it to get the actual images out. 
def gs(t): #Outputs the possible grid states for a tangle.
    state = []
    niv = halves(t)
    for i in range(len(niv[0])):
        k = len(alpha_betas(t)[1]) - i
        while k >= 0:
            for j in range(len(niv[0][i])):
                for h in range(len(niv[1][k])):
                    if set(image(niv[0][i][j])).isdisjoint(domain(niv[1][k][h])):
                        print("yup",i,k,j,h,[niv[0][i][j], niv[1][k][h]])
                        state.append([niv[0][i][j], niv[1][k][h]])
                    else:
                        print("nope")
            k -=1
    return state




def ralg_diff_generator(b): #Reverse of algebra. Need to factor in the modular relationships, Same as del-. Introduces a crossing between two strands that do not have one and sums.
    diff = []
    for i in range(len(b)):
        j = i + 1
        while j < len(b):
            d = []
            for k in range(len(b)):
                d.append(b[k])
            if  not(bsc(b[i], b[j])):
                num1 = (list(b[i])[0], list(b[j])[1])
                num2 = (list(b[j])[0], list(b[i])[1])
                d.remove(b[i])
                d.append(num1)
                d.remove(b[j])
                d.append(num2)
                f2add(d, diff)
            j += 1
        
    return diff

def switch(b,c): #Takes a grid state and reveses that shit like dm. Slide up into my dm's mwah. Edit this to take in a grid state hot and fresh from the method.
    #TODO figure out wat the hell this method is, then mod by F2.
    switch = []
    for i in range(len(b)):
        for j in range(len(c)):
            d = []
            e = []
            for k in range(len(b)):
                d.append(b[k])
            for k in range(len(c)):
                e.append(c[k])
            num1 = (list(b[i])[0], list(c[j])[0])
            num2 = (list(b[i])[1], list(c[j])[1])
            d.remove(b[i])
            d.append(num1)
            e.remove(c[j])
            e.append(num2)
            switch.append([d,e])
        j +=1
    return switch

def orange_on_black_mod(b,i,j,u): #Returns true if the ith and jth strands of bijection b at position u (i.e. strands go from u-0.5 to u) both cross each other as well as an orange strand.
    #for i in range(len(b)):
    if bsc(b[i],b[j]):
        k = int(min(list(b[i])[0], list(b[j])[0]) + 0.5)
        while k < max(list(b[i])[0], list(b[j])[0]):
            if (2*u)%2 == 0:
                for ass in range(int(min(list(b[i])[1], list(b[j])[1])), int(max(list(b[i])[1], list(b[j])[1]))):
                    if _matrix[get_entry(int(u-1),k)][get_entry(int(u),ass)] != 0:
                        return True
            else:
                for ass in range(int(min(list(b[i])[1], list(b[j])[1])), int(max(list(b[i])[1], list(b[j])[1]))):
                    if _matrix[get_entry(int(u-0.5),k)][get_entry(int(u+0.5),ass)] != 0:
                        return True

            k = k+1
    return False



def d_plus_half(b, u): #Takes a half of a grid state and gives d_plus. b is the bijection, u is the position of the bijections range (so that A3-->B3 has u=3.5)
    diff = []
    for i in range(len(b)):
        j = i + 1
        while j < len(b):
            d = []
            for k in range(len(b)):
                d.append(b[k])
            if  bsc(b[i], b[j]):
                num1 = (list(b[i])[0], list(b[j])[1])
                num2 = (list(b[j])[0], list(b[i])[1])
                d.remove(b[i])
                d.append(num1)
                d.remove(b[j])
                d.append(num2)
                if not (alg_diff_generator_modulo(b,i,j) or orange_on_black_mod(b,i,j,u) or f2in(d , diff)):     #NOTE: Double pointed homology will require modifying the modulo.
                    diff.append(d)
                elif not (alg_diff_generator_modulo(b,i,j) or orange_on_black_mod(b,i,j,u)) and f2in(d , diff):
                    diff.remove(f2in_permutation_finder(d,diff)[1])
            j += 1
    return diff

def alg_diff_generator_modulo_stupid_plus(b,i,j,u): 
    bool = False
    arr = []

    if (2*u)%2 == 1:
        u = u + 0.5

    for l in range(int(min(b[i][1], b[j][1])+0.5),int(max(b[i][1], b[j][1])+0.5)):
        arr.append(l)
        for m in range(3):
            if (0 <= u < tangles-1 and 0 <= l < boundary_points) and (0 <= u-1 < tangles-1 and 0 <= l-1+m < boundary_points):
                if _matrix[int(get_entry(u,l))][int(get_entry(u-1,l-1+m))]:
                    bool = True
    for k in range(len(b)):
        if b[i][1] < b[k][1] < b[j][1] and (b[k][0] > max(b[j][0],b[i][0]) or  b[k][0] < min(b[j][0],b[i][0])):
            bool = True
        if min(b[i][1], b[j][1]) < b[k][1] < max(b[i][1], b[j][1]):
            arr.remove(l)
    if len(arr) > 0:
        bool = True

    if bool:
        return bool
    else:
        return alg_diff_generator_modulo(b,i,j)

def alg_diff_generator_modulo_stupid_minus(b,i,j,u):
    bool = False
    arr = []
    for l in range(int(min(b[i][0], b[j][0])+0.5),int(max(b[i][0], b[j][0])+0.5)):
        arr.append(l)
        for m in range(3):
            if (0 <= u-1 < tangles-1 and 0 <= l < boundary_points) and (0 <= u < tangles-1 and 0 <= l-1+m < boundary_points):
                if _matrix[get_entry(u-1,l)][get_entry(u,l-1+m)] != 0:
                  bool = True

    for k in range(len(b)):
        if b[i][0] < b[k][0] < b[j][0] and (b[k][1] > max(b[j][1],b[i][1]) or  b[k][1] < min(b[j][1],b[i][1])):
            bool = True
        if min(b[i][0], b[j][0]) < b[k][0] < max(b[i][0], b[j][0]):
            arr.remove(l)
    if len(arr) > 0:
        bool = True

    if bool:
        return bool
    else:
        return alg_diff_generator_modulo(b,i,j)



def d_plus_half_stupid(b,u):
    diff = []
    for i in range(len(b)):
        j = i + 1
        while j < len(b):
            d = []
            for k in range(len(b)):
                d.append(b[k])
            if  bsc(b[i], b[j]):
                num1 = (list(b[i])[0], list(b[j])[1])
                num2 = (list(b[j])[0], list(b[i])[1])
                d.remove(b[i])
                d.append(num1)
                d.remove(b[j])
                d.append(num2)
                if not (alg_diff_generator_modulo_stupid_plus(b,i,j,u) or orange_on_black_mod(b,i,j,u) or f2in(d , diff)):     #NOTE: Double pointed homology will require modifying the modulo.
                    diff.append(d)
                elif not (alg_diff_generator_modulo_stupid_plus(b,i,j,u) or orange_on_black_mod(b,i,j,u)) and f2in(d , diff):
                    diff.remove(f2in_permutation_finder(d,diff)[1])
            j += 1
    return diff


def d_minus_half_stupid(b, u):
    diff = []
    for i in range(len(b)):
        j = i + 1
        while j < len(b):
            d = []
            for k in range(len(b)):
                d.append(b[k])
            if not(bsc(b[i], b[j])):
                num1 = (list(b[i])[0], list(b[j])[1])
                num2 = (list(b[j])[0], list(b[i])[1])
                d.remove(b[i])
                d.append(num1)
                d.remove(b[j])
                d.append(num2)
                if not(alg_diff_generator_modulo_stupid_minus(b,i,j,u) or orange_on_black_mod(b,i,j,u) or f2in(d , diff)):
                    diff.append(d)
                elif not(alg_diff_generator_modulo_stupid_minus(b,i,j,u) or orange_on_black_mod(b,i,j,u)) and f2in(d , diff):
                    diff.remove(f2in_permutation_finder(d , diff)[1])
            j += 1
        
    return diff



def d_minus_half(b, u): #Takes a half of a grid state and gives d_minus. b is the bijection, u is the position of the bijections range (so that A3-->B3 has u=3.5)
    diff = []
    for i in range(len(b)):
        j = i + 1
        while j < len(b):
            d = []
            for k in range(len(b)):
                d.append(b[k])
            if not(bsc(b[i], b[j])):
                num1 = (list(b[i])[0], list(b[j])[1])
                num2 = (list(b[j])[0], list(b[i])[1])
                d.remove(b[i])
                d.append(num1)
                d.remove(b[j])
                d.append(num2)
                if not(alg_diff_generator_modulo(b,i,j,u) or orange_on_black_mod(b,i,j,u) or f2in(d , diff)):
                    diff.append(d)
                elif not(alg_diff_generator_modulo(b,i,j,u) or orange_on_black_mod(b,i,j,u)) and f2in(d , diff):
                    diff.remove(f2in_permutation_finder(d , diff)[1])
            j += 1
        
    return diff

def to_simple_strands(b): #Given a grid state returns it like simple strands ie [((0, -0.5), (0.5, 4.5))] maps to [(-0.5,4.5)]
    left = []
    for i in range(len(b)):
        left.append((b[i][0][1], b[i][1][1]))
    return left

def from_simple_strands(b,u): #u is the rightmost position. Given [(-0.5,4.5),(0.5,0.5)], 0.5 recover [((0, -0.5), (0.5, 4.5)), ((0,0.5),(0.5,0.5))].
    starbs = []
    for i in range(len(b)):
        starbs.append(((u-0.5,b[i][0]),(u,b[i][1])))
    return starbs

#FOR GENERATORS (I.e. GRID STATES, TWO HALVES)

def d_plus_generator(b):   #Takes two halves. Grid state is inputted like [[((0, -0.5), (0.5, 4.5))], [((0.5, 0.5), (1, 0.5)), ((0.5, 5.5), (1, 4.5)), ((0.5, -0.5), (1, 3.5))]]
    autre = []
    arr = d_plus_half(to_simple_strands(b[1]),b[1][0][1][0])
    for i in range(len(arr)):
        f2add(from_simple_strands(arr[i],b[1][0][1][0]) , autre)
    for j in range(len(autre)):
        autre[j] = [b[0],autre[j]]
    return autre


def d_minus_generator(b):  #Takes two halves. Grid state is inputted like [[((0, -0.5), (0.5, 4.5))], [((0.5, 0.5), (1, 0.5)), ((0.5, 5.5), (1, 4.5)), ((0.5, -0.5), (1, 3.5))]]
#Output is array of two-halves.
    autre = []
    arr = d_minus_half(to_simple_strands(b[0]),b[1][0][0][0])
    for i in range(len(arr)):
        f2add(from_simple_strands(arr[i],b[1][0][0][0]) , autre)    
    for j in range(len(autre)):
        autre[j] = [autre[j],b[1]]
    return autre

def d_m_modulo(left,right,i,j,u):
    #jth strand of right and ith strand of left.
    #method will tell me if the strands do the weird mod thing defined in petkova paper (3rd and 4th mod relations, the first two are coded as the stupid plus and minus mods for the stupid plus and stupid minus methods). thanks.
    if left[i][1] < right[j][0]:
        for k in range(len(left)):
            if left[k][1] > left[i][1] and left[k][0] > left[i][0] and left[k][1] < right[j][0]:
                return True
        for k in range(len(right)):
            if left[i][1] < right[k][0] < right[j][0]:
                if bsc(right[k], right[j]):
                    return True
        l = left[i][1]+0.5
        m = left[i][1]-0.5
        while l < boundary_points:
            m = l - 1
            while m < l+2 and m > left[i][0]: 
                if(_matrix[int(get_entry(u,l))][int(get_entry(u-1,m))]) != 0:
                    return True
                m += 1
            l += 1

        l = right[j][1]+0.5
        m = right[j][1]-0.5
        while l < boundary_points:
            m = l - 1 
            while m < l+2 and m < right[j][0]: 
                if(_matrix[int(get_entry(u,l))][int(get_entry(u-1,m))]) != 0:
                    return True
                m += 1
            l += 1

    elif left[i][1] > right[j][0]: 
        for k in range(len(left)):
            if left[k][1] < left[i][1] and left[k][0] < left[i][0] and left[k][1] > right[j][0]:
                return True
        for k in range(len(right)):
            if left[i][1] > right[k][0] > right[j][0]:
                if bsc(right[k], right[j]):
                    return True

        l = left[i][1]-0.5
        m = left[i][1]+0.5
        while l > -1:
            m = l + 1
            while m > l-2 and m < left[i][0]: 
                if(_matrix[int(get_entry(u,l))][int(get_entry(u-1,m))]) != 0:
                    return True
                m -= 1
            l -= 1

        l = right[j][1]-0.5
        m = right[j][1]+0.5
        while l > -1:
            m = l + 1 
            while m > l-2 and m > right[j][0]: 
                if(_matrix[int(get_entry(u,l))][int(get_entry(u-1,m))]) != 0:
                    return True
                m -= 1
            l -= 1

    else:
        print("ERROR: Invalid grid state. One of our beta curves is associated with more than one alpha curve.")

def d_m_generator(b): #Takes two halves. Grid state is inputted like [[((0, -0.5), (0.5, 4.5))], [((0.5, 0.5), (1, 0.5)), ((0.5, 5.5), (1, 4.5)), ((0.5, -0.5), (1, 3.5))]]
    #Outputs an array of two-halves (ie array of grid states.)
    diff = []
    u = list(list(b[0][0])[0])[0] + 1
    left = to_simple_strands(b[0])
    right = to_simple_strands(b[1])
    #Now left and right are arrays of grid states
    leftern = []
    rightern = []
    leftern = d_plus_half_stupid(left,u-0.5)
    rightern = d_minus_half_stupid(right,u)


        #Tricky cuz u have to permute b[0] and b[1] not b.
    for i in range(len(leftern)):
        d = [from_simple_strands(leftern[i], u-0.5) , from_simple_strands(right,u)]
        f3add(d, diff)

    for i in range(len(rightern)):
        d = [from_simple_strands(left,u-0.5),from_simple_strands(rightern[i], u)]
        f3add(d,diff)

    for i in range(len(left)):
        for j in range(len(right)):
            if not d_m_modulo(left,right,i,j,u):
                lefterino = left
                righterino = right
                list(lefterino[i])[1] = right[j][0]
                list(righterino[j])[0] = left[i][1]
                lefterino = tuple(lefterino)
                righterino = tuple(righterino)

                d = [from_simple_strands(lefterino,u-0.5),from_simple_strands(righterino,u)]
                f3add(d, diff)

                #TODO Do the f3-in and add to diff the value [from_simple_strands(lefterino,u-0.5),f_s_s(righterino,u)]


    return diff



'''
#Given two bijections l and r (for the left and right grid states of a generator of CT(T) for a given elementary tangle), spits out the differential computed using Heegard diagrams 
def grid_state_differential_heegard_generator(l,r):
    for i in range(len(all_bijections)):

NOTE: When counting rectangles, make sure to use alpha_betas to see which circles there actually are. Heegard methods only give you the positions of the Xs and Os, not how big the squares are
'''



'''
for i in range(tangles):
    print(i)
    print(alpha_betas(i))
'''
#print(gs(1)) #[[((0, -0.5), (0.5, 4.5)), ((0,0.5),(0.5, 0.5))], [((0.5, 5.5), (1, 4.5)), ((0.5, -0.5), (1, 3.5))]]
print(boundary_points)
print(d_m_generator([[((0, -0.5), (0.5, 9.5)), ((0,5.5),(0.5,7.5))], [((0.5, 0.5), (1, 0.5)), ((0.5, 5.5), (1, 4.5)), ((0.5, -0.5), (1, 3.5))]]))




