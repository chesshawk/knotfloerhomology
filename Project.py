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

sequence = [-1,2] #This is where you tell the program what crossings are present in the knot's tangle decomposition
#If any element in the above array is zero I will find you and I will kill you. 

#supposing notation is that i corresponds to a short crossing between i and i + 1

#if max = i then number of boundary points is 2(i+1)
boundary_points = 2*(get_max(sequence)+1)
print(boundary_points)

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



#Print the matrix
for i in range((tangles - 1) * boundary_points):
    for j in range((tangles - 1) * boundary_points):
        if _matrix[i][j] != 0:
            print (get_back(i), get_back(j), _matrix[i][j])  

#Need to see the convention for the crossings but here we will give each crossing a coordinate. Check the convention...
crossing_coord = []
for i in range(len(sequence)):
    if sequence[i] > 0:
        crossing_coord.append((caps - 1 + i + 0.75, boundary_points - sequence[i] - 0.5))
    else:
        crossing_coord.append((caps - 1 + i + 0.25, boundary_points + sequence[i] - 0.5))
print (crossing_coord)
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
            bij = list(itertools.zip_longest(b[i],b[j]))
            bijlist.append(bij)
    return bijlist
 
def all_bijections (a):
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
    return domain
def image (a): #Similar to domain method
    image = []
    for i in range(len(a)):
        image.append(list(a[i])[1])
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
            product.append((list(a[i])[0],list(b[i])[1]))
    return product            


#We consider a sum of non-zero elements in the algebera to be given by a list/array of bijections. Here is a method to multiply sums.
# Remember, each bijection is a list of tuples! So the arguments here are lists of lists of tuples! How fun!

def mul_sum(a,b):
    prod_sum = []
    for i in range(len(a)):
        for j in range(len(b)):
            if alg_mult(a[i], b[j]) != 0:
                prod_sum.append(alg_mult(a[i],b[j]))
    return prod_sum

"""
What needs to be done is to continue to work with the alg_mult method in order to take into account the modular relationships as discussed.
"""


def is_idempotent_generator(a): #for a generator
    for i in range(len(a)):    
        if list(a[i])[0] != list(a[i])[1]:
            return False
    return True

def is_idempotent(a): #for a list of bijections 
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
#if you give me some stupid shit with no crossings I WILL BE ANGRY
#This is still slow as balls but I'm trying

def sign_sequence(i):
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


def allowable_grid_state(b):
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

def alg_diff_modulo(b,i,j):
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

def alg_diff(b): #Differential for a single element. Need to code in the modular relations. 
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
                if not alg_diff_modulo(b,i,j):
                    diff.append(d)
            j += 1
    return diff

print (alg_diff([(-0.5,3.5), (0.5, 1.5), (2.5, 0.5)]))   

print("rt")
for i in range(tangles-1):
    print(i,sign_sequence(i))

