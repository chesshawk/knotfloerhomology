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
sequence = [1] #This is where you tell the program what crossings are present in the knot's tangle decomposition
#If any element in the above array is zero I will find you and I will kill you. 

#supposing notation is that i corresponds to a short crossing between i and i + 1

#if max = i then number of boundary points is 2(i+1)
boundary_points = 2*(get_max(sequence)+1)
print(boundary_points)

tangles = len(sequence) + boundary_points
print(tangles)

#further tangle decomposition

crossing_tangles = len(sequence) #number of tangles with a cross

caps = len(sequence) + 1

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
for i in range(boundary_points):
    s.append(i)
def findsubsets(S,m):
    return set(itertools.permutations(S,m))
a = []
for i in range(len(s) + 1):
    a.append(list(findsubsets(s,i)))

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
"""
need code to parse bijections better
"""
# Coding the multiplcation of algebra elements. I.e. the concatenations.

def domain (a): #For domain need to input a list of tuples or all_bijections then 2 indices. Returns a set. 
    domain = []
    for i in range(len(a)):
        domain.append(list(a[i])[0])
    domain = set(domain)
    return domain
def image (a): #Similar to domain method
    image = []
    for i in range(len(a)):
        image.append(list(a[i])[1])
    image = set(image)
    return image
"""
rewrite the above methods so that you do not need to index the full bijections list
"""

print (all_bijections(a)[2][4])
print (domain(all_bijections(a)[2][4]))
print (image(all_bijections(a)[2][4]))

def alg_mult (a, b): #imput is two elements of all_bij a[number of elements][bijection]
    product = []
    if len(a) != len(b):
        return 0
    elif image(a) != domain(b):
        return 0
    else:
        for i in range(len(a)):
            product.append((list(a[i])[0],list(b[i])[1]))
    return product
               
print (all_bijections(a)[2][0])
print (all_bijections(a)[2][1])
print (alg_mult(all_bijections(a)[2][0], all_bijections(a)[2][1]))
print (all_bijections(a)[2].index(alg_mult(all_bijections(a)[2][0], all_bijections(a)[2][1])))

def get_index(element): #Inputs an element in the algebra (list of tuples) and outputs the index in the list all bijections. Returns a tuple.
    return (len(element),all_bijections(a)[len(element)].index(element))
print(get_index(alg_mult(all_bijections(a)[2][0], all_bijections(a)[2][1])))
print (get_index([(0,1),(1,2),(2,3)]))
def get_element(x, y): #Given a tuple of positive integers returns the corresponding bijection.
    if x >= len(all_bijections(a)) or y >= len(all_bijections(a)[x]):
        print ("nope")
    else:
        return all_bijections(a)[x][y]
print (get_element(3,91))



    



    
    