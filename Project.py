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

sequence = [1,-2] #This is where you tell the program what crossings are present in the knot's tangle decomposition
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
    print(num_cross, right_cross, left_cross)
    return (num_cross - right_cross,(left_cross - right_cross)/2)


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
        while i < 2+u:
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
        while i < 2+wall:
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
                right_coords.append((t-1, boundary_points - 0.5 - i))
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

    for fuck in range(boundary_points):
        if i == tangles-2:
            if ss[fuck] == -1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, fuck)][get_entry(i,k)] != 0:
                        cappair = k
                left[fuck][min(cappair,fuck)] = "O"
            elif ss[fuck] == 1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, fuck)][get_entry(i,k)] != 0:
                        cappair = k
                left[fuck][min(cappair,fuck)] = "X"
        else:
            if sign_sequence(i+1)[fuck] == 0 and ss[fuck] == -1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, fuck)][get_entry(i,k)] != 0:
                        cappair = k
                left[fuck][min(cappair,fuck)] = "O"
            elif sign_sequence(i+1)[fuck] == 0 and ss[fuck] == 1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, fuck)][get_entry(i,k)] != 0:
                        cappair = k
                left[fuck][min(cappair,fuck)] = "X"

            elif ss[fuck] == 1 and abs(crosspoint - fuck) > 0.5:
                left[fuck][fuck] = "X"
            elif ss[fuck] == -1 and abs(crosspoint - fuck) > 0.5:
                left[fuck][fuck] = "O"
            elif ss[fuck] == 1 and crosspoint - fuck == 0.5:
                left[fuck][fuck+1] = "X"
            elif ss[fuck] == -1 and crosspoint - fuck == 0.5:
                left[fuck][fuck+1] = "O"
            elif ss[fuck] == 1 and crosspoint - fuck == -0.5:
                left[fuck][fuck-1] = "X"
            elif ss[fuck] == -1 and crosspoint - fuck == -0.5:
                left[fuck][fuck-1] = "O"
                


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

    for fuck in range(boundary_points):
        if i != 0:
            if sign_sequence(i-1)[fuck] == 0 and ss[fuck] == -1:
                cappair = -1
                for k in range(boundary_points):
                    #print("k=",k)
                    if _matrix[get_entry(i, fuck)][get_entry(i,k)] != 0:
                        cappair = k
                left[fuck][min(cappair,fuck)] = "X"
            elif sign_sequence(i-1)[fuck] == 0 and ss[fuck] == 1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, fuck)][get_entry(i,k)] != 0:
                        cappair = k
                left[fuck][min(cappair,fuck)] = "O"

            elif ss[fuck] == 1 and abs(crosspoint - fuck) > 0.5:
                left[fuck][fuck] = "O"
            elif ss[fuck] == -1 and abs(crosspoint - fuck) > 0.5:
                left[fuck][fuck] = "X"
            elif ss[fuck] == 1 and crosspoint - fuck == 0.5:
                left[fuck][fuck+1] = "O"
            elif ss[fuck] == -1 and crosspoint - fuck == 0.5:
                left[fuck][fuck+1] = "X"
            elif ss[fuck] == 1 and crosspoint - fuck == -0.5:
                left[fuck][fuck-1] = "O"
            elif ss[fuck] == -1 and crosspoint - fuck == -0.5:
                left[fuck][fuck-1] = "X"
        else:
            if ss[fuck] == -1:
                cappair = -1
                for k in range(boundary_points):
                    #print("k=",k)
                    if _matrix[get_entry(i, fuck)][get_entry(i,k)] != 0:
                        cappair = k
                left[fuck][min(cappair,fuck)] = "X"
            elif ss[fuck] == 1:
                cappair = -1
                for k in range(boundary_points):
                    if _matrix[get_entry(i, fuck)][get_entry(i,k)] != 0:
                        cappair = k
                left[fuck][min(cappair,fuck)] = "O"

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
'''
#Next thing to do is to generate all the gridstates as bijections. Can do it given an individual tangle?
def halves(t): #Inputs a tangle and returns the halves of each tangle and their bijections.
    left_half = []
    right_half = []
    if t != tangles:
        for i in range(min(len(alpha_betas(t)[0]), len(alpha_betas(t)[1])) + 1):
            left_half.append(list(modified_sub_bijections_list(list(findsubsets(alpha_betas(t)[0],i)),list(findsubsets(alpha_betas(t)[1],i)))))
        for k in range(min(len(alpha_betas(t)[1]), len(alpha_betas(t)[2])) + 1):
            right_half.append(modified_sub_bijections_list(list(findsubsets(alpha_betas(t)[1],k)),list(findsubsets(alpha_betas(t)[2],k))))
        return[left_half, right_half]
    else:
        for i in range(min(len(alpha_betas(t)[0]), len(alpha_betas(t)[1])) + 1):
            left_half.append(modified_sub_bijections_list(list(findsubsets(alpha_betas(t)[0],i)),list(findsubsets(alpha_betas(t)[1],i))))
        return[left_half]
#The indexing is thus... First left or right, sie of bijection and then the bijection itself then the elements of the bijection... Need to list it to get the actual images out. 
def gs(t): #Outputs the possible grid states for a tangle.
    state = []
    for i in range(len(halves(t)[0])):
        k = len(alpha_betas(t)[1]) - i
        while k >= 0:
            for j in range(len(halves(t)[0][i])):
                for h in range(len(halves(t)[1][k])):
                    if set(image(halves(t)[0][i][j])).isdisjoint(domain(halves(t)[1][k][h])):
                        print("yup")
                        state.append([halves(t)[0][i][j], halves(t)[1][k][h]])
            k -=1
    return state
print(gs(1)[2])

'''


'''
#Given two bijections l and r (for the left and right grid states of a generator of CT(T) for a given elementary tangle), spits out the differential computed using Heegard diagrams 
def grid_state_differential_heegard_generator(l,r):
    for i in range(len(all_bijections)):
'''



print("caps",caps)
print("tangles",tangles)
print("ryounday",boundary_points)

for i in range(tangles):
    print(i)
    print(alpha_betas(i))


#print_heegard()






