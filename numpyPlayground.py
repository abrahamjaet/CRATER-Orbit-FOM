import numpy as np 
import orbitDict 
import matplotlib.pyplot as plt 

##### array addition/subtration #######
# arr = np.array([1,1,1])
# print(arr)
# arr = arr - 1
# print(arr) 

######## Orbit simulation testing ########### 

## before ##
arr_points = np.array([[1,1,1],[2,2,2],[3,3,3]])
arr_sat_positions = np.array([[0,1,2],[3,4,5],[6,7,8]])
arr_sc_to_point = arr_points - arr_sat_positions
r_moon = np.array([0,0,0])
r_moonToPoint = arr_points - r_moon 
r_moonToSat = arr_sat_positions - r_moon 
angle1 = np.arccos((np.dot(arr_sc_to_point[0],r_moonToPoint[0])/(np.linalg.norm(arr_sc_to_point[0])*np.linalg.norm(r_moonToPoint[0]))))
angle1_2 = np.arccos((np.dot(arr_sc_to_point[1],r_moonToPoint[1])/(np.linalg.norm(arr_sc_to_point[1])*np.linalg.norm(r_moonToPoint[1]))))
angle1_3 = np.arccos((np.dot(arr_sc_to_point[2],r_moonToPoint[2])/(np.linalg.norm(arr_sc_to_point[2])*np.linalg.norm(r_moonToPoint[2]))))
angle2 = np.arccos((np.dot(r_moonToSat[0],r_moonToPoint[0])/(np.linalg.norm(r_moonToSat[0])*np.linalg.norm(r_moonToPoint[0]))))

## using vectors ## 
# a = np.sum(arr_sc_to_point*r_moonToPoint,axis=1) # this is the dot product equivalent in numpy for row-wise dot product 

################# numerator test #####################

# print("numberator test") # passed 

# print(np.sum(arr_sc_to_point*r_moonToPoint,axis=1))
# print(np.dot(arr_sc_to_point[0],r_moonToPoint[0]))
# print(np.dot(arr_sc_to_point[1],r_moonToPoint[1]))
# print(np.dot(arr_sc_to_point[2],r_moonToPoint[2]))
part1 = np.sum(arr_sc_to_point*r_moonToPoint,axis=1)

################ denominator test part 1 #####################

# print('denominator part 1 test') # passed 
part2 = np.sum(np.abs(arr_sc_to_point)**2,axis=-1)**(1./2)
# print('denominator part 1 test') # passed 
# print(part2)
# print(np.linalg.norm(arr_sc_to_point[0]))
# print(np.linalg.norm(arr_sc_to_point[1]))
# print(np.linalg.norm(arr_sc_to_point[2]))

################ denominator test part 2 #####################


part3 = np.sum(np.abs(r_moonToPoint)**2,axis=1)**(1./2)
# print('denominator part 2 test') # passed 
# print(part3)
# print(np.linalg.norm(r_moonToPoint[0]))
# print(np.linalg.norm(r_moonToPoint[1]))
# print(np.linalg.norm(r_moonToPoint[2]))

################ denominator test #####################
part4 = part2*part3 
# print('denominator test')
# print(part4)
# print(np.linalg.norm(arr_sc_to_point[0])*np.linalg.norm(r_moonToPoint[0]))
# print(np.linalg.norm(arr_sc_to_point[1])*np.linalg.norm(r_moonToPoint[1]))
# print(np.linalg.norm(arr_sc_to_point[2])*np.linalg.norm(r_moonToPoint[2]))


################ division test #####################
part5 = part1 / part4
# print(part5)
# print('division test')
# print(np.dot(arr_sc_to_point[0],r_moonToPoint[0])/(np.linalg.norm(arr_sc_to_point[0])*np.linalg.norm(r_moonToPoint[0])))
# print(np.dot(arr_sc_to_point[1],r_moonToPoint[1])/(np.linalg.norm(arr_sc_to_point[1])*np.linalg.norm(r_moonToPoint[1])))
# print(np.dot(arr_sc_to_point[2],r_moonToPoint[2])/(np.linalg.norm(arr_sc_to_point[2])*np.linalg.norm(r_moonToPoint[2])))

################ arccos test #####################
part6 = np.arccos(part5)
# print('arccos test')
# print(part6)
# print(np.arccos(np.dot(arr_sc_to_point[0],r_moonToPoint[0])/(np.linalg.norm(arr_sc_to_point[0])*np.linalg.norm(r_moonToPoint[0]))))
# print(np.arccos(np.dot(arr_sc_to_point[1],r_moonToPoint[1])/(np.linalg.norm(arr_sc_to_point[1])*np.linalg.norm(r_moonToPoint[1]))))
# print(np.arccos(np.dot(arr_sc_to_point[2],r_moonToPoint[2])/(np.linalg.norm(arr_sc_to_point[2])*np.linalg.norm(r_moonToPoint[2]))))



################ integration test #####################
print(part6)
print(angle1)
print(angle1_2)
print(angle1_3)