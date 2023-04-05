import matplotlib.pyplot as plt
import numpy as np
import orbitDict
import time 
#from csltk.utilities import System
import math 

def reorder(lst, newIdx):
    pos = newIdx #lst.index(first)
    return lst[pos:] + lst[:pos]

start = time.time()

## inputs are a design 
## takes in list of orbits with corresponding number of sats per orbit 
numOrbits = 1 
# ax = plt.axes(projection = '3d')
# sys = System(mu=0.01215058560962404, lstar=389703.2648292776, tstar=382981.2891290545)
# ax = sys.plot_system()
ax = plt.axes(projection ="3d") 
###############################################
################# Constants ###################
###############################################

TU = 382981.2891290545 # seconds/TU
LU = 389703 # 1 LU (distance from Earth to Moon in km)
moon_rad = 1740/LU # radius of the Moon in [km]
moonPos = np.array([0.98784942,0,0]) # LU 
lat_spacing = 10 # [deg]
max_points = 12 
latitude = np.linspace(90,-90,19)
gap = 360/max_points 
synodicPeriod = 29.523 # days 
synodicS = synodicPeriod*24*60*60 # s 
sToHr = (1/(60*60))

################################################
############## Grid point Calcluation ##########
################################################

grid_point_coordinates = np.genfromtxt('FOMLoad.csv', delimiter=',') # x, y, z, r,theta,phi,lat,long
grid_point_coordinates[:,0] = grid_point_coordinates[:,0] + 0.98784942 
grid_point_coordinates = grid_point_coordinates[:,0:3]
totalPoints = len(grid_point_coordinates[:,0])

# plot grid points 
ax = plt.axes(projection ="3d") 
ax.scatter3D(grid_point_coordinates[:,0],grid_point_coordinates[:,1],grid_point_coordinates[:,2])
# plt.show()

#################################################
############# Orbit Simulation #################
#################################################


ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
filename = 'trajectories/orbit1.dat'
orb1 = orbitDict.chosenOrbits.load(filename)
orbits = [orb1]
numSats = [1]

periodTU = orb1.T # TU 
periodDays = periodTU*TU*(1/60)*(1/60)*(1/24) # [days]
loops = round(synodicS/100)
rows = totalPoints
cols = loops
rotate = 360/(synodicS/100) # deg to rotate per loop 
coverageMain = np.zeros((rows,cols)) # each row corresponds to a grid point, each col corresponds to a time 
counter = 0 
# go through each orbit


a = 1

for orb in orbits:
    
    coverage = np.zeros((rows,cols)) # rows = points, cols = time 
    timeCounterCoverage = 0 # restart time counter goes with loops 
    counterLoops = 0 # restarting orbits counter 
   
    if not numSats[counter] == 0:
        for i in range(loops): 
            
            # allocate sat position 
            if i >= (len(orb.x)*counterLoops + len(orb.x)):
                counterLoops = counterLoops + 1
            i = i - len(orb.x)*counterLoops
            sat_xPos = orb.x[i]
            sat_yPos = orb.y[i]
            sat_zPos = orb.z[i]
            
            ############# new stuff ##################

            r_spacecraft = [sat_xPos, sat_yPos, sat_zPos] # vector from center of earth to satellite
            r_spacecraftToPoint = np.zeros((len(grid_point_coordinates[:,0]), 3))
            r_spacecraftToPoint[:,0] = grid_point_coordinates[:,0] - sat_xPos
            r_spacecraftToPoint[:,1] = grid_point_coordinates[:,1] - sat_yPos
            r_spacecraftToPoint[:,2] = grid_point_coordinates[:,2] - sat_zPos
            r_moonToPoint = np.zeros((len(grid_point_coordinates[:,0]), 3))
            r_moonToPoint[:,0] = grid_point_coordinates[:,0] - moonPos[0]
            r_moonToPoint[:,1] = grid_point_coordinates[:,1] - moonPos[1]
            r_moonToPoint[:,2] = grid_point_coordinates[:,2] - moonPos[2]
            r_moonToSat = [r_spacecraft[0]-moonPos[0],r_spacecraft[1]-moonPos[1],r_spacecraft[2]-moonPos[2]] # vector from moon center to sat 
            
            r_moonToSat = np.array(r_moonToSat)
            
            
            

            
            part1 = np.sum(np.array(r_spacecraftToPoint)*r_moonToPoint,axis=1)
            part2 = np.sum(np.abs(r_spacecraftToPoint)**2,axis=-1)**(1./2)
            part3 = np.sum(np.abs(r_moonToPoint)**2,axis=1)**(1./2)
            part4 = part2*part3 
            part5 = part1 / part4
            angle1 = np.arccos(part5)
            

            part1B = np.sum(r_moonToSat*r_moonToPoint,axis=1)
            part2B = np.sum(np.abs(r_moonToSat)**2,axis=-1)**(1./2)
            part3B = np.sum(np.abs(r_moonToPoint)**2,axis=1)**(1./2)
            part4B = part2B*part3B 
            part5B = part1B / part4B
            angle2 = np.arccos(part5B)
            


            idxAngle1 = np.where(angle1 > np.pi/2) # gives us the indices that are
            idxAngle2 = np.where(angle2 < np.pi/2) # gives us the indices that are
            idx = np.intersect1d(idxAngle1,idxAngle2)
            
            
            coverage[idx,timeCounterCoverage] = 1
            # print(coverage[idx,timeCounterCoverage])
            
            # go through each grid point 
            # print(coverage[idx,timeCounterCoverage])

            
           
            
 
            
            timeCounterCoverage = timeCounterCoverage+1 

        ################ Phasing multiple satellites  ########################
        
        satellites = numSats[counter]
        
        satIDXstep = round(loops/satellites)
        for a in range(satellites):
            for b in range(totalPoints): 
                point = coverage[b,:]
                point = reorder(list(point),a*satIDXstep)
                for c in range(loops): 
                    if point[c] == 1: 
                        coverageMain[b,c] = 1
       
        if np.count_nonzero(coverageMain) == loops*totalPoints: 
            print('working')
            break          
                
    counter = counter + 1
         


#  def reorder(lst, newIdx):
#     pos = newIdx #lst.index(first)
#     return lst[pos:] + lst[:pos]       
    
print('Calculating FOM...')
###################################################
############## calculate FOM ######################     done once per design
###################################################

coverageMainNP = np.array(coverageMain)
if np.count_nonzero(coverageMainNP) == loops * totalPoints:
    percentCoverage = np.zeros(totalPoints) + 100  # % covered
    maxCoverageGap = np.zeros(totalPoints)  # longest coverage gap by each point
    meanCoverageGap = np.zeros(totalPoints)  # average coverage gap for each point
    timeAvgGap = np.zeros(totalPoints)  # time avg gap for each point
    meanResponseTime = np.zeros(totalPoints)
else:
    percentCoverage = np.zeros(totalPoints)  # % covered
    maxCoverageGap = np.zeros(totalPoints)  # longest coverage gap by each point
    meanCoverageGap = np.zeros(totalPoints)  # average coverage gap for each point
    timeAvgGap = np.zeros(totalPoints)  # time avg gap for each point
    meanResponseTime = np.zeros(totalPoints)
    for i in range(totalPoints):
        arr = coverageMain[i, :]
        checkArr = arr
        checkArr = np.append(checkArr,1)  # adding a 1 to the start and end so the diff function works if theres no LOS to start
        checkArr = np.insert(checkArr, 0, 1)
        # find the indices where the array changes value from 0 to 1
        zero_indices = np.where(checkArr == 1)[0]  # Find the indices of all zeros
        zero_diffs = np.diff(zero_indices)  # Compute the differences between adjacent indices

        maxCoverageGap[i] = (max(np.concatenate(([0], zero_diffs)) - 1))*100  # Find the maximum number of repeated zeros
        print("MCG ",maxCoverageGap[i])
        #print("Max cov gap ", maxCoverageGap)
        zero_gaps = np.concatenate(([0], zero_diffs)) - 1

        ### Getting Mean Coverage Gap ########
        zero_list = np.delete(zero_gaps, [0])  # deleting the first appended element because it always comes to -1 due to line 28, i wont lose any data from this
        # print("Zeroslist ", zero_list)
        numGaps = np.count_nonzero(zero_list)
        meanCoverageGap[i] = (np.sum(zero_list) / numGaps)*100
        #print("mean cog gap ", meanCoverageGap)

        ######### GETTING PERCENT COVERAGE ############
        percentCoverage[i] = (np.count_nonzero(arr) / (len(arr))) * 100
        #print("percent coverage ", percentCoverage)
        ######## TIME AVG GAP #########
        timeAvgGap[i] = ((np.sum(zero_list ** 2)) / len(arr)) * 100
        #print("Time Avg Gap", timeAvgGap)

        ####### MEAN RESPONSE TIME #######
        # using formula for nth trinagular numbers
        # zero array * ((n + (n+1) / 2)
        mask = zero_list != 0
        zero_list_MRT = zero_list[mask]
        MRTvecA = zero_list_MRT + 1  # this is the (n + 1) step
        MRTvecB = MRTvecA * zero_list_MRT  # this is the (n * (n+1)) step
        MRTvecC = MRTvecB / 2  # final MRT vec # this the (n * (n+1)) /2 step

        meanResponseTime[i] = (sum(MRTvecC) / len(arr))*100
        #print("MRT ", meanResponseTime)

print("Length PERC COV ",len(percentCoverage),"Length MAX GAP ",len(maxCoverageGap),"Length MEAN GAP ",len(meanCoverageGap),"MRT ",len(meanCoverageGap), "Length TIME AVG GAP ",len(timeAvgGap))
###############################################
#################### Score / Return ####################
###############################################

# return

# weightPC = 0.2
# weightMaxCG = 0.15
# weightMeanCG = 0.1
# weightTAR = 0.25
# weightMRT = 0.3
# score = sum(percentCoverage)*weightPC + (1/sum(meanCoverageGap))*weightMeanCG + (1/sum(timeAvgGap))*weightTAR + (1/sum(maxCoverageGap))*weightMaxCG
# + (1/sum(meanResponseTime))*weightMRT

###############################################
#################### EXTRA ####################
###############################################

# ax.plot(orb.x,orb.y,orb.z)
# ax.set_aspect('auto')
for i in range(totalPoints):
     print('Point:', i, ', Percent Coverage:', round(percentCoverage[i]), '%, Max Coverage Gap:',
           round(maxCoverageGap[i]) * sToHr, 'hr, Mean Coverage Gap:', round(meanCoverageGap[i]) * sToHr,
           'hr, Time avg Gap:', round(timeAvgGap[i]) * sToHr, 'hr, Mean Response Time:',
           round(meanResponseTime[i]) * sToHr, 'hr')

############################
avgPercentCoverage = np.mean(percentCoverage)
avgMaxCoverageGap = np.mean(maxCoverageGap)
avgMeanCoverageGap = np.mean(meanCoverageGap)
avgTimeAvgGap = np.mean(timeAvgGap)
avgMeanResponseTime = np.mean(meanResponseTime)

print("PERC COV ",avgPercentCoverage,"MAX GAP ",avgMaxCoverageGap,"MEAN GAP ",avgMeanCoverageGap,"MRT ",avgMeanResponseTime, "TIME AVG GAP ",avgTimeAvgGap)


# # ax.scatter(coordinates[:][])
# # time = (coverage.count(1))*100
# # print(time, 'in Seconds')
end = time.time()
print('run time: ', end-start)

# plt.show()

