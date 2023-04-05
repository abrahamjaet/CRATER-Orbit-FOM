import matplotlib.pyplot as plt
import numpy as np
import orbitDict
import time 
from csltk.utilities import System
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
            # sat_xPos = orb.x
            # sat_yPos = orb.y
            # sat_zPos = orb.z

            r_spacecraft = [sat_xPos, sat_yPos, sat_zPos] # vector from center of earth to satellite
            # print(r_spacecraft)
            # time.sleep(20)
            r_spacecraftToPoint = [grid_point_coordinates[:,0] - sat_xPos, grid_point_coordinates[:,1] - sat_yPos, grid_point_coordinates[:,2] - sat_zPos] # vector sat to point
            r_moonToPoint = [grid_point_coordinates[:,0]-moonPos[0],grid_point_coordinates[:,1]-moonPos[1],grid_point_coordinates[:,2]-moonPos[2]]# vector from center of moon to point 
            r_moonToSat = [r_spacecraft[0]-moonPos[0],r_spacecraft[1]-moonPos[1],r_spacecraft[2]-moonPos[2]] # vector from moon center to sat 
            # print(r_moonToSat)
            r_moonToSat = np.array(r_moonToSat)
            # print(r_moonToSat)
            r_moonToSat = r_moonToSat.reshape(3,1)
            # print(r_moonToSat)
            # time.sleep(10)
            print(r_spacecraft)
            print(len(r_spacecraft))
            print(len(grid_point_coordinates[:,0]))
            print(r_spacecraftToPoint)
            print(len(r_spacecraftToPoint))

            time.sleep(20)
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
            print(angle1)
            # print(idxAngle1)
            coverage[idx,timeCounterCoverage] = 1

            # go through each grid point 
            # print(coverage[idx,timeCounterCoverage])
            time.sleep(1)
            
 
            
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
if np.count_nonzero(coverageMain) == loops*totalPoints:
    percentCoverage = np.zeros(totalPoints) + 100 # % covered 
    maxCoverageGap = np.zeros(totalPoints) # longest coverage gap by each point 
    meanCoverageGap = np.zeros(totalPoints) # average coverage gap for each point 
    timeAvgGap = np.zeros(totalPoints) # time avg gap for each point 
    meanResponseTime = np.zeros(totalPoints) 
else: 
    percentCoverage = np.zeros(totalPoints) # % covered 
    maxCoverageGap = np.zeros(totalPoints) # longest coverage gap by each point 
    meanCoverageGap = np.zeros(totalPoints) # average coverage gap for each point 
    timeAvgGap = np.zeros(totalPoints) # time avg gap for each point 
    meanResponseTime = np.zeros(totalPoints) 
    for i in range(totalPoints):

        percentCoverage[i] = (np.count_nonzero(coverageMain[i,:])/loops)*100 # percent of time each grid point is seen [%]
        
        if np.count_nonzero(coverageMain[i,:]) == loops:
            #print('all 1')
            timeAvgGap[i] = 0 
            meanCoverageGap[i] = 0 
            maxCoverageGap[i] = 0
            meanResponseTime[i] = 0 
        elif np.count_nonzero(coverageMain[i,:]) == 0:
            #print('all 0')
            meanCoverageGap[i] = loops*100 # mean coverage gap for grid point i [s]
            timeAvgGap[i] = loops*100 # [s]
            maxCoverageGap[i] = loops*100 # seconds of maximum coverage gap for grid point i [s]
            meanResponseTime[i] = loops*100 # s
        else:
            #print('mix')
            counterMCG = 0 # counter consecutive 0s, coverage gap 
            coverageGaps = [] # length of coverage gaps  
            meanRT = []
            for j in range(loops):
                # if coverage[i,j] == 1: 
                #     print(i,coverage[i,j])
                if coverageMain[i,j] == 0: # if uncovered 
                    counterMCG = counterMCG + 1 # +1 counter of gap 
                    meanRT.append(counterMCG)  
                else: 
                    if not counterMCG == 0:
                        coverageGaps.append(counterMCG)
                        counterMCG = 0  
            if not counterMCG == 0:
                coverageGaps.append(counterMCG)               
            numGaps = len(coverageGaps)
            if numGaps == 0: 
                numGaps = 1
            meanCoverageGap[i] = (sum(coverageGaps)/numGaps)*100 # mean coverage gap for grid point i [s]
            timeAvgGap[i] = (sum(np.array(coverageGaps)**2)/loops)*100 # [s]
            maxCoverageGap[i] = max(coverageGaps)*100 # seconds of maximum coverage gap for grid point i [s]
            meanResponseTime[i] = (sum(meanRT)/loops)*100 # s
            
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
    print('Point:',i,', Percent Coverage:',round(percentCoverage[i]),'%, Max Coverage Gap:',round(maxCoverageGap[i])*sToHr,'hr, Mean Coverage Gap:',round(meanCoverageGap[i])*sToHr,'hr, Time avg Gap:',round(timeAvgGap[i])*sToHr,'hr, Mean Response Time:',round(meanResponseTime[i])*sToHr,'hr')
    #print('Percent Coverage:',round(percentCoverage[i]),'%, Max Coverage Gap:',round(maxCoverageGap[i]*sToHr),'hr, Mean Coverage Gap:',round(meanCoverageGap[i]*sToHr),'hr, Time avg Gap:',round(timeAvgGap[i]*sToHr),'hr, Mean Response Time:',round(meanResponseTime[i]*sToHr),'hr')



############################


# # ax.scatter(coordinates[:][])
# # time = (coverage.count(1))*100
# # print(time, 'in Seconds')
# end = time.time()
# print('run time:',(end-start)/60,'min')

# plt.show()

