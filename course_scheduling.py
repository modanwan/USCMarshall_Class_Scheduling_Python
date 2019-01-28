import pandas as pd
import numpy as np
import gurobipy as grb

mod = grb.Model()


#%% Prepare data
data = pd.read_excel("570_example_data.xlsx",sheet_name=None,index_col=0)

courseL = list(data['Course_Prop'].index)
roomL = list(data['Room_Prop'].index)
slotL = list(data['SlotID_Time'].index)

course2L = []
for course in courseL:
    if (data['Course_Prop'].loc[course,'Dur_Slots']==2):
        course2L.append(course)

timeL = sorted(list(set(data['SlotID_Time'].loc[:,'Time'])))
wdayL = ['Mon','Tue','Wed','Thu','Fri']

first_slot = slotL[0]
last_slot = slotL[len(slotL)-1]

n_slot = len(timeL)
n_wday = len(wdayL)

ttl_course = len(data['Course_Prop']) 
ttl_imp    = sum(data['Course_Prop'].loc[:,'Importance'])

ttl_course_time_pref = sum(np.max(data['Course_Time_Pref'].values, axis=1))
ttl_prof_time_pref   = 0
for  course in courseL:
    ttl_prof_time_pref += max(data['Prof_Time_Pref'].loc[data['Course_Prop'].loc[course,'Professor']])

pair_slotL = [0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15]
cons_slotL = [5, 7, 14, 16, 23, 25, 32, 34, 41, 43]
        
count=0
#%% Create decision variables
# Main DV
x={}
for course in courseL:
    for room in roomL:
        for slot in slotL:
            x[course,room,slot] = mod.addVar(lb=0, vtype=grb.GRB.BINARY)
         
# Auxilary DV to set constraints to (0,2) for 2-slot courses
x_2slot_c={}
x_same_room={}
x_ith_pair_slot={}
x_ith_cons_slot={}
for course in courseL:
    if(data['Course_Prop'].loc[course,'Dur_Slots']==2):
        x_2slot_c[course] = mod.addVar(lb=0, vtype=grb.GRB.BINARY)
        for room in roomL:
            x_same_room[course,room] = mod.addVar(lb=0, vtype=grb.GRB.BINARY)
            if(data['Course_Prop'].loc[course,'Pair']):
                for i in pair_slotL:
                    x_ith_pair_slot[course,room,i] = mod.addVar(lb=0, vtype=grb.GRB.BINARY)
            else:
                for i in cons_slotL:
                    x_ith_cons_slot[course,room,i] = mod.addVar(lb=0, vtype=grb.GRB.BINARY)

# Upper bound to spread scheduling across rooms
U_room = mod.addVar(lb=0, vtype=grb.GRB.INTEGER)

# Upper bound to spread scheduling across week days
U_wday = mod.addVar(lb=0, vtype=grb.GRB.INTEGER)

                    
#%% Create basic constraints
# Course slots = 1 or 2
for course in courseL:
    # 1-slot course:
    if(data['Course_Prop'].loc[course,'Dur_Slots']==1):
        # Take 1/0 roomroom-slot
        mod.addConstr(sum(sum(x[course,room,slot] for slot in slotL)for room in roomL)<=1)
    # 2=slot course: 
    elif(data['Course_Prop'].loc[course,'Dur_Slots']==2):
        # Take 2/0 room-slot
        mod.addConstr(sum(sum(x[course,room,slot] for slot in slotL)for room in roomL)==2*x_2slot_c[course])
        # NO same time; diff rooms
        for slot in slotL:
            mod.addConstr(sum(x[course,room,slot] for room in roomL)<=1)
        # One course should be allocated in one room
        for room in roomL:
            mod.addConstr(sum(x[course,room,slot] for slot in slotL)==2*x_same_room[course,room])
        # 2-slot course scheduling: pair or consecutive
        if(data['Course_Prop'].loc[course,'Pair']): # pair
            for room in roomL:
                for i in pair_slotL:
                    mod.addConstr(x[course,room,'slot_'+str(i)]+x[course,room,'slot_'+str(i+n_slot*2)]==2*x_ith_pair_slot[course,room,i])
                # block Friday slots, cannot be paired
                mod.addConstr(sum(x[course,room,'slot_'+str(j)] for j in range(n_slot*4,n_slot*5))==0)
                mod.addConstr(sum(x[course,room,'slot_'+str(d*n_slot+7)]+x[course,room,'slot_'+str(d*n_slot+8)] for d in range(n_wday-1))==0)
                
        else: # consecutive; no overnight
            for room in roomL:
                # pick one of the consecutive slots (daily last 4 slots)
                for i in cons_slotL:
                    mod.addConstr(x[course,room,'slot_'+str(i)]+x[course,room,'slot_'+str(i+1)]==2*x_ith_cons_slot[course,room,i])
                # block non-consecutive slots (daily first 5 slots)
                mod.addConstr(sum(sum(x[course,room,'slot_'+str(w*n_slot+j)] for j in range(5))for w in range(n_wday))==0)
                
# one course at same room same time
for room in roomL:
    for slot in slotL:
        mod.addConstr(sum(x[course,room,slot] for course in courseL)<=1)
        

#%% Create Room-Course constraints
# Room vs. Course conflicts (1=Yes, 0=No)
course_Room_conflictM = pd.DataFrame(np.zeros((len(courseL),len(roomL))),index=courseL,columns=roomL)

# Room capacity
for course in courseL:
    for room in roomL:
        if(data['Course_Prop'].loc[course,'Enrollment']>data['Room_Prop'].loc[room,'Capacity']):
            course_Room_conflictM.loc[course,room] += 1

# Room functionality
for course in courseL:
    for room in roomL:
        if(np.product(list(data['Course_Prop'].loc[course,'F_Computer':'F_Disabled'] <= data['Room_Prop'].loc[room,'F_Computer':'F_Disabled']))==0):
            course_Room_conflictM.loc[course,room] += 10
            
# Convert conflit matrix into constraints
for course in courseL:
    for room in roomL:
        if(course_Room_conflictM.loc[course,room]>0):
            mod.addConstr(sum(x[course,room,slot] for slot in slotL)==0)


#%% Create Slot-Course constraints
# Slot vs. Course conflicts (1=Yes, 0=No)
course_Slot_conflictM = pd.DataFrame(np.zeros((len(courseL),len(slotL))),index=courseL,columns=slotL)

# Professor unavailability (Time)
for course in courseL:
    course_Slot_conflictM.loc[course,:] += 1-data['Prof_Slot_Avail'].loc[data['Course_Prop'].loc[course,'Professor'],first_slot:last_slot]

# Convert conflit matrix into constraints
for course in courseL:
    for slot in slotL:
        if(course_Slot_conflictM.loc[course,slot]>0):
            mod.addConstr(sum(x[course,room,slot] for room in roomL)==0)
            

#%% Create Room-Slot constraints
# Slot vs. Room conflicts (1=Yes, 0=No)
room_Slot_conflictM = pd.DataFrame(np.zeros((len(roomL),len(slotL))),index=roomL,columns=slotL)

# Room unavailability (Time)
room_Slot_conflictM = 1-data['Room_Slot_Avail'].loc[:,first_slot:last_slot]

# Convert conflit matrix into constraints
for room in roomL:
    for slot in slotL:
        if(room_Slot_conflictM.loc[room,slot]>0):
            mod.addConstr(sum(x[course,room,slot] for course in courseL)==0)
            

#%% Create Paired courses constraints
# List of paired courses with time conflicts
paired_Conflict_CourseL = data['Conflict_Pair']

# Adding courses taught by one professor
for i in range(len(courseL)-1):
    for j in range(1+i,len(courseL)):
        if(data['Course_Prop'].loc[courseL[i],'Professor']==data['Course_Prop'].loc[courseL[j],'Professor']):
            paired_Conflict_CourseL = paired_Conflict_CourseL.append({'Course_ID1':courseL[i], 'Course_ID2':courseL[j]}, ignore_index=True)

# Convert conflit list into constraints
for slot in slotL:
    for i in range(len(paired_Conflict_CourseL)):
        course_i,course_j = paired_Conflict_CourseL.loc[i]
        mod.addConstr(sum(x[course_i,room,slot]+x[course_j,room,slot] for room in roomL)<=1)
            

#%% Spread the scheduling across time slots and rooms
for room in roomL:
    mod.addConstr(sum(sum(x[course,room,slot] for slot in slotL) for course in courseL)<=U_room)

for d in range(n_wday):
    mod.addConstr(sum(sum(sum(x[course,room,'slot_'+str(d*n_slot+i)] for i in range(n_slot)) for room in roomL) for course in courseL)<=U_wday)


#%% Set up objective
# Course scheduling completion (scheduled units/total units) <MAX>
obj_schedule_completion = sum(sum(sum(x[course,room,slot]/data['Course_Prop'].loc[course,'Dur_Slots'] for slot in slotL) for room in roomL) for course in courseL)/ttl_course

# Total importance scheduled <MAX>
obj_course_importance = sum(sum(sum(x[course,room,slot]/data['Course_Prop'].loc[course,'Dur_Slots']*data['Course_Prop'].loc[course,'Importance'] for slot in slotL) for room in roomL) for course in courseL)/ttl_imp

# Course time preference
obj_course_time_pref = 0
for course in courseL:
    course_time_pref = data['Course_Time_Pref'].loc[course]
    if(sum(course_time_pref)):
        course_time_prefL = pd.Series([course_time_pref[0],
                                       course_time_pref[1],course_time_pref[1],
                                       course_time_pref[2],course_time_pref[2],course_time_pref[2],course_time_pref[2],
                                       course_time_pref[3],course_time_pref[3]]*5,
                                       index = slotL)
        obj_course_time_pref += sum(sum(x[course,room,slot]/data['Course_Prop'].loc[course,'Dur_Slots']*course_time_prefL.loc[slot] for slot in slotL) for room in roomL)

obj_course_time_pref /= ttl_course_time_pref

# Professor time preference
obj_prof_time_pref = 0
for course in courseL:
    prof_time_pref = data['Prof_Time_Pref'].loc[data['Course_Prop'].loc[course,'Professor']]
    if(sum(prof_time_pref)):
        prof_time_prefL = pd.Series([prof_time_pref[0],
                                     prof_time_pref[1],prof_time_pref[1],
                                     prof_time_pref[2],prof_time_pref[2],prof_time_pref[2],prof_time_pref[2],
                                     prof_time_pref[3],prof_time_pref[3]]*5,
                                     index = slotL)
        obj_prof_time_pref += sum(sum(x[course,room,slot]/data['Course_Prop'].loc[course,'Dur_Slots']*prof_time_prefL.loc[slot] for slot in slotL) for room in roomL)

obj_prof_time_pref /= ttl_prof_time_pref
        
# Entropy <MIN>
obj_schedule_balance = -(U_room+U_wday)

# Final weighted objective
obj = 10000 * obj_schedule_completion \
    + 1000  * obj_course_importance \
    + 100   * obj_course_time_pref \
    + 10    * obj_prof_time_pref \
    + 0.001 * obj_schedule_balance

mod.setObjective(obj, sense=grb.GRB.MAXIMIZE)
            
#%% Run optimization & output           
mod.setParam('TimeLimit',60)
mod.optimize()

# Output optimized solution
opt_x={}
opt_cour={}
for course in courseL:
    opt_cour[course]=[]
    for room in roomL:
        for slot in slotL:
            if (x[course,room,slot].x):
                opt_x[room,slot] = course
                opt_cour[course].append(str(room))
                opt_cour[course].append(str(data['SlotID_Time'].loc[slot,'Wday']+" "+data['SlotID_Time'].loc[slot,'Time']))

table = pd.DataFrame(index=roomL, columns=slotL)
for key,value in opt_x.items():
    table.loc[key]=value
    
table.to_excel('output_schedule.xlsx')

#%% Verifications

## Total scheduled courses
#sum(sum(sum(x[course,room,slot].x/data['Course_Prop'].loc[course,'Dur_Slots'] for slot in slotL) for room in roomL) for course in courseL)
#ttl_course
#
## Total scheduled importances
#sum(sum(sum(x[course,room,slot].x/data['Course_Prop'].loc[course,'Dur_Slots']*data['Course_Prop'].loc[course,'Importance'] for slot in slotL) for room in roomL) for course in courseL)
#ttl_imp
#
## Weekday UBound
#for d in range(n_wday):
#    print(sum(sum(sum(x[course,room,'slot_'+str(d*n_slot+i)].x for i in range(n_slot)) for room in roomL) for course in courseL))
#U_wday.x
#
## Room UBound
#for room in roomL:
#    print(sum(sum(x[course,room,slot].x for slot in slotL) for course in courseL))
#U_room.x
#
## Course time pref
#opt_course_time_pref = 0
#for course in courseL:
#    course_time_pref = data['Course_Time_Pref'].loc[course]
#    if(sum(course_time_pref)):
#        course_time_prefL = pd.Series([course_time_pref[0],
#                                       course_time_pref[1],course_time_pref[1],
#                                       course_time_pref[2],course_time_pref[2],course_time_pref[2],course_time_pref[2],
#                                       course_time_pref[3],course_time_pref[3]]*5,
#                                       index = slotL)
#        opt_course_time_pref += sum(sum(x[course,room,slot].x/data['Course_Prop'].loc[course,'Dur_Slots']*course_time_prefL.loc[slot] for slot in slotL) for room in roomL)
#ttl_course_time_pref
#
## Professor time preference
#opt_prof_time_pref = 0
#for course in courseL:
#    prof_time_pref = data['Prof_Time_Pref'].loc[data['Course_Prop'].loc[course,'Professor']]
#    if(sum(prof_time_pref)):
#        prof_time_prefL = pd.Series([prof_time_pref[0],
#                                     prof_time_pref[1],prof_time_pref[1],
#                                     prof_time_pref[2],prof_time_pref[2],prof_time_pref[2],prof_time_pref[2],
#                                     prof_time_pref[3],prof_time_pref[3]]*5,
#                                     index = slotL)
#        count+=1
#        opt_prof_time_pref += sum(sum(x[course,room,slot].x/data['Course_Prop'].loc[course,'Dur_Slots']*prof_time_prefL.loc[slot] for slot in slotL) for room in roomL)
#ttl_prof_time_pref


