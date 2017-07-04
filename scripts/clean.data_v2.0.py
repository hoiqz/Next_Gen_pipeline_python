import datetime
import re
import sys

def getTimeDifferenceFromNow(start,end):
    timeDiff=start-end
    return timeDiff.total_seconds() / 60

########## iterate the file first time to get the <4% baseline value in 30mins##
##############window####
###this version checks the next 10 seconds before declaring it as an episode
file = sys.argv[1]
#print(file)
filereader=open(file)
filewriter=open(file+".out","w")
toskip= 10*60 # delete the first 10 mins of the readings
skip=0 # skip the first two lines becuase i assume it is the machine extra info
baselinecounter=0
baseline_in_sec=30*60 # time in secs to calculate baseline value
totalOSreading=0
OSreadingCounter=0
equipmentStartTime=0
episodescounter=0
for line in filereader:
    if skip<toskip:
        skip=skip+1
        next
    elif baselinecounter < baseline_in_sec:
        baselinecounter =baselinecounter +1
        line=line.replace("\n","")
        line=line.replace("\"","")
        line=line.replace("\t","")
        #line=line.replace(" ","")
        line=re.sub("\s+","",line)
#        print(line)
        linearray=line.split(",")
        unformatted_date=linearray[0]+linearray[1]
#        print (linearray[2])
        if linearray[2]=='' or float(linearray[2]) < 50:
#            print ("weird OS reading %s , skipping"%linearray[2])
            next
        else:
            totalOSreading=totalOSreading+float(linearray[2])
            OSreadingCounter=OSreadingCounter+1
    else :
        break

print ("total OS reading %.2f"%totalOSreading)
print ("total counter %d"%OSreadingCounter)
meanOSreading=float(totalOSreading/OSreadingCounter)
fourpercentOS=0.96*meanOSreading
filewriter.write("mean reading %.2f\n"%meanOSreading)
filewriter.write("low reading %.2f\n"%fourpercentOS)
filereader.close()

##### iterate the file again now doing the main counting ####
filereader=open(file)
skip=0 # skip the first two lines becuase i assume it is the machine extra info
findlowreading=1
findbaselinereading=0
start_session="0"
start_epi_reading=0.00
episode_start=0 # 0 = no ,1 = yes
ten_second_checker=1 # this checker is to count the number of seconds below
#baseline

for line in filereader:
    if skip<toskip:
        skip=skip+1
        next
        # note here that we do not need to do the 30 mins baseline cal again.
        # justproceed
    else:
        line=line.replace("\n","")
        line=line.replace("\"","")
        line=line.replace("\t","")
        line=line.replace(" ","")
        print(line)
        linearray=line.split(",")
        unformatted_date=linearray[0]+linearray[1]
        print(unformatted_date)
        if linearray[2]=='' or float(linearray[2]) < 50:
            next
        else:
            tmpstart=datetime.datetime.strptime(unformatted_date, '%m/%d/%y%H:%M:%S')
            print (tmpstart)
            OSreading=float(linearray[2])
    #        print (OSreading)
            if equipmentStartTime == 0:
                equipmentStartTime=tmpstart
                filewriter.write("measurement start time: %s\n"%equipmentStartTime)
            if findlowreading==1 and OSreading <= fourpercentOS:
                if start_session == '0':
                    start_session =str(tmpstart)
                    start_epi_reading=OSreading
                if episode_start ==1 and ten_second_checker ==10 and start_session:
                    #means that episode has lasted 10 seconds
                    episodescounter=episodescounter+1
                    #print ("Episode at ",tmpstart," with reading %d"%OSreading)
                    filewriter.write("Episode started at %s with reading %.1f below 4percent baseline of %.1f\n"%(start_session,start_epi_reading,fourpercentOS))
                    #filewriter.write("Episode started at %s with reading %.1f below 4percent baseline of %.1f\n"%(str(tmpstart),OSreading,fourpercentOS))
                    findbaselinereading=1
                    findlowreading=0
                    ten_second_checker=1 # reset the episde checker to 1
                else:
                    # if anyy thing else, jus add a sec to episode time
                    episode_start =1
                    ten_second_checker = ten_second_checker +1
            elif findlowreading==1 and OSreading > fourpercentOS:
                if episode_start==1:
                    episode_start=0
                    ten_second_checker=1 # reset the episde checker to 1
                    start_session ='0'

            elif findbaselinereading==1 and OSreading >=meanOSreading:
    #            print ("OS reading return to baseline at ",tmpstart," with reading %d"%OSreading)
                filewriter.write("OS reading return to baseline at %s with reading %.1f\n"%(str(tmpstart),OSreading))
                findlowreading=1
                findbaselinereading=0
                start_session='0'
            else:
                next
filewriter.write("Total Episodes = %d\n"%episodescounter)

filereader.close()
filewriter.close()

   
