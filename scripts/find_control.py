import re
import sys
infile=sys.argv[1]
infile2=sys.argv[2]
case=open(infile,"r")
case_control=open("case_control_pair.tsv","w")
alreadychose=[]

def ethnicity_code(eth):
    if eth =='1':
        return 'Chinese'
    elif eth =='2':
        return 'Malay'
    elif eth =='3':
        return 'Indian'
    else:
        return 'Others'

def gender_code(gen):
    if gen=='1':
        return 'Female'
    elif gen=='2':
        return 'Male'
    else:
        return 'Unknown'

for line in case:
    found_3=0
    match=re.search(r"^Metformin",line)
    if match:
        header=line
        print (header)
        continue
    patientinfo=line.split("\t")
    case_id=patientinfo[1]
    case_age=int(patientinfo[7])
    case_gender=patientinfo[9]
    case_eth=patientinfo[8]
#    print ("%s %s %s" %(case_age,case_gender,case_eth))
    
    control=open(infile2,"r")
    for con_line in control:
        control_match=re.search(r"^study_id",con_line)
        if control_match:
            continue
        control_info=con_line.split("\t")
        con_id=control_info[0]
        if control_info[9] =='':
            continue
        con_age=int(control_info[9])

        con_eth=control_info[10]
        con_eth=ethnicity_code(con_eth)
        con_gender=control_info[11]
        con_gender=gender_code(con_gender)

#        print ("controls :%s %s %s" %(con_age,con_gender,con_eth))
   #     if case_id != con_id and case_gender == con_gender: 
        if case_id != con_id and case_gender == con_gender and con_age<=case_age+5 and con_age >= case_age-5 and con_eth ==case_eth and con_id not in alreadychose:
            print ("%s %s %s %s match %s %s %s %s"%(case_id,case_age,case_gender,case_eth,con_id,con_age,con_gender,con_eth))
            alreadychose.append(con_id)
            case_control.write(con_line)
#            case_control.write("%s %s %s %s match %s %s %s %s \n"%(case_id,case_age,case_gender,case_eth,con_id,con_age,con_gender,con_eth))

            found_3+=1
        if found_3==3:
            break


    control.close


