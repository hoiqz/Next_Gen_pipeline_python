__author__ = 'HOI0920H'

# run this code to simulate a BCL folder structure
# run this code to simulate a fastq folder structure

#import string and random for use in generating the random name
import string
import random

import os
import datetime
time=datetime.date.today()
yymmdd=time.strftime("%y%m%d")
print(yymmdd)
# print (time.strftime("%y%m%d"))

# randomletter_size=4
def random_gen(size=4,char=string.ascii_uppercase+string.digits):
    return "".join(random.choice(char) for _ in range(size))

#test the random gen
print (random_gen())
print (random_gen(10))
uid=random_gen()
exp_name=yymmdd+"_CRUNextSeq_" + uid+"_FC"
print(exp_name)
currentwd=os.getcwd()
print(currentwd+"\\"+ exp_name)
base= currentwd+ '\\' + exp_name
intensities_bylane=base+ "\\Data\\Intensities"+"\\"+"L00"

if not os.path.exists(intensities_bylane):
    for num in range(1,10):
        print (intensities_bylane+str(num))
        os.makedirs(intensities_bylane+str(num))


intensities_bybasecall=base+ "\\Data\\Intensities"+"\\"+"BaseCalls"
intensities_bybasecall_bylane=base+ "\\Data\\Intensities"+"\\"+"BaseCalls"+"\\"+"L00"

if not os.path.exists(intensities_bybasecall_bylane):
    for num in range(1,10):
        print ("creating directories: %s%d" %(intensities_bybasecall_bylane,num))
        os.makedirs((intensities_bybasecall_bylane+str(num)))
        os.makedirs(intensities_bybasecall_bylane+str(num)+"\\C1.1")
        open(intensities_bybasecall_bylane+str(num)+"\\C1.1\\test.bcl","a").close()
        #make dummy fastq files
        open(intensities_bybasecall+"\\sample"+str(num)+".fastq","w").close()







