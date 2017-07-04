__author__ = 'HOI0920H'
import itertools
import sys
import re
input=sys.argv[1]

def find_in_file(file,columns):
    filepointer=open(file,"r")
    # print("looking for columns:")
    # print(columns)
    number_of_patient_samples=0
    for line in filepointer:
        # line.strip()
        # print (line)
        if re.match(r'^sample_id',line):
            continue
        line_list=line.strip().split(",")
        # print ("line split into %s" %line_list)

        found_zero="false"
        for ele in columns:
            # print ("line list value %s" %line_list[ele])
            if line_list[ele] == "0":
                found_zero="true"
                # print("found a zero!! breaking!")
                break

        if found_zero != "true":
            number_of_patient_samples+=1

    if number_of_patient_samples >0:
        # print ("%s samples falls under columns %s" %(number_of_patient_samples , convert_to_header(columns)))
        #string this up for R
        stringfy_intersection= "%s=%s" %(convert_to_header(columns),number_of_patient_samples)
        print(stringfy_intersection)
    return stringfy_intersection

def convert_to_header(array):
    options = {
           1 : "SC",
           2 : "VF",
           3 : "SM",
           4 : "RNA",
           5 : "Blood",
           6 : "Plasma",
           7 : "Serum",
           8 : "Urine(Short)",
           9 : "Urine(Tall)",
           10: "Urine_present"
           }
    VennDiagram={
        1 : "1",
        2 : "2",
        3 : "3",
        4 : "4",
    }
    string_array=[]
    for item in array:
        # print (item)
        # string_array.append(options[item])
        # string_array.append(str(item))
        string_array.append(VennDiagram[item])
    combine_header= "".join(string_array)
    print(combine_header)
    return combine_header


all_overlaps=[]
outfile=open(sys.argv[1]+".overlap","w")
stuff = list(range(1,5))
counter=0
max_size=len(stuff)+2
for L in range(1, max_size):
    for subset in itertools.combinations(stuff, L):
        # print(subset)
        #samples will be returned a string of the overlap
        samples = find_in_file(input, subset)
        all_overlaps.append(samples)

    # counter+=1

# print (counter)
final_string=",".join(all_overlaps)
print (final_string)
outfile.write(final_string)
outfile.close()




