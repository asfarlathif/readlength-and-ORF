"""TO READ A MULTIFASTA FILE AND ANALYZE THE RECORDS AND DNA SEQUENCES IN IT"""

def maxmin(record_dict): 
    """TO IDENTIFY THE MAXIMUM AND MINIMUM LENGTH SEQUENCES IN THE FILE AND RECORD ITS LENGTHS"""
    max_min = input("max OR min:") 
#Maximum or Minimum sequence length choice
    if max_min== 'max':
        record=(max(record_dict,key=record_dict.get)) #getting the ID of the first sequence with maximum length
    elif max_min== 'min':
        record=(min(record_dict,key=record_dict.get)) #getting the ID of the first sequence with minimum length
    else:
        print("INVALID INPUT")
        return
    seq_list=[record] #list of IDs of the sequences with maximum/minimum length
    no_of_seq=1
#checking whether there are more sequneces with same maximum/minimum length
    for a in record_dict: 
        if record_dict[a]!=record_dict[record]:
            continue
        else:
            if a==record:
                continue
            else:
                no_of_seq+=1
                seq_list.append(a)
    list_dict={}
    for a in seq_list: #creating dictionary of the final list of sequences with maximum/minimum length
        list_dict[a]=record_dict[a]
#printing the results
    if no_of_seq==1:
        if max_min=='max':
            print("There is one sequence in the record with maximum length")
            print("Longest sequence in the file is:")
            print("%s and its length is" % record, record_dict[record])
        else:
            print("There is one sequence in the record with minimum length")
            print("Shortest sequence in the file is:")
            print("%s and its length is" % record, record_dict[record])
    else:
        if max_min=='max':
            print("There are %d sequences in the record with maximum length" % len(list_dict))
            print("Longest sequences in the file are:")
            for a in list_dict:
                print("%s and its length is" % a, record_dict[a])
        else:
            print("There are %d sequences in the record with minimum length" % len(list_dict))
            print("Shortest sequences in the file are:")
            for a in list_dict:
                print("%s and its length is" % a, record_dict[a])
    x=input("Do you want to check again?") #to call the function again
    if x =='yes':
        maxmin(record_dict)
    else:
        return    
def startstop_codon(dna,frame):
    """ TO CHECK THE PRESENCE OF START CODON AND ITS POSITION"""
    start_codons=['ATG','atg'] ; stop_codons=['TAA','TAG','TGA','taa','tag','tga'] #referece list of start and stop codons
    for i in range(frame,len(dna),3):
        codon1=dna[i:i+3]
        if codon1 in start_codons:
            position1=i#getting the index of the start codon
            ORF_start=position1+1
            for j in range(position1,len(dna),3):
                codon2=dna[j:j+3]
                if codon2 in stop_codons:
                    position2=j#getting the index of the stop codon
                    ORF_stop=position2+1
                    break #terminating the loop when a stop codon is found
            break
    try:
        return [len(dna[position1:position2])+3,ORF_start]
    except UnboundLocalError:
        None 

def finding_ORF(frame,records):
    """TO IDENTIFY THE ORF IN THE READING FRAMES AND 
    IDENTIFY THE LONGEST ORF IN EACH SEQUENCE AND THE POSITION OF THE START CODON"""
    sequences=[reads.seq for reads in records] #creating alist of all the sequence reads as elements of the list
    ORF_lengths={}
    positions={}
    max_lenth_ORFs=[]
    for read in range(len(records)): #loop for accessing each sequence read from the list
        seq=str(sequences[read])
        a=startstop_codon(seq,frame) #calling the codon reading function
        if a== None: #assigning value for reads with no codons
            a = [0,0]
        ORF_lengths[records[read].id]=a #creating a dictionary of ORF lengths(value) and sequence IDs(key)
    max_len=max(ORF_lengths,key=ORF_lengths.get) #getting the ID of the Sequence read containing the maximum length ORF
    #printing the result
    print("Longest ORF in Frame %d is of length %d and its ID is %s." % (frame+1,ORF_lengths[max_len][0],max_len))
    print("Starting position is: %d" %ORF_lengths[max_len][1])
    
def main():
    from Bio import SeqIO
    filename=input('Enter the filename:')
    records=list(SeqIO.parse(filename,"fasta")) #parsing the file inputs into a list
    no_of_records=len(records) #no of records in the fasta file
    print("Number of the records in the given file is %d" % no_of_records)
    record_dict={}
    #creating a dictionary with ID as key and length of the sequences as value
    for a in records:
        record_dict[a.id]=len(a.seq)
    #calling the maximum?/minimum function
    if input("go to maximum or minimum?") =='yes':
        maxmin(record_dict)
    #calling the ORF function
    frame=0
    while frame<3: #looping three times for each frame
        finding_ORF(frame,records)
        frame+=1


