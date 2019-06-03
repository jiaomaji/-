import io
import re
import math
import pandas
import os
import time 
start = time.clock()
def Exchange_P(m):
    for k in range(-2,6):
        for A in ("a","g","c","t"):
            sum=0
            for B in ("a","g","c","t"):
                sum+=m[A+B][str(k)]
            for B in ("a","g","c","t"):
                m[A+B][str(k)]/=sum#将统计总数转换成概率
def Exchange_p(m):
    
    for A in ("a","g","c","t"):
        sum=0       
        for B in ("a","g","c","t"):
            sum+=m[A+B]
        for B in ("a","g","c","t"):
            m[A+B]/=sum#转换
def pp(m):
    print("****   ",end="")
    for A in ("a","g","c","t"):
        for B in ("a","g","c","t"):
            print(A+B+"     ",end="")
    print("\n")
    print("  ",end="")
    for A in ("a","g","c","t"):
        for B in ("a","g","c","t"):
            print('%.3f'%m[A+B],end="")
            print("   ",end="")        #制表
def PP_1(m):
     print("****   ",end="")
     for A in ("a","g","c","t"):
            for B in ("a","g","c","t"):
                print(A+B+"      ",end="")
     print("\n")
     for i in range(-2,6):
        print(i,end="")
        print("  ",end="")
        for A in ("a","g","c","t"):
            for B in ("a","g","c","t"):
                print('%.3f'%m[A+B][str(i)],end="")
                print("   ",end="")
        print("\n")      #制表
def deal_(file):
    print(file)
    os.system("cls")
    file.readline()
    line=file.readline().replace("..",",")
    num=len(line)
    line=line[8:(num-2)]
    donor=line.split(",")
    DNA_line=""
    for line in file.readlines():
        DNA_line+=line.strip()
    n=len(DNA_line)
    j=0
    for i in donor:
        if (j%2):
            for k in range(int(i)-2,int(i)+6):
                if k<n:
                    matrix_1[DNA_line[k]+DNA_line[k-1]][str(k-int(i))]+=1;
        j+=1
    for i in range(1,len(DNA_line)):
        if DNA_line[i] in ("a","g","c","t")and DNA_line[i-1] in ("a","g","c","t"):
            t=DNA_line[i]+DNA_line[i-1]
            matrix_2[t]+=1 
    file.close()  #求概率
def deal_training_set(fp):
    global confirm__#匹配的donor位点总数
    global new__#新的donor
    global lost__
    global pa
    global pg
    global pc
    global pt
    global First_p
    fp.readline()
    line=fp.readline().replace("..",",")
    num=len(line)
    line=line[8:(num-2)]
    donor=line.split(",")
    new_donor=[]
    DNA_line=""
    for line in fp.readlines():
        DNA_line+=line.strip()
    for i in range(3,len(DNA_line)-5):#打分
        if DNA_line[i-3] in ("a","g","c","t"): 
            sum=First_p[DNA_line[i-3]]
            for var in range(i-2,i+6):
                A=DNA_line[var]
                B=DNA_line[var-1]
                if A in ("a","g","c","t"):
                    if B in("a","g","c","t"): 
                        print(matrix_1[A+B][str(var-i)])
                        print(matrix_2[A+B])
                        print(matrix_1[A+B][str(var-i)]/matrix_2[A+B])
                        if not matrix_1[A+B][str(var-i)] or not matrix_2[A+B]:
                            continue
                        else:
                            sum+=math.log(matrix_1[A+B][str(var-i)]/matrix_2[A+B])
                    else:
                        i=var+3
                        break
                else:
                    i=var+4
                    break
            if sum>c:
                new_donor.append(str(i))
        confirm=0
        for var in new_donor:
            if str(var) in donor:
                confirm+=1
        confirm__+=confirm
        new__+=(len(donor)-confirm)
        lost__+=(len(new_donor)-confirm) #打分     
def Training_set_read():
    s = os.sep
    root = "G:" + s + "AA" + s 
    for rt, dirs, files in os.walk(root):
        for f in files:
            file=open(root+f)
            deal_(file) 
def score_training_set():
    s = os.sep
    root = "G:" + s + "AA" + s 
    for rt, dirs, files in os.walk(root):
        for f in files:
            file=open(root+f)
            deal_training_set(file)   #打分
            file.close()

def vv(var,n):
    global ssum
    if (n==9):
        f_ans.writelines(var)
        f_ans.write(" ")
        a=pre_mark(var)
        f_ans.write(a)
        if float(a)>=5:
            f_ans.write("*************************************")
        f_ans.write(" \n")
    else:
        for i in ("a","g","c","t"):
            vv(var+i,n+1)
    ssum+=1#建立一个数据库
def pre_mark(line):

    sum=First_p[line[0]]
    for i in range(1,8):
        if not matrix_1[line[i]+line[i-1]][str(i-3)] or not matrix_2[line[i]+line[i-1]]:
            continue
        else :
            sum+=math.log(matrix_1[line[i]+line[i-1]][str(i-3)]/matrix_2[line[i]+line[i-1]])
    return (str(sum))#给这个库汇中的所有可能打一个分


#正文
matrix_1={"aa":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"ag":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"ac":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"at":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"ga":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"gg":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"gc":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"gt":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"ca":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"cg":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"cc":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"ct":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"ta":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"tt":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"tc":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0},"tg":{"-3":0,"-2":0,"-1":0,"0":0,"1":0,"2":0,"3":0,"4":0,"5":0}}
matrix_2={"aa":0,"ag":0,"ac":0,"at":0,"tt":0,"ta":0,"tg":0,"tc":0,"ca":0,"cg":0,"cc":0,"ct":0,"ga":0,"gg":0,"gc":0,"gt":0}
Training_set_read()
confirm__=0 #匹配的donor位点总数
new__=0 #新的donor
lost__=0 #丢失的donor位点  
pa=matrix_1["aa"]["-2"]+matrix_1["ta"]["-2"]+matrix_1["ca"]["-2"]+matrix_1["ga"]["-2"]#pa为在一位置的前景概率
pt=matrix_1["at"]["-2"]+matrix_1["gt"]["-2"]+matrix_1["ct"]["-2"]+matrix_1["tt"]["-2"]
pc=matrix_1["ac"]["-2"]+matrix_1["gc"]["-2"]+matrix_1["cc"]["-2"]+matrix_1["tc"]["-2"]
pg=matrix_1["ag"]["-2"]+matrix_1["gg"]["-2"]+matrix_1["cg"]["-2"]+matrix_1["tg"]["-2"]
ba=matrix_2["aa"]+matrix_2["ta"]+matrix_2["ca"]+matrix_2["ga"]#ba为在一位置的背景概率
bt=matrix_2["at"]+matrix_2["tt"]+matrix_2["ct"]+matrix_2["gt"]
bc=matrix_2["ac"]+matrix_2["tc"]+matrix_2["cc"]+matrix_2["gc"]
bg=matrix_2["ag"]+matrix_2["tg"]+matrix_2["cg"]+matrix_2["gg"]
First_p={"A":0,"G":0,"C":0,"T":0}   #S(X)的前半部分
First_p["a"]=math.log(pa/ba)
First_p["g"]=math.log(pg/bg)
First_p["c"]=math.log(pc/bc)
First_p["t"]=math.log(pt/bt)
ssum=0
c=5
Exchange_P(matrix_1)#统计矩转换概率
Exchange_p(matrix_2)#统计矩阵转换概率
PP_1(matrix_1)#概率矩阵打印
print("\n\*******************************************************************************************")
pp(matrix_2)#打印概率矩阵
#score_training_set()
#print(confirm__)
#print(new__)
#print(lost__)
f_ans=open("G://DD//A.TXT","w")
var=""
vv(var,0)
f_ans.write(str(ssum))
f_ans.close()






























elapsed = (time.clock() - start)
print("Time used:",elapsed)
